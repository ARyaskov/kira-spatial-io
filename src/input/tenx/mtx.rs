//! MatrixMarket ingestion to canonical CSR with two passes over the input.

use std::io::BufRead;
use std::mem::size_of;
use std::path::Path;

use crate::config::{DuplicatePolicy, LoadConfig};
use crate::determinism::float::ensure_f32_finite_nonneg;
use crate::error::{IoPathExt, SpatialIoError};
use crate::input::util::open_text_maybe_gz;
use crate::model::csr::{BinsCsr, Indptr};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct MtxStream {
    pub n_rows: u32,
    pub n_cols: u32,
    pub nnz: u64,
}

#[derive(Clone, Copy, Debug)]
struct MtxTriplet {
    row0: u32,
    col0: u32,
    val: f32,
}

pub(crate) fn ensure_budget(
    bytes_needed: u64,
    cfg: &LoadConfig,
    context: &str,
) -> Result<(), SpatialIoError> {
    let budget_bytes = (cfg.memory_budget_mb as u64)
        .checked_mul(1024)
        .and_then(|v| v.checked_mul(1024))
        .ok_or_else(|| {
            SpatialIoError::MemoryLimitExceeded("memory budget conversion overflow".to_string())
        })?;

    if bytes_needed > budget_bytes {
        return Err(SpatialIoError::MemoryLimitExceeded(format!(
            "{context}: required {bytes_needed} bytes > budget {budget_bytes} bytes"
        )));
    }

    Ok(())
}

pub(crate) fn build_csr_from_mtx_path(
    matrix_path: &Path,
    cfg: &LoadConfig,
    old_to_new_bins: &[u32],
    feat_old_to_new: &[u32],
) -> Result<BinsCsr, SpatialIoError> {
    let n_bins = old_to_new_bins.len();
    let n_genes = feat_old_to_new.len();

    let mut counts = vec![0_u32; n_bins];
    let header_a = {
        let reader = open_text_maybe_gz(matrix_path)?;
        stream_mtx_triplets(reader, matrix_path, |triplet| {
            let barcode_old = triplet.col0 as usize;
            let bin_new = *old_to_new_bins.get(barcode_old).ok_or_else(|| {
                SpatialIoError::DimensionMismatch(format!(
                    "barcode index out of range for mapping: {barcode_old}"
                ))
            })? as usize;

            let gene_old = triplet.row0 as usize;
            let gene_new = *feat_old_to_new.get(gene_old).ok_or_else(|| {
                SpatialIoError::DimensionMismatch(format!(
                    "feature index out of range for mapping: {gene_old}"
                ))
            })?;
            if gene_new as usize >= n_genes {
                return Err(SpatialIoError::InvalidCsr(format!(
                    "gene_id out of range after remap: {gene_new}"
                )));
            }

            counts[bin_new] = counts[bin_new]
                .checked_add(1)
                .ok_or_else(|| {
                    SpatialIoError::InvalidCsr("row nnz count overflow while building CSR".to_string())
                })?;
            Ok(())
        })?
    };
    validate_header_dims(header_a, n_genes, n_bins)?;

    let final_bytes_upper = csr_bytes(n_bins, header_a.nnz)?;
    ensure_budget(
        final_bytes_upper,
        cfg,
        "final CSR would exceed budget (upper bound)",
    )?;

    let mut indptr64 = Vec::with_capacity(n_bins + 1);
    indptr64.push(0_u64);
    for &count in &counts {
        let next = indptr64
            .last()
            .copied()
            .and_then(|v| v.checked_add(count as u64))
            .ok_or_else(|| SpatialIoError::InvalidCsr("indptr overflow".to_string()))?;
        indptr64.push(next);
    }
    let nnz_upper = *indptr64
        .last()
        .ok_or_else(|| SpatialIoError::InvalidCsr("indptr construction failed".to_string()))?;

    let write_ptr_bytes = vec_bytes(n_bins, size_of::<u64>())?;
    ensure_peak_within_115(
        final_bytes_upper,
        write_ptr_bytes,
        cfg,
        "pass B temporary buffers",
    )?;

    let nnz_upper_usize = usize_from_u64(nnz_upper, "nnz does not fit usize")?;
    let mut indices = vec![0_u32; nnz_upper_usize];
    let mut data = vec![0_f32; nnz_upper_usize];
    let mut write_ptr = indptr64[..n_bins].to_vec();
    let mut insertion_order = vec![0_u32; nnz_upper_usize];
    let mut next_seq: u32 = 0;

    let reader = open_text_maybe_gz(matrix_path)?;
    let header_b = stream_mtx_triplets(reader, matrix_path, |triplet| {
        let bin_new = old_to_new_bins[triplet.col0 as usize] as usize;
        let gene_new = feat_old_to_new[triplet.row0 as usize];

        let pos_u64 = write_ptr[bin_new];
        write_ptr[bin_new] = write_ptr[bin_new]
            .checked_add(1)
            .ok_or_else(|| SpatialIoError::InvalidCsr("write pointer overflow".to_string()))?;

        let pos = usize_from_u64(pos_u64, "CSR index position overflow")?;
        if pos >= indices.len() {
            return Err(SpatialIoError::InvalidCsr(
                "write position out of bounds".to_string(),
            ));
        }

        indices[pos] = gene_new;
        data[pos] = triplet.val;
        insertion_order[pos] = next_seq;
        next_seq = next_seq.wrapping_add(1);
        Ok(())
    })?;
    if header_a != header_b {
        return Err(SpatialIoError::UnsupportedFormat(
            "matrix.mtx header changed between passes".to_string(),
        ));
    }

    for i in 0..n_bins {
        if write_ptr[i] != indptr64[i + 1] {
            return Err(SpatialIoError::InvalidCsr(format!(
                "row write count mismatch for row {i}: {} != {}",
                write_ptr[i],
                indptr64[i + 1]
            )));
        }
    }

    let max_row_nnz = counts.iter().copied().max().unwrap_or(0) as usize;
    let row_pairs_tmp = vec_bytes(max_row_nnz, size_of::<(u32, u32, f32)>())?;
    let indptr_tmp = vec_bytes(n_bins + 1, size_of::<u64>())?;
    ensure_peak_within_115(
        final_bytes_upper,
        row_pairs_tmp
            .checked_add(indptr_tmp)
            .ok_or_else(|| SpatialIoError::InvalidCsr("temporary size overflow".to_string()))?,
        cfg,
        "pass C temporary buffers",
    )?;

    let mut indptr2 = Vec::with_capacity(n_bins + 1);
    indptr2.push(0_u64);
    let mut write_cursor: usize = 0;

    for row in 0..n_bins {
        let start = usize_from_u64(indptr64[row], "row start offset overflow")?;
        let end = usize_from_u64(indptr64[row + 1], "row end offset overflow")?;

        let mut pairs: Vec<(u32, u32, f32)> = (start..end)
            .map(|i| (indices[i], insertion_order[i], data[i]))
            .collect();
        pairs.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

        let mut idx = 0;
        while idx < pairs.len() {
            let gene_id = pairs[idx].0;
            if gene_id as usize >= n_genes {
                return Err(SpatialIoError::InvalidCsr(format!(
                    "gene_id out of bounds after sort: {gene_id}"
                )));
            }

            let mut sum = pairs[idx].2;
            let mut dup_count = 1;
            idx += 1;
            while idx < pairs.len() && pairs[idx].0 == gene_id {
                sum += pairs[idx].2;
                dup_count += 1;
                idx += 1;
            }
            if dup_count > 1 && matches!(cfg.duplicate_policy, DuplicatePolicy::Error) {
                return Err(SpatialIoError::InvalidCsr(format!(
                    "duplicate (bin,gene)={row},{gene_id} with policy=Error"
                )));
            }
            ensure_f32_finite_nonneg(sum)?;

            indices[write_cursor] = gene_id;
            data[write_cursor] = sum;
            write_cursor += 1;
        }

        indptr2.push(write_cursor as u64);
    }

    indices.truncate(write_cursor);
    data.truncate(write_cursor);

    let csr = BinsCsr {
        indptr: Indptr::from_u64(indptr2),
        indices,
        data,
        n_bins: n_bins as u32,
        n_genes: n_genes as u32,
        nnz: write_cursor as u64,
        normalized: false,
    };

    validate_csr_invariants(&csr)?;
    Ok(csr)
}

fn validate_header_dims(
    header: MtxStream,
    n_rows: usize,
    n_cols: usize,
) -> Result<(), SpatialIoError> {
    if header.n_rows as usize != n_rows {
        return Err(SpatialIoError::DimensionMismatch(format!(
            "MTX n_rows {} != feature rows {}",
            header.n_rows, n_rows
        )));
    }
    if header.n_cols as usize != n_cols {
        return Err(SpatialIoError::DimensionMismatch(format!(
            "MTX n_cols {} != barcode count {}",
            header.n_cols, n_cols
        )));
    }
    Ok(())
}

fn stream_mtx_triplets<R, F>(
    mut r: R,
    path: &Path,
    mut on_triplet: F,
) -> Result<MtxStream, SpatialIoError>
where
    R: BufRead,
    F: FnMut(MtxTriplet) -> Result<(), SpatialIoError>,
{
    let mut line = String::new();

    let header = loop {
        line.clear();
        let n = r.read_line(&mut line).io_path(path)?;
        if n == 0 {
            return Err(SpatialIoError::UnsupportedFormat(
                "matrix.mtx missing dimension header".to_string(),
            ));
        }
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('%') {
            continue;
        }

        let mut parts = trimmed.split_ascii_whitespace();
        let n_rows = parse_u32(parts.next().unwrap_or_default(), "n_rows")?;
        let n_cols = parse_u32(parts.next().unwrap_or_default(), "n_cols")?;
        let nnz = parse_u64(parts.next().unwrap_or_default(), "nnz")?;
        if parts.next().is_some() {
            return Err(SpatialIoError::UnsupportedFormat(
                "invalid MatrixMarket dimension line".to_string(),
            ));
        }
        break MtxStream {
            n_rows,
            n_cols,
            nnz,
        };
    };

    let mut seen_triplets: u64 = 0;
    loop {
        line.clear();
        let n = r.read_line(&mut line).io_path(path)?;
        if n == 0 {
            break;
        }

        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('%') {
            continue;
        }

        let mut parts = trimmed.split_ascii_whitespace();
        let row1 = parse_u32(parts.next().unwrap_or_default(), "row")?;
        let col1 = parse_u32(parts.next().unwrap_or_default(), "col")?;
        let val = parse_f32(parts.next().unwrap_or_default(), "value")?;
        ensure_f32_finite_nonneg(val)?;

        if row1 == 0 || col1 == 0 {
            return Err(SpatialIoError::InvalidCsr(
                "MatrixMarket indices are 1-based; found 0".to_string(),
            ));
        }

        let row0 = row1 - 1;
        let col0 = col1 - 1;
        if row0 >= header.n_rows || col0 >= header.n_cols {
            return Err(SpatialIoError::DimensionMismatch(format!(
                "triplet index out of range: row {row1} col {col1}"
            )));
        }

        on_triplet(MtxTriplet { row0, col0, val })?;
        seen_triplets = seen_triplets.checked_add(1).ok_or_else(|| {
            SpatialIoError::InvalidCsr("nnz counter overflow".to_string())
        })?;
    }

    if seen_triplets != header.nnz {
        return Err(SpatialIoError::InvalidCsr(format!(
            "nnz mismatch: header {} but parsed {seen_triplets}",
            header.nnz
        )));
    }

    Ok(header)
}

fn validate_csr_invariants(csr: &BinsCsr) -> Result<(), SpatialIoError> {
    let expected_indptr_len = csr
        .n_bins
        .checked_add(1)
        .ok_or_else(|| SpatialIoError::InvalidCsr("n_bins overflow".to_string()))?
        as usize;
    if csr.indptr.len() != expected_indptr_len {
        return Err(SpatialIoError::InvalidCsr(format!(
            "indptr length {} != n_bins + 1 {expected_indptr_len}",
            csr.indptr.len(),
        )));
    }

    if csr.indptr.first().unwrap_or(1) != 0 {
        return Err(SpatialIoError::InvalidCsr(
            "indptr[0] must be 0".to_string(),
        ));
    }

    let nnz = csr.nnz;
    if csr.indptr.last().unwrap_or(0) != nnz {
        return Err(SpatialIoError::InvalidCsr(format!(
            "indptr[last] {} != nnz {nnz}",
            csr.indptr.last().unwrap_or(0)
        )));
    }

    if csr.indices.len() != nnz as usize || csr.data.len() != nnz as usize {
        return Err(SpatialIoError::InvalidCsr(
            "indices/data length mismatch with nnz".to_string(),
        ));
    }

    for row in 0..(csr.n_bins as usize) {
        let start = csr.indptr.get(row) as usize;
        let end = csr.indptr.get(row + 1) as usize;
        if start > end {
            return Err(SpatialIoError::InvalidCsr(format!(
                "indptr is not monotonic at row {row}"
            )));
        }

        let mut prev: Option<u32> = None;
        for idx in start..end {
            let gene = csr.indices[idx];
            if gene >= csr.n_genes {
                return Err(SpatialIoError::InvalidCsr(format!(
                    "gene index out of bounds in row {row}: {gene}"
                )));
            }
            if let Some(p) = prev
                && gene <= p
            {
                return Err(SpatialIoError::InvalidCsr(format!(
                    "row {row} indices must be strictly increasing"
                )));
            }
            prev = Some(gene);
            ensure_f32_finite_nonneg(csr.data[idx])?;
        }
    }

    Ok(())
}

fn ensure_peak_within_115(
    final_bytes: u64,
    temp_bytes: u64,
    cfg: &LoadConfig,
    context: &str,
) -> Result<(), SpatialIoError> {
    let peak = final_bytes
        .checked_add(temp_bytes)
        .ok_or_else(|| SpatialIoError::MemoryLimitExceeded("peak memory overflow".to_string()))?;
    if final_bytes >= 1024 * 1024 {
        let allowed_peak = final_bytes
            .checked_add(final_bytes / 100 * 15)
            .ok_or_else(|| {
                SpatialIoError::MemoryLimitExceeded("allowed peak memory overflow".to_string())
            })?;

        if peak > allowed_peak {
            return Err(SpatialIoError::MemoryLimitExceeded(format!(
                "{context}: temporary buffers exceed 1.15x final CSR size"
            )));
        }
    }

    ensure_budget(peak, cfg, context)
}

fn csr_bytes(n_bins: usize, nnz: u64) -> Result<u64, SpatialIoError> {
    let nnz_usize = usize_from_u64(nnz, "nnz does not fit usize")?;

    let indptr = vec_bytes(n_bins + 1, size_of::<u64>())?;
    let indices = vec_bytes(nnz_usize, size_of::<u32>())?;
    let data = vec_bytes(nnz_usize, size_of::<f32>())?;

    indptr
        .checked_add(indices)
        .and_then(|v| v.checked_add(data))
        .ok_or_else(|| SpatialIoError::MemoryLimitExceeded("CSR size overflow".to_string()))
}

fn vec_bytes(len: usize, elem_size: usize) -> Result<u64, SpatialIoError> {
    (len as u64)
        .checked_mul(elem_size as u64)
        .ok_or_else(|| SpatialIoError::MemoryLimitExceeded("byte-size overflow".to_string()))
}

fn usize_from_u64(v: u64, msg: &str) -> Result<usize, SpatialIoError> {
    usize::try_from(v).map_err(|_| SpatialIoError::InvalidCsr(msg.to_string()))
}

fn parse_u32(s: &str, field: &str) -> Result<u32, SpatialIoError> {
    s.parse::<u32>().map_err(|_| {
        SpatialIoError::UnsupportedFormat(format!("invalid {field} in matrix.mtx: {s}"))
    })
}

fn parse_u64(s: &str, field: &str) -> Result<u64, SpatialIoError> {
    s.parse::<u64>().map_err(|_| {
        SpatialIoError::UnsupportedFormat(format!("invalid {field} in matrix.mtx: {s}"))
    })
}

fn parse_f32(s: &str, field: &str) -> Result<f32, SpatialIoError> {
    s.parse::<f32>().map_err(|_| {
        SpatialIoError::UnsupportedFormat(format!("invalid {field} in matrix.mtx: {s}"))
    })
}

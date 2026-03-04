use std::mem::size_of;

use hdf5::Dataset;

use crate::config::LoadConfig;
use crate::determinism::float::ensure_f32_finite_nonneg;
use crate::error::SpatialIoError;
use crate::input::tenx::mtx::ensure_budget;
use crate::model::csr::BinsCsr;

pub(crate) fn build_csr_from_h5(
    file: &hdf5::File,
    cfg: &LoadConfig,
    old_to_new_bins: &[u32],
    feat_old_to_new: &[u32],
) -> Result<BinsCsr, SpatialIoError> {
    let shape_ds = file.dataset("/matrix/shape").map_err(|_| {
        SpatialIoError::UnsupportedFormat("missing /matrix/shape dataset".to_string())
    })?;
    let shape = read_shape(&shape_ds)?;
    if shape.len() != 2 {
        return Err(SpatialIoError::UnsupportedFormat(
            "unexpected /matrix/shape rank".to_string(),
        ));
    }

    let n_features = shape[0] as usize;
    let n_barcodes = shape[1] as usize;

    if n_features != feat_old_to_new.len() {
        return Err(SpatialIoError::DimensionMismatch(format!(
            "H5 shape features {} != loaded features {}",
            n_features,
            feat_old_to_new.len()
        )));
    }
    if n_barcodes != old_to_new_bins.len() {
        return Err(SpatialIoError::DimensionMismatch(format!(
            "H5 shape barcodes {} != loaded barcodes {}",
            n_barcodes,
            old_to_new_bins.len()
        )));
    }

    let indptr_ds = file.dataset("/matrix/indptr").map_err(|_| {
        SpatialIoError::UnsupportedFormat("missing /matrix/indptr dataset".to_string())
    })?;
    let indices_ds = file.dataset("/matrix/indices").map_err(|_| {
        SpatialIoError::UnsupportedFormat("missing /matrix/indices dataset".to_string())
    })?;
    let data_ds = file.dataset("/matrix/data").map_err(|_| {
        SpatialIoError::UnsupportedFormat("missing /matrix/data dataset".to_string())
    })?;

    let indptr_old = read_u64_like(&indptr_ds, "/matrix/indptr")?;
    if indptr_old.len() != n_barcodes + 1 {
        return Err(SpatialIoError::UnsupportedFormat(format!(
            "unexpected H5 sparse layout: indptr len {} != n_barcodes+1 {}",
            indptr_old.len(),
            n_barcodes + 1
        )));
    }

    let nnz = *indptr_old
        .last()
        .ok_or_else(|| SpatialIoError::UnsupportedFormat("empty /matrix/indptr".to_string()))?
        as usize;

    let final_csr_bytes = ((n_barcodes + 1) as u64)
        .checked_mul(size_of::<u64>() as u64)
        .and_then(|v| v.checked_add((nnz as u64) * size_of::<u32>() as u64))
        .and_then(|v| v.checked_add((nnz as u64) * size_of::<f32>() as u64))
        .ok_or_else(|| {
            SpatialIoError::MemoryLimitExceeded("final CSR size overflow".to_string())
        })?;
    ensure_budget(final_csr_bytes, cfg, "final CSR would exceed budget")?;

    let mut counts_new = vec![0_u64; n_barcodes];
    for old_bin in 0..n_barcodes {
        let c = indptr_old[old_bin + 1]
            .checked_sub(indptr_old[old_bin])
            .ok_or_else(|| SpatialIoError::InvalidCsr("indptr is not monotonic".to_string()))?;
        counts_new[old_to_new_bins[old_bin] as usize] = c;
    }

    let mut indptr_new = Vec::with_capacity(n_barcodes + 1);
    indptr_new.push(0_u64);
    for c in &counts_new {
        let next = indptr_new
            .last()
            .copied()
            .and_then(|v| v.checked_add(*c))
            .ok_or_else(|| SpatialIoError::InvalidCsr("indptr overflow".to_string()))?;
        indptr_new.push(next);
    }

    let nnz_new = *indptr_new
        .last()
        .ok_or_else(|| SpatialIoError::InvalidCsr("new indptr is empty".to_string()))?
        as usize;
    if nnz_new != nnz {
        return Err(SpatialIoError::InvalidCsr(format!(
            "nnz mismatch after bin remap: {} != {}",
            nnz_new, nnz
        )));
    }

    let mut indices = vec![0_u32; nnz];
    let mut data = vec![0_f32; nnz];
    let mut write_ptr = indptr_new[..n_barcodes].to_vec();

    let budget_bytes = (cfg.memory_budget_mb as u64)
        .checked_mul(1024)
        .and_then(|v| v.checked_mul(1024))
        .ok_or_else(|| {
            SpatialIoError::MemoryLimitExceeded("memory budget conversion overflow".to_string())
        })?;
    let fixed_temp = (n_barcodes as u64) * size_of::<u64>() as u64;
    let free_for_chunk = budget_bytes
        .saturating_sub(final_csr_bytes)
        .saturating_sub(fixed_temp);
    let bytes_per_elem = (size_of::<u32>() + size_of::<f32>()) as u64;
    let chunk_elems = usize::try_from((free_for_chunk / bytes_per_elem).max(1)).unwrap_or(1);

    for old_bin in 0..n_barcodes {
        let start = indptr_old[old_bin] as usize;
        let end = indptr_old[old_bin + 1] as usize;
        let new_bin = old_to_new_bins[old_bin] as usize;

        let mut offset_in_row = 0usize;
        let mut cursor = start;
        while cursor < end {
            let chunk_end = (cursor + chunk_elems).min(end);

            let idx_chunk = read_u32_like_range(&indices_ds, cursor, chunk_end, "/matrix/indices")?;
            let data_chunk = read_f32_like_range(&data_ds, cursor, chunk_end, "/matrix/data")?;
            if idx_chunk.len() != data_chunk.len() {
                return Err(SpatialIoError::DimensionMismatch(
                    "indices/data chunk lengths mismatch".to_string(),
                ));
            }

            let out_base = write_ptr[new_bin] as usize + offset_in_row;
            for (i, (&gene_old_u32, &value)) in idx_chunk.iter().zip(data_chunk.iter()).enumerate()
            {
                let gene_old = gene_old_u32 as usize;
                let gene_new = *feat_old_to_new.get(gene_old).ok_or_else(|| {
                    SpatialIoError::DimensionMismatch(format!(
                        "feature index out of mapping range: {}",
                        gene_old
                    ))
                })?;
                ensure_f32_finite_nonneg(value)?;

                indices[out_base + i] = gene_new;
                data[out_base + i] = value;
            }

            offset_in_row += idx_chunk.len();
            cursor = chunk_end;
        }

        write_ptr[new_bin] = write_ptr[new_bin]
            .checked_add((end - start) as u64)
            .ok_or_else(|| SpatialIoError::InvalidCsr("write pointer overflow".to_string()))?;
    }

    for row in 0..n_barcodes {
        if write_ptr[row] != indptr_new[row + 1] {
            return Err(SpatialIoError::InvalidCsr(format!(
                "row write mismatch at {}: {} != {}",
                row,
                write_ptr[row],
                indptr_new[row + 1]
            )));
        }
    }

    let mut indptr2 = Vec::with_capacity(n_barcodes + 1);
    indptr2.push(0_u64);
    let mut write_cursor = 0usize;

    for row in 0..n_barcodes {
        let start = indptr_new[row] as usize;
        let end = indptr_new[row + 1] as usize;

        let mut pairs: Vec<(u32, f32)> = (start..end).map(|i| (indices[i], data[i])).collect();
        pairs.sort_unstable_by_key(|(g, _)| *g);

        let mut i = 0usize;
        while i < pairs.len() {
            let gene = pairs[i].0;
            if gene as usize >= n_features {
                return Err(SpatialIoError::InvalidCsr(format!(
                    "gene index out of range after remap: {}",
                    gene
                )));
            }

            let mut sum = pairs[i].1;
            i += 1;
            while i < pairs.len() && pairs[i].0 == gene {
                sum += pairs[i].1;
                i += 1;
            }
            ensure_f32_finite_nonneg(sum)?;

            indices[write_cursor] = gene;
            data[write_cursor] = sum;
            write_cursor += 1;
        }

        indptr2.push(write_cursor as u64);
    }

    indices.truncate(write_cursor);
    data.truncate(write_cursor);

    validate_csr(BinsCsr {
        indptr: indptr2,
        indices,
        data,
        n_bins: n_barcodes as u32,
        n_genes: n_features as u32,
        nnz: write_cursor as u64,
        normalized: false,
    })
}

fn validate_csr(csr: BinsCsr) -> Result<BinsCsr, SpatialIoError> {
    if csr.indptr.len() != csr.n_bins as usize + 1 {
        return Err(SpatialIoError::InvalidCsr(
            "indptr length mismatch".to_string(),
        ));
    }
    if csr.indptr.first().copied().unwrap_or(1) != 0 {
        return Err(SpatialIoError::InvalidCsr(
            "indptr[0] must be 0".to_string(),
        ));
    }
    if csr.indptr.last().copied().unwrap_or(0) != csr.nnz {
        return Err(SpatialIoError::InvalidCsr(
            "indptr[last] must equal nnz".to_string(),
        ));
    }

    for row in 0..csr.n_bins as usize {
        let s = csr.indptr[row] as usize;
        let e = csr.indptr[row + 1] as usize;
        if s > e {
            return Err(SpatialIoError::InvalidCsr(format!(
                "indptr not monotonic at row {}",
                row
            )));
        }

        let mut prev: Option<u32> = None;
        for idx in s..e {
            let g = csr.indices[idx];
            if g >= csr.n_genes {
                return Err(SpatialIoError::InvalidCsr(format!(
                    "gene idx out of range {}",
                    g
                )));
            }
            if let Some(p) = prev
                && g <= p
            {
                return Err(SpatialIoError::InvalidCsr(format!(
                    "row {} indices are not strictly increasing",
                    row
                )));
            }
            prev = Some(g);
            ensure_f32_finite_nonneg(csr.data[idx])?;
        }
    }

    Ok(csr)
}

fn read_shape(ds: &Dataset) -> Result<Vec<u64>, SpatialIoError> {
    if let Ok(v) = ds.read_raw::<u64>() {
        return Ok(v);
    }
    if let Ok(v) = ds.read_raw::<u32>() {
        return Ok(v.into_iter().map(|x| x as u64).collect());
    }
    if let Ok(v) = ds.read_raw::<i64>() {
        return v
            .into_iter()
            .map(|x| {
                if x < 0 {
                    Err(SpatialIoError::UnsupportedFormat(
                        "negative value in /matrix/shape".to_string(),
                    ))
                } else {
                    Ok(x as u64)
                }
            })
            .collect();
    }
    if let Ok(v) = ds.read_raw::<i32>() {
        return v
            .into_iter()
            .map(|x| {
                if x < 0 {
                    Err(SpatialIoError::UnsupportedFormat(
                        "negative value in /matrix/shape".to_string(),
                    ))
                } else {
                    Ok(x as u64)
                }
            })
            .collect();
    }

    Err(SpatialIoError::UnsupportedFormat(
        "unsupported /matrix/shape dtype".to_string(),
    ))
}

fn read_u64_like(ds: &Dataset, path: &str) -> Result<Vec<u64>, SpatialIoError> {
    if let Ok(v) = ds.read_raw::<u64>() {
        return Ok(v);
    }
    if let Ok(v) = ds.read_raw::<u32>() {
        return Ok(v.into_iter().map(|x| x as u64).collect());
    }
    if let Ok(v) = ds.read_raw::<i64>() {
        return v
            .into_iter()
            .map(|x| {
                if x < 0 {
                    Err(SpatialIoError::UnsupportedFormat(format!(
                        "negative value in {}",
                        path
                    )))
                } else {
                    Ok(x as u64)
                }
            })
            .collect();
    }
    if let Ok(v) = ds.read_raw::<i32>() {
        return v
            .into_iter()
            .map(|x| {
                if x < 0 {
                    Err(SpatialIoError::UnsupportedFormat(format!(
                        "negative value in {}",
                        path
                    )))
                } else {
                    Ok(x as u64)
                }
            })
            .collect();
    }

    Err(SpatialIoError::UnsupportedFormat(format!(
        "unsupported dtype in {}",
        path
    )))
}

fn read_u32_like_range(
    ds: &Dataset,
    start: usize,
    end: usize,
    path: &str,
) -> Result<Vec<u32>, SpatialIoError> {
    if let Ok(v) = ds.read_raw::<u32>() {
        return Ok(slice_range(&v, start, end, path)?.to_vec());
    }
    if let Ok(v) = ds.read_raw::<u64>() {
        let mut out = Vec::with_capacity(end.saturating_sub(start));
        for &x in slice_range(&v, start, end, path)? {
            if x > u32::MAX as u64 {
                return Err(SpatialIoError::UnsupportedFormat(format!(
                    "value too large for u32 in {}",
                    path
                )));
            }
            out.push(x as u32);
        }
        return Ok(out);
    }
    if let Ok(v) = ds.read_raw::<i32>() {
        let mut out = Vec::with_capacity(end.saturating_sub(start));
        for &x in slice_range(&v, start, end, path)? {
            if x < 0 {
                return Err(SpatialIoError::UnsupportedFormat(format!(
                    "negative value in {}",
                    path
                )));
            }
            out.push(x as u32);
        }
        return Ok(out);
    }

    Err(SpatialIoError::UnsupportedFormat(format!(
        "unsupported dtype in {}",
        path
    )))
}

fn read_f32_like_range(
    ds: &Dataset,
    start: usize,
    end: usize,
    path: &str,
) -> Result<Vec<f32>, SpatialIoError> {
    if let Ok(v) = ds.read_raw::<f32>() {
        return Ok(slice_range(&v, start, end, path)?.to_vec());
    }
    if let Ok(v) = ds.read_raw::<f64>() {
        let mut out = Vec::with_capacity(end.saturating_sub(start));
        for &x in slice_range(&v, start, end, path)? {
            let f = x as f32;
            ensure_f32_finite_nonneg(f)?;
            out.push(f);
        }
        return Ok(out);
    }
    if let Ok(v) = ds.read_raw::<i32>() {
        let mut out = Vec::with_capacity(end.saturating_sub(start));
        for &x in slice_range(&v, start, end, path)? {
            let f = x as f32;
            ensure_f32_finite_nonneg(f)?;
            out.push(f);
        }
        return Ok(out);
    }

    Err(SpatialIoError::UnsupportedFormat(format!(
        "unsupported dtype in {}",
        path
    )))
}

fn slice_range<'a, T>(
    data: &'a [T],
    start: usize,
    end: usize,
    path: &str,
) -> Result<&'a [T], SpatialIoError> {
    data.get(start..end).ok_or_else(|| {
        SpatialIoError::DimensionMismatch(format!(
            "slice out of bounds in {}: {}..{} with len {}",
            path,
            start,
            end,
            data.len()
        ))
    })
}

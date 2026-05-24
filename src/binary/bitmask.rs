//! Fixed-width `u64` LE serialization for `BitVec`.

use bitvec::order::Lsb0;
use bitvec::vec::BitVec;

/// Packs bits into LSB0-ordered little-endian `u64` words; trailing bits in the last word are zero.
pub fn tissue_mask_to_u64_words(bits: &BitVec) -> Vec<u64> {
    let n = bits.len();
    let n_words = n.div_ceil(64);
    let mut words = vec![0_u64; n_words];
    for (i, b) in bits.iter().enumerate() {
        if *b {
            words[i / 64] |= 1_u64 << (i % 64);
        }
    }
    words
}

/// Reconstructs a `BitVec` from platform-independent `u64` words and an exact bit count.
pub fn tissue_mask_from_u64_words(words: &[u64], n_bits: usize) -> BitVec {
    let mut out: BitVec<usize, Lsb0> = BitVec::with_capacity(n_bits);
    out.resize(n_bits, false);
    for i in 0..n_bits {
        let w = words[i / 64];
        if ((w >> (i % 64)) & 1) == 1 {
            out.set(i, true);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use bitvec::prelude::*;

    #[test]
    fn roundtrip_empty() {
        let bits: BitVec = bitvec![];
        let words = tissue_mask_to_u64_words(&bits);
        assert!(words.is_empty());
        let back = tissue_mask_from_u64_words(&words, 0);
        assert_eq!(back.len(), 0);
    }

    #[test]
    fn roundtrip_aligned() {
        let mut bits: BitVec = BitVec::repeat(false, 128);
        bits.set(0, true);
        bits.set(63, true);
        bits.set(64, true);
        bits.set(127, true);
        let words = tissue_mask_to_u64_words(&bits);
        assert_eq!(words.len(), 2);
        let back = tissue_mask_from_u64_words(&words, 128);
        assert_eq!(back, bits);
    }

    #[test]
    fn roundtrip_unaligned() {
        let mut bits: BitVec = BitVec::repeat(false, 130);
        bits.set(0, true);
        bits.set(65, true);
        bits.set(129, true);
        let words = tissue_mask_to_u64_words(&bits);
        assert_eq!(words.len(), 3);
        let back = tissue_mask_from_u64_words(&words, 130);
        assert_eq!(back, bits);
    }
}

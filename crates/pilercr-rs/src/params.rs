use crate::{dp::DpParams, kmer::FilterParams};

pub const DEFAULT_MEMORY_BUDGET: usize = 512 * 1024 * 1024;

pub fn estimated_memory(sequence_length: usize, params: &FilterParams) -> f64 {
    let index_entry_size = 8 + 2 * std::mem::size_of::<usize>();
    let width = params.tube_offset as f64 + params.seed_diffs as f64;
    let tubes = ((sequence_length as f64 + width - 1.0) / params.tube_offset as f64).floor() + 1.0;
    2.0 * index_entry_size as f64 * 4f64.powi(params.word_size as i32)
        + tubes * 12.0
        + sequence_length as f64
}

pub fn get_params(
    sequence_length: usize,
    diameter: usize,
    min_hit_length: usize,
    min_identity: f64,
    memory_budget: usize,
) -> Result<(FilterParams, DpParams), String> {
    let bound = (i32::MAX / 8) as usize;
    if sequence_length > bound || diameter == 0 || diameter > bound {
        return Err("sequence length or diameter exceeds supported coordinates".into());
    }
    if min_hit_length <= 4 || min_hit_length > bound {
        return Err("minimum hit length must be greater than 4 and fit DP coordinates".into());
    }
    if !min_identity.is_finite() || !(0.0..=1.0).contains(&min_identity) {
        return Err("minimum identity must be finite and between 0 and 1".into());
    }
    if memory_budget == 0 {
        return Err("memory budget must be positive (bytes)".into());
    }
    let minimum_word = ((diameter as f64).log(4.0) - 10f64.log(4.0) + 0.5) as isize;
    let mut seed_length = min_hit_length as isize;
    let mut seed_diffs = (min_hit_length as f64 * (1.0 - min_identity) + 0.5) as usize;
    loop {
        for word_size in (minimum_word.max(1) as usize..=13).rev() {
            if seed_length + 1 - word_size as isize * (seed_diffs as isize + 1) <= 0 {
                continue;
            }
            let filter = FilterParams {
                word_size,
                seed_length: seed_length as usize,
                seed_diffs,
                tube_offset: seed_diffs + 32,
                diameter,
            };
            if estimated_memory(sequence_length, &filter) <= memory_budget as f64
                && diameter as f64 / 4f64.powi(word_size as i32) <= 10.0
            {
                return Ok((
                    filter,
                    DpParams {
                        min_hit_length,
                        min_identity,
                        max_hit_length: 192.max(min_hit_length),
                    },
                ));
            }
        }
        if seed_length >= (min_hit_length / 4) as isize {
            seed_length -= 4;
        } else if seed_diffs > 0 {
            seed_diffs -= 1;
        } else {
            return Err("failed to find filter parameters within the memory budget".into());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reference_defaults_and_memory_selection() {
        let (filter, dp) = get_params(1_000_000, 512, 16, 0.94, DEFAULT_MEMORY_BUDGET).unwrap();
        assert_eq!(
            (
                filter.word_size,
                filter.seed_length,
                filter.seed_diffs,
                filter.tube_offset
            ),
            (8, 16, 1, 33)
        );
        assert_eq!(dp, DpParams::default());
        let (filter, _) = get_params(1000, 512, 100, 1.0, DEFAULT_MEMORY_BUDGET).unwrap();
        assert_eq!(filter.word_size, 11);
        assert!(estimated_memory(1000, &filter) <= DEFAULT_MEMORY_BUDGET as f64);
    }

    #[test]
    fn rejects_invalid_and_impossible_parameters() {
        for (length, identity, memory) in [
            (4, 0.94, 1000),
            (16, f64::NAN, 1000),
            (16, 1.01, 1000),
            (16, 0.94, 0),
            (16, 0.94, 1),
        ] {
            assert!(get_params(1000, 512, length, identity, memory).is_err());
        }
        assert!(get_params(0, 1, 16, 0.94, DEFAULT_MEMORY_BUDGET).is_ok());
    }
}

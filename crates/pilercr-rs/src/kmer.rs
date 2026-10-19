use std::collections::{HashMap, VecDeque};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FilterParams {
    pub word_size: usize,
    pub seed_length: usize,
    pub seed_diffs: usize,
    pub tube_offset: usize,
    pub diameter: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FilterHit {
    pub q_from: usize,
    pub q_to: usize,
    pub diagonal: isize,
}

pub fn encode_word(word: &[u8]) -> Option<usize> {
    if word.len() > usize::BITS as usize / 2 {
        return None;
    }
    word.iter().try_fold(0, |code, base| {
        let letter = match base {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => return None,
        };
        Some((code << 2) | letter)
    })
}

#[derive(Default, Clone, Copy)]
struct TubeState {
    q_lo: usize,
    q_hi: usize,
    count: usize,
}

struct Tubes {
    states: Vec<TubeState>,
    hits: Vec<FilterHit>,
    length: usize,
    offset: usize,
    word_size: usize,
    min_count: usize,
    max_distance: usize,
}

impl Tubes {
    fn save(&mut self, index: usize) -> Result<(), String> {
        let state = self.states[index % self.states.len()];
        self.hits.try_reserve(1).map_err(|e| e.to_string())?;
        self.hits.push(FilterHit {
            q_from: state.q_lo,
            q_to: state.q_hi + self.word_size,
            diagonal: self.length as isize - (index * self.offset) as isize,
        });
        Ok(())
    }

    fn hit(&mut self, index: usize, q: usize) -> Result<(), String> {
        let slot = index % self.states.len();
        let state = self.states[slot];
        if state.count == 0 || (q > state.q_hi && q - state.q_hi > self.max_distance) {
            if state.count >= self.min_count {
                self.save(index)?;
            }
            self.states[slot] = TubeState {
                q_lo: q,
                q_hi: q,
                count: 1,
            };
        } else {
            self.states[slot].count = state
                .count
                .checked_add(1)
                .ok_or("tube hit count overflow")?;
            self.states[slot].q_hi = q;
        }
        Ok(())
    }

    fn end(&mut self, q: usize) -> Result<(), String> {
        if q == 0 {
            return Ok(());
        }
        let index = q / self.offset;
        let slot = index % self.states.len();
        if self.states[slot].count >= self.min_count {
            self.save(index)?;
        }
        self.states[slot].count = 0;
        Ok(())
    }
}

pub fn filter(sequence: &[u8], params: &FilterParams) -> Result<Vec<FilterHit>, String> {
    let FilterParams {
        word_size: k,
        seed_length,
        seed_diffs,
        tube_offset,
        diameter,
    } = *params;
    if k == 0 || k >= usize::BITS as usize / 2 {
        return Err("word_size must be positive and 4^word_size must fit usize".into());
    }
    if tube_offset == 0 || tube_offset < seed_diffs {
        return Err("tube_offset must be positive and at least seed_diffs".into());
    }
    if diameter < k {
        return Err("diameter must be at least word_size".into());
    }
    let min_count = seed_diffs
        .checked_add(1)
        .and_then(|errors| k.checked_mul(errors))
        .and_then(|covered| seed_length.checked_add(1)?.checked_sub(covered))
        .filter(|&count| count > 0)
        .ok_or("seed parameters must give a positive, representable word threshold")?;
    let width = tube_offset
        .checked_add(seed_diffs)
        .ok_or("tube width overflow")?;
    let n = sequence.len();
    let tube_count = n
        .checked_add(width - 1)
        .and_then(|span| (span / tube_offset).checked_add(1))
        .ok_or("tube count overflow")?;
    let max_hit_index = n.checked_mul(2).ok_or("diagonal index overflow")? / tube_offset;
    let max_index = max_hit_index.max(tube_count - 1);
    if [n, seed_length, seed_diffs, tube_offset, diameter, width]
        .iter()
        .any(|&value| value > isize::MAX as usize)
        || max_index
            .checked_mul(tube_offset)
            .is_none_or(|value| value > isize::MAX as usize)
    {
        return Err("parameters or diagonal exceed signed coordinate range".into());
    }
    let window = diameter - k + 1;
    if n < k || n - k < window || !sequence.windows(k).any(|word| encode_word(word).is_some()) {
        return Ok(Vec::new());
    }
    let mut states = Vec::new();
    states
        .try_reserve_exact(tube_count)
        .map_err(|e| e.to_string())?;
    states.resize(tube_count, TubeState::default());
    let mut tubes = Tubes {
        states,
        hits: Vec::new(),
        length: n,
        offset: tube_offset,
        word_size: k,
        min_count,
        max_distance: seed_length - k,
    };
    let mut buckets: HashMap<usize, VecDeque<usize>> = HashMap::new();
    let mut ticker = width;
    for end in 0..=n - k {
        if end >= window {
            let start = end - window;
            if let Some(code) = encode_word(&sequence[start..start + k]) {
                let bucket = buckets.get_mut(&code).ok_or("missing FIFO bucket")?;
                if bucket.pop_front() != Some(start) {
                    return Err("inconsistent FIFO bucket".into());
                }
                for &q in bucket.iter() {
                    let diagonal_index = n - start + q;
                    let index = diagonal_index / tube_offset;
                    tubes.hit(index, q)?;
                    if diagonal_index % tube_offset < seed_diffs {
                        let preceding = if index == 0 {
                            tube_count - 1
                        } else {
                            index - 1
                        };
                        tubes.hit(preceding, q)?;
                    }
                }
                if bucket.is_empty() {
                    buckets.remove(&code);
                }
            }
        }
        ticker -= 1;
        if ticker == 0 {
            tubes.end(end)?;
            ticker = tube_offset;
        }
        if let Some(code) = encode_word(&sequence[end..end + k]) {
            if !buckets.contains_key(&code) {
                buckets.try_reserve(1).map_err(|e| e.to_string())?;
            }
            let bucket = buckets.entry(code).or_default();
            bucket.try_reserve(1).map_err(|e| e.to_string())?;
            bucket.push_back(end);
        }
    }
    tubes.end(n - 1)?;
    for index in 0..tube_count {
        if tubes.states[index].count >= min_count {
            tubes.save(index)?;
            tubes.states[index].count = 0;
        }
    }
    Ok(tubes.hits)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params() -> FilterParams {
        FilterParams {
            word_size: 2,
            seed_length: 4,
            seed_diffs: 0,
            tube_offset: 2,
            diameter: 6,
        }
    }

    fn triples(hits: Vec<FilterHit>) -> Vec<(usize, usize, isize)> {
        hits.into_iter()
            .map(|h| (h.q_from, h.q_to, h.diagonal))
            .collect()
    }

    #[test]
    fn word_encoding() {
        assert_eq!(encode_word(b"AcgT"), Some(27));
        assert_eq!(encode_word(b""), Some(0));
        assert_eq!(encode_word(b"ACNT"), None);
        assert_eq!(encode_word(b"U"), None);
        assert_eq!(encode_word(&[255]), None);
        assert_eq!(
            encode_word(&vec![b'T'; usize::BITS as usize / 2]),
            Some(usize::MAX)
        );
        assert_eq!(encode_word(&vec![b'A'; usize::BITS as usize / 2 + 1]), None);
    }

    #[test]
    fn short_and_ambiguous_inputs() {
        for sequence in [b"".as_slice(), b"A", b"AAAAAA", b"NNNNNNNNNNNNNN"] {
            assert!(filter(sequence, &params()).unwrap().is_empty());
        }
        let mut p = params();
        p.diameter = 100;
        assert!(filter(b"ACGTACGTACGT", &p).unwrap().is_empty());
    }

    #[test]
    fn fifo_and_final_flush_order() {
        assert_eq!(
            triples(filter(b"AAAAAAAAAAAA", &params()).unwrap()),
            vec![(2, 10, 12), (4, 11, 10), (1, 8, 0),]
        );
    }

    #[test]
    fn overlapping_tubes_and_ambiguous_gaps() {
        let p = FilterParams {
            seed_length: 5,
            seed_diffs: 1,
            ..params()
        };
        for (sequence, expected) in [
            (
                b"AAAAAAAAAAAA".as_slice(),
                vec![(4, 11, 12), (1, 9, 0), (2, 11, -2)],
            ),
            (
                b"AaAaNaaAaAaAaAaAa",
                vec![(1, 4, 1), (6, 13, 1), (5, 16, 17), (1, 15, -1)],
            ),
            (
                b"AAAAAANNNNNNNNNNAAAAAAAAAAAA",
                vec![
                    (1, 6, 0),
                    (2, 6, -2),
                    (20, 27, 28),
                    (17, 25, 0),
                    (18, 27, -2),
                ],
            ),
        ] {
            assert_eq!(triples(filter(sequence, &p).unwrap()), expected);
        }
        let p = FilterParams {
            word_size: 1,
            seed_length: 2,
            seed_diffs: 1,
            tube_offset: 1,
            diameter: 3,
        };
        assert_eq!(
            triples(filter(b"AaNaaRYaaA", &p).unwrap()),
            vec![(1, 2, -1), (1, 2, 0), (3, 4, 10), (4, 5, 0), (3, 5, -1),]
        );
    }

    #[test]
    fn rejects_invalid_and_overflowing_parameters() {
        let bad = [
            FilterParams {
                word_size: 0,
                ..params()
            },
            FilterParams {
                word_size: usize::MAX,
                ..params()
            },
            FilterParams {
                tube_offset: 0,
                ..params()
            },
            FilterParams {
                seed_diffs: 3,
                ..params()
            },
            FilterParams {
                seed_length: 0,
                ..params()
            },
            FilterParams {
                diameter: 1,
                ..params()
            },
            FilterParams {
                seed_length: usize::MAX,
                ..params()
            },
            FilterParams {
                seed_diffs: usize::MAX,
                tube_offset: usize::MAX,
                ..params()
            },
            FilterParams {
                tube_offset: usize::MAX,
                ..params()
            },
            FilterParams {
                diameter: usize::MAX,
                ..params()
            },
        ];
        for p in bad {
            assert!(filter(b"", &p).is_err(), "{p:?}");
        }
    }
}

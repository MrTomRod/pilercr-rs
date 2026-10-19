use crate::filter::Trapezoid;
use crate::seq::Sequence;
use rayon::prelude::*;

#[derive(Clone, Debug, PartialEq)]
pub struct DpParams {
    pub min_hit_length: usize,
    pub min_identity: f64,
    pub max_hit_length: usize,
}

impl Default for DpParams {
    fn default() -> Self {
        Self {
            min_hit_length: 16,
            min_identity: 0.94,
            max_hit_length: 192,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct DpHit {
    pub a_lo: usize,
    pub a_hi: usize,
    pub b_lo: usize,
    pub b_hi: usize,
    pub low_diagonal: isize,
    pub high_diagonal: isize,
    pub score: i32,
    pub error: f32,
}

const MAX_GAP: isize = 5;
const DIFF_COST: i32 = 3;
const BLOCK_COST: i32 = 15;

struct Row {
    origin: isize,
    direction: isize,
    values: Vec<i32>,
}

impl Row {
    fn new(origin: isize, direction: isize) -> Self {
        Self {
            origin,
            direction,
            values: Vec::new(),
        }
    }

    fn get(&self, column: isize) -> i32 {
        self.values[((column - self.origin) * self.direction) as usize]
    }
}

struct End {
    a: isize,
    b: isize,
    low: isize,
    high: isize,
    score: i32,
}

fn forward(bases: &[u8], mid: isize, mut lo: isize, mut hi: isize, limit: isize) -> End {
    let n = bases.len() as isize;
    lo = lo.max(0);
    hi = hi.min(n);
    let mut row = Row::new(lo, 1);
    for _ in lo..=hi {
        row.values.push(0);
    }
    let initial_hi = hi;
    hi = (hi + MAX_GAP).min(n);
    for j in initial_hi + 1..=hi {
        row.values.push(row.get(j - 1) - DIFF_COST);
    }
    let mut end = End {
        a: lo,
        b: mid,
        low: mid - hi,
        high: mid - lo,
        score: 0,
    };
    for i in mid..(mid + limit - 1).min(n) {
        if lo > hi {
            break;
        }
        let mut next = Row::new(lo, 1);
        let mut v = row.get(lo);
        let mut c = v - DIFF_COST;
        next.values.push(c);
        for j in lo + 1..=hi {
            let t = c;
            c = v;
            v = row.get(j);
            if bases[i as usize] == bases[(j - 1) as usize] {
                c += 4;
            }
            c = c.max(v).max(t) - DIFF_COST;
            next.values.push(c);
            if c >= end.score {
                end.score = c;
                end.b = i + 1;
                end.a = j;
            }
        }
        let mut j = hi + 1;
        if j <= n {
            if bases[i as usize] == bases[(j - 1) as usize] {
                v += 4;
            }
            v = v.max(c) - DIFF_COST;
            next.values.push(v);
            if v > end.score {
                end.score = v;
                end.b = i + 1;
                end.a = j;
            }
            j += 1;
            while j <= n {
                v -= DIFF_COST;
                if v < end.score - BLOCK_COST {
                    break;
                }
                next.values.push(v);
                j += 1;
            }
        }
        hi = j - 1;
        while lo <= hi && next.get(lo) < end.score - BLOCK_COST {
            lo += 1;
        }
        while lo <= hi && next.get(hi) < end.score - BLOCK_COST {
            hi -= 1;
        }
        end.high = end.high.max(i + 1 - lo);
        end.low = end.low.min(i + 1 - hi);
        row = next;
    }
    end
}

fn reverse(
    bases: &[u8],
    top: isize,
    column: isize,
    bot: isize,
    mut drop: i32,
    limit: isize,
) -> End {
    let mut lo = column.max(0);
    let mut hi = column.min(bases.len() as isize);
    let mut row = Row::new(hi, -1);
    for _ in lo..=hi {
        row.values.push(0);
    }
    let initial_lo = lo;
    lo = (lo - MAX_GAP).max(0);
    for j in (lo..initial_lo).rev() {
        row.values.push(row.get(j + 1) - DIFF_COST);
    }
    let mut end = End {
        a: lo,
        b: top,
        low: top - hi,
        high: top - lo,
        score: 0,
    };
    if top - 1 <= bot {
        drop = BLOCK_COST;
    }
    for i in ((top - limit).max(0)..top).rev() {
        if lo > hi {
            break;
        }
        let mut next = Row::new(hi, -1);
        let mut v = row.get(hi);
        let mut c = v - DIFF_COST;
        next.values.push(c);
        for j in (lo..hi).rev() {
            let t = c;
            c = v;
            v = row.get(j);
            if bases[i as usize] == bases[j as usize] {
                c += 4;
            }
            c = c.max(v).max(t) - DIFF_COST;
            next.values.push(c);
            if c >= end.score {
                end.score = c;
                end.b = i;
                end.a = j;
            }
        }
        let mut j = lo - 1;
        if j >= 0 {
            if bases[i as usize] == bases[j as usize] {
                v += 4;
            }
            v = v.max(c) - DIFF_COST;
            next.values.push(v);
            if v > end.score {
                end.score = v;
                end.b = i;
                end.a = j;
            }
            j -= 1;
            while j >= 0 {
                v -= DIFF_COST;
                if v < end.score - drop {
                    break;
                }
                next.values.push(v);
                j -= 1;
            }
        }
        lo = j + 1;
        while lo <= hi && next.get(lo) < end.score - drop {
            lo += 1;
        }
        while lo <= hi && next.get(hi) < end.score - drop {
            hi -= 1;
        }
        if i == bot {
            drop = BLOCK_COST;
        }
        end.high = end.high.max(i - lo);
        end.low = end.low.min(i - hi);
        row = next;
    }
    end
}

pub fn sort_hits(hits: &mut [DpHit]) {
    crate::sort::reference_sort_by_key(hits, |hit| hit.a_lo);
}

fn same_contig(sequence: &Sequence, hit: &DpHit) -> bool {
    let Some((index, _)) = sequence.global_to_local(hit.a_lo) else {
        return false;
    };
    [hit.a_hi - 1, hit.b_lo, hit.b_hi - 1].iter().all(|&pos| {
        sequence
            .global_to_local(pos)
            .is_some_and(|(other, _)| other == index)
    })
}

struct Candidate {
    hits: Vec<DpHit>,
    covered: Vec<usize>,
}

/// Computes the independent alignment work for one trap in the sorted list.
///
/// Coverage marks only reference later traps; whether this trap is itself
/// covered is decided by the sequential pass in [`align_traps`].
fn align_candidate(
    sequence: &Sequence,
    sorted: &[Trapezoid],
    current: usize,
    params: &DpParams,
    min_len: isize,
    limit: isize,
) -> Candidate {
    let mut candidate = Candidate {
        hits: Vec::new(),
        covered: Vec::new(),
    };
    let trap = sorted[current];
    if trap.left < 8 || trap.bottom == trap.top {
        return candidate;
    }
    let n = sequence.bases.len() as isize;
    let mut pending = vec![(trap.bottom as isize, trap.top as isize)];
    while let Some((bottom, top)) = pending.pop() {
        let mid = (bottom + top) / 2;
        let lo = mid - trap.right;
        let hi = mid - trap.left;
        if hi < 0 || lo > n as isize {
            continue;
        }
        let lend = forward(&sequence.bases, mid, lo, hi, limit);
        let mut x = 1;
        let hend = loop {
            let end = reverse(
                &sequence.bases,
                lend.b,
                lend.a,
                mid + MAX_GAP,
                BLOCK_COST + 2 * x * DIFF_COST,
                limit,
            );
            if !(end.b > mid + x as isize * MAX_GAP && end.score < lend.score) {
                break end;
            }
            x += 1;
        };
        if lend.b - hend.b >= min_len && lend.a - hend.a >= min_len {
            let indel = ((hend.a - hend.b) - (lend.a - lend.b)).abs();
            let error = (0.25
                - (hend.score as f64 - indel as f64) / (4.0 * (lend.b - hend.b) as f64))
                as f32;
            if error as f64 <= 1.0 - params.min_identity {
                for j in current + 1..sorted.len() {
                    let t = sorted[j];
                    if t.bottom as isize >= lend.b {
                        break;
                    }
                    let width = t.right - t.left + 1;
                    let height = (t.top - t.bottom + 1) as isize;
                    let left = t.left.max(hend.low);
                    let right = t.right.min(hend.high);
                    if left > right {
                        continue;
                    }
                    let overlap = if t.top as isize > lend.b {
                        lend.b - t.bottom as isize + 1
                    } else {
                        height
                    };
                    if ((right - left + 1) as f64 / width as f64) * (overlap as f64 / height as f64)
                        > 0.99
                    {
                        candidate.covered.push(j);
                    }
                }
                let hit = DpHit {
                    a_lo: hend.a as usize,
                    a_hi: lend.a as usize,
                    b_lo: hend.b as usize,
                    b_hi: lend.b as usize,
                    low_diagonal: -hend.high,
                    high_diagonal: -hend.low,
                    score: hend.score,
                    error,
                };
                let reference_length = (hit.a_hi - hit.a_lo + hit.b_hi - hit.b_lo + 2) / 2;
                if (16..=128).contains(&reference_length) && same_contig(sequence, &hit) {
                    candidate.hits.push(hit);
                }
            }
        }
        let upper_bottom = lend.b + MAX_GAP;
        if top - upper_bottom > min_len {
            pending.push((upper_bottom, top));
        }
        let lower_top = hend.b - MAX_GAP;
        if lower_top - bottom > min_len && lower_top < top - MAX_GAP {
            pending.push((bottom, lower_top));
        }
    }
    candidate
}

pub fn align_traps(
    sequence: &Sequence,
    traps: &[Trapezoid],
    params: &DpParams,
) -> Result<Vec<DpHit>, String> {
    let bound = (i32::MAX / 8) as usize;
    if params.min_hit_length == 0
        || params.max_hit_length < params.min_hit_length
        || params.max_hit_length > bound
        || !params.min_identity.is_finite()
        || !(0.0..=1.0).contains(&params.min_identity)
    {
        return Err("invalid DP length or identity parameters".into());
    }
    let n = sequence.bases.len();
    if n > bound {
        return Err("sequence is too long for DP coordinates".into());
    }
    if traps.iter().any(|t| {
        t.bottom > t.top
            || t.top > n
            || t.left > t.right
            || t.left < -(bound as isize)
            || t.right > bound as isize
    }) {
        return Err("invalid trapezoid coordinates".into());
    }
    let mut previous_end = 0;
    for contig in &sequence.contigs {
        let end = contig
            .start
            .checked_add(contig.len)
            .ok_or("contig end overflow")?;
        if contig.start < previous_end || end > n {
            return Err("invalid contig coordinates".into());
        }
        previous_end = end;
    }
    let mut sorted: Vec<_> = traps.to_vec();
    sorted.sort_by_key(|t| t.bottom);
    let min_len = params.min_hit_length as isize;
    let limit = params.max_hit_length as isize;
    let candidates: Vec<_> = (0..sorted.len())
        .into_par_iter()
        .map(|current| align_candidate(sequence, &sorted, current, params, min_len, limit))
        .collect();
    let mut covered = vec![false; sorted.len()];
    let mut hits = Vec::new();
    for (current, candidate) in candidates.into_iter().enumerate() {
        if covered[current] {
            continue;
        }
        for j in candidate.covered {
            covered[j] = true;
        }
        hits.extend(candidate.hits);
    }
    Ok(hits)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::Contig;

    fn sequence(bases: &[u8]) -> Sequence {
        Sequence {
            bases: bases.to_vec(),
            contigs: vec![Contig {
                start: 0,
                len: bases.len(),
                label: "test".into(),
            }],
        }
    }

    fn repeat() -> (Sequence, Trapezoid) {
        let motif = b"ACGTTGCATGTCAGTACGATCGTA";
        let mut bases = motif.to_vec();
        bases.extend_from_slice(b"GGGCCCCAAAATTTGG");
        let bottom = bases.len();
        bases.extend_from_slice(motif);
        (
            sequence(&bases),
            Trapezoid {
                bottom,
                top: bases.len(),
                left: bottom as isize,
                right: bottom as isize,
            },
        )
    }

    fn align_with_threads(sequence: &Sequence, traps: &[Trapezoid], threads: usize) -> Vec<DpHit> {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        pool.install(|| align_traps(sequence, traps, &DpParams::default()).unwrap())
    }

    #[test]
    fn exact_self_repeat_and_half_open_ends() {
        let (seq, trap) = repeat();
        let hits = align_traps(&seq, &[trap], &DpParams::default()).unwrap();
        assert_eq!(hits.len(), 1);
        let hit = &hits[0];
        assert_eq!(
            (hit.a_lo, hit.a_hi, hit.b_lo, hit.b_hi, hit.score, hit.error),
            (0, 24, 40, 64, 24, 0.0)
        );
        assert!(hit.low_diagonal <= -40 && hit.high_diagonal >= -40);
        assert!(hit.low_diagonal < hit.high_diagonal);
    }

    #[test]
    fn identity_and_length_parameters() {
        let (mut seq, trap) = repeat();
        seq.bases[52] = b'T';
        let hits = align_traps(&seq, &[trap], &DpParams::default()).unwrap();
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].score, 20);
        assert!((hits[0].error - 1.0 / 24.0).abs() < 1e-6);
        let strict = DpParams {
            min_identity: 1.0,
            ..DpParams::default()
        };
        assert!(align_traps(&seq, &[trap], &strict).unwrap().is_empty());
        let long = DpParams {
            min_hit_length: 25,
            ..DpParams::default()
        };
        assert!(align_traps(&seq, &[trap], &long).unwrap().is_empty());
        let short = DpParams {
            max_hit_length: 16,
            ..DpParams::default()
        };
        assert!(
            align_traps(&seq, &[trap], &short)
                .unwrap()
                .iter()
                .all(|h| h.b_hi - h.b_lo <= 16)
        );
    }

    #[test]
    fn rejects_cross_contig_and_padding() {
        let (mut seq, trap) = repeat();
        seq.contigs = vec![
            Contig {
                start: 0,
                len: 24,
                label: "a".into(),
            },
            Contig {
                start: 40,
                len: 24,
                label: "b".into(),
            },
        ];
        assert!(
            align_traps(&seq, &[trap], &DpParams::default())
                .unwrap()
                .is_empty()
        );
        seq.contigs.truncate(1);
        assert!(
            align_traps(&seq, &[trap], &DpParams::default())
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn empty_short_and_diagonal_gate() {
        assert!(
            align_traps(&sequence(b""), &[], &DpParams::default())
                .unwrap()
                .is_empty()
        );
        let short = Trapezoid {
            bottom: 8,
            top: 12,
            left: 8,
            right: 8,
        };
        assert!(
            align_traps(&sequence(b"ACGTACGTACGT"), &[short], &DpParams::default())
                .unwrap()
                .is_empty()
        );
        let (seq, mut trap) = repeat();
        trap.left = 7;
        assert!(
            align_traps(&seq, &[trap], &DpParams::default())
                .unwrap()
                .is_empty()
        );
        trap.left = 100;
        trap.right = 100;
        assert!(
            align_traps(&seq, &[trap], &DpParams::default())
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn recursion_finds_separated_repeats() {
        let first = b"ACGTTGCATGTCAGTACGATCGTA";
        let second = b"TGCACCTGATCGTAGCATTCAGGT";
        let mut bases = first.to_vec();
        bases.extend_from_slice(&[b'A'; 30]);
        bases.extend_from_slice(second);
        bases.extend_from_slice(&[b'G'; 20]);
        bases.extend_from_slice(first);
        bases.extend_from_slice(&[b'C'; 30]);
        bases.extend_from_slice(second);
        let trap = Trapezoid {
            bottom: 98,
            top: 176,
            left: 98,
            right: 98,
        };
        let hits = align_traps(&sequence(&bases), &[trap], &DpParams::default()).unwrap();
        assert_eq!(hits.len(), 2);
        assert_eq!(
            (hits[0].a_lo, hits[0].a_hi, hits[0].b_lo, hits[0].b_hi),
            (0, 24, 98, 122)
        );
        assert_eq!(
            (hits[1].a_lo, hits[1].a_hi, hits[1].b_lo, hits[1].b_hi),
            (54, 78, 152, 176)
        );
        for hit in hits {
            assert_eq!(
                (hit.low_diagonal, hit.high_diagonal, hit.score, hit.error),
                (-105, -92, 24, 0.0)
            );
        }
    }

    #[test]
    fn covered_traps_and_reference_byte_equality() {
        let (seq, trap) = repeat();
        let single = align_traps(&seq, &[trap], &DpParams::default()).unwrap();
        assert_eq!(
            align_traps(&seq, &[trap, trap], &DpParams::default()).unwrap(),
            single
        );
        let (mut seq, trap) = repeat();
        seq.bases[10] = b'N';
        seq.bases[50] = b'N';
        assert_eq!(
            align_traps(&seq, &[trap], &DpParams::default()).unwrap()[0].score,
            24
        );
        seq.bases[50] = b'n';
        assert_eq!(
            align_traps(&seq, &[trap], &DpParams::default()).unwrap()[0].score,
            20
        );
    }

    #[test]
    fn nested_covered_trap_is_skipped_in_every_pool() {
        let (seq, trap) = repeat();
        let mut nested = trap;
        nested.top = trap.bottom + 20;
        let nested_alone = align_traps(&seq, &[nested], &DpParams::default()).unwrap();
        assert_eq!(nested_alone.len(), 1);
        let expected = align_traps(&seq, &[trap], &DpParams::default()).unwrap();
        assert_eq!(nested_alone, expected);
        // Equal bottoms keep their input order: whichever duplicate comes first
        // covers the other, so every pool returns the same single hit.
        for traps in [[trap, nested], [nested, trap]] {
            for threads in [1, 2, 4] {
                assert_eq!(
                    align_with_threads(&seq, &traps, threads),
                    expected,
                    "nested trap must be covered with {threads} thread(s)"
                );
            }
        }
    }

    #[test]
    fn rayon_pools_preserve_recursive_hit_order() {
        let first = b"ACGTTGCATGTCAGTACGATCGTA";
        let second = b"TGCACCTGATCGTAGCATTCAGGT";
        let mut bases = first.to_vec();
        bases.extend_from_slice(&[b'A'; 30]);
        bases.extend_from_slice(second);
        bases.extend_from_slice(&[b'G'; 20]);
        bases.extend_from_slice(first);
        bases.extend_from_slice(&[b'C'; 30]);
        bases.extend_from_slice(second);
        let seq = sequence(&bases);
        let trap = Trapezoid {
            bottom: 98,
            top: 176,
            left: 98,
            right: 98,
        };
        let traps = [trap, trap];
        let single = align_with_threads(&seq, &traps, 1);
        assert_eq!(single.len(), 2);
        for threads in [2, 4, 8] {
            assert_eq!(
                align_with_threads(&seq, &traps, threads),
                single,
                "hit order changed with {threads} threads"
            );
        }
        assert_eq!(
            single,
            align_traps(&seq, &traps, &DpParams::default()).unwrap()
        );
    }

    #[test]
    fn draft_acceptance_keeps_reference_inclusive_length() {
        let motif: Vec<_> = (0..128)
            .map(|i| b"ACGTTGCATGTCAGTACGATCGTA"[i % 24])
            .collect();
        for (length, accepted) in [(127, true), (128, false)] {
            let mut bases = motif[..length].to_vec();
            bases.extend_from_slice(&[b'!'; 32]);
            let bottom = bases.len();
            bases.extend_from_slice(&motif[..length]);
            let trap = Trapezoid {
                bottom,
                top: bases.len(),
                left: bottom as isize,
                right: bottom as isize,
            };
            let hits = align_traps(&sequence(&bases), &[trap], &DpParams::default()).unwrap();
            assert_eq!(!hits.is_empty(), accepted);
        }
    }

    #[test]
    fn insertion_matches_reference_score_and_error() {
        let (mut seq, mut trap) = repeat();
        seq.bases.insert(52, b'C');
        seq.contigs[0].len += 1;
        trap.top += 1;
        trap.left -= 1;
        trap.right += 1;
        let hits = align_traps(&seq, &[trap], &DpParams::default()).unwrap();
        assert_eq!(hits.len(), 1);
        assert_eq!(
            hits[0],
            DpHit {
                a_lo: 0,
                a_hi: 24,
                b_lo: 40,
                b_hi: 65,
                low_diagonal: -48,
                high_diagonal: -34,
                score: 21,
                error: 0.05,
            }
        );
    }

    #[test]
    fn invalid_parameters_and_coordinates() {
        let (seq, trap) = repeat();
        for params in [
            DpParams {
                min_hit_length: 0,
                ..DpParams::default()
            },
            DpParams {
                max_hit_length: 15,
                ..DpParams::default()
            },
            DpParams {
                min_identity: f64::NAN,
                ..DpParams::default()
            },
            DpParams {
                min_identity: -0.1,
                ..DpParams::default()
            },
            DpParams {
                min_identity: 1.1,
                ..DpParams::default()
            },
        ] {
            assert!(align_traps(&seq, &[trap], &params).is_err());
        }
        for bad in [
            Trapezoid {
                top: seq.bases.len() + 1,
                ..trap
            },
            Trapezoid {
                bottom: trap.top + 1,
                ..trap
            },
            Trapezoid {
                left: 42,
                right: 40,
                ..trap
            },
        ] {
            assert!(align_traps(&seq, &[bad], &DpParams::default()).is_err());
        }
    }
}

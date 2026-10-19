use crate::kmer::{FilterHit, FilterParams};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Trapezoid {
    pub bottom: usize,
    pub top: usize,
    pub left: isize,
    pub right: isize,
}

#[derive(Clone, Copy)]
struct Trap {
    bottom: i128,
    top: i128,
    left: i128,
    right: i128,
}

pub fn merge_hits(
    sequence: &[u8],
    hits: &[FilterHit],
    params: &FilterParams,
) -> Result<Vec<Trapezoid>, String> {
    if hits.is_empty() {
        return Ok(Vec::new());
    }
    if params.word_size == 0 || params.tube_offset == 0 {
        return Err("word_size and tube_offset must be positive".into());
    }
    if hits
        .iter()
        .any(|hit| hit.q_from > sequence.len() || hit.q_to > sequence.len())
    {
        return Err("filter hit coordinates are outside the sequence".into());
    }

    let bin_width = params.tube_offset as i128 + params.seed_diffs as i128 - 1;
    let left_padding = 2 + bin_width;
    let bottom_padding = params.word_size as i128 + 2;
    let mut sorted: Vec<_> = hits.iter().collect();
    sorted.sort_by_key(|hit| hit.q_from);
    let mut active: Vec<Trap> = Vec::new();
    let mut finished = Vec::new();

    for hit in sorted {
        let diagonal = -(hit.diagonal as i128);
        if diagonal <= params.seed_diffs as i128 {
            continue;
        }
        let bottom = hit.q_from as i128;
        let top = hit.q_to as i128;
        let threshold = bottom - bottom_padding;
        let mut index = 0;
        loop {
            if index == active.len() {
                active.push(Trap {
                    bottom,
                    top,
                    left: diagonal,
                    right: diagonal + bin_width,
                });
                break;
            }
            if active[index].top < threshold {
                finished.push(active.remove(index));
            } else if diagonal > active[index].right + 2 {
                index += 1;
            } else if diagonal >= active[index].left - left_padding {
                active[index].right = active[index].right.max(diagonal + bin_width);
                active[index].left = active[index].left.min(diagonal);
                active[index].top = active[index].top.max(top);

                if index > 0 && active[index - 1].right + 2 >= active[index].left {
                    let trap = active.remove(index);
                    let previous = &mut active[index - 1];
                    previous.right = trap.right;
                    previous.bottom = previous.bottom.min(trap.bottom);
                    previous.top = previous.top.max(trap.top);
                } else if index + 1 < active.len()
                    && active[index + 1].left - 2 <= active[index].right
                {
                    let next = active.remove(index + 1);
                    let trap = &mut active[index];
                    trap.right = next.right;
                    trap.bottom = trap.bottom.min(next.bottom);
                    trap.top = trap.top.max(next.top);
                }
                break;
            } else {
                active.insert(
                    index,
                    Trap {
                        bottom,
                        top,
                        left: diagonal,
                        right: diagonal + bin_width,
                    },
                );
                break;
            }
        }
    }

    finished.extend(active);
    finished.reverse();
    finished
        .into_iter()
        .map(|trap| {
            Ok(Trapezoid {
                bottom: usize::try_from(trap.bottom)
                    .map_err(|_| "trapezoid bottom is out of range")?,
                top: usize::try_from(trap.top).map_err(|_| "trapezoid top is out of range")?,
                left: isize::try_from(trap.left)
                    .map_err(|_| "trapezoid left diagonal is out of range")?,
                right: isize::try_from(trap.right)
                    .map_err(|_| "trapezoid right diagonal is out of range")?,
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params() -> FilterParams {
        FilterParams {
            word_size: 4,
            seed_length: 12,
            seed_diffs: 1,
            tube_offset: 3,
            diameter: 100,
        }
    }

    fn hit(q_from: usize, q_to: usize, diagonal: isize) -> FilterHit {
        FilterHit {
            q_from,
            q_to,
            diagonal: -diagonal,
        }
    }

    fn trap(bottom: usize, top: usize, left: isize, right: isize) -> Trapezoid {
        Trapezoid {
            bottom,
            top,
            left,
            right,
        }
    }

    #[test]
    fn diagonal_boundaries_are_inclusive() {
        let sequence = [b'A'; 100];
        for (diagonal, expected) in [
            (15, vec![trap(20, 35, 10, 18)]),
            (16, vec![trap(21, 35, 16, 19), trap(20, 30, 10, 13)]),
            (5, vec![trap(20, 35, 5, 13)]),
            (4, vec![trap(20, 30, 10, 13), trap(21, 35, 4, 7)]),
        ] {
            assert_eq!(
                merge_hits(
                    &sequence,
                    &[hit(20, 30, 10), hit(21, 35, diagonal)],
                    &params()
                )
                .unwrap(),
                expected
            );
        }
    }

    #[test]
    fn vertical_boundary_uses_word_size_plus_two() {
        let sequence = [b'A'; 100];
        assert_eq!(
            merge_hits(&sequence, &[hit(10, 20, 10), hit(26, 35, 10)], &params()).unwrap(),
            vec![trap(10, 35, 10, 13)]
        );
        assert_eq!(
            merge_hits(&sequence, &[hit(10, 20, 10), hit(27, 35, 10)], &params()).unwrap(),
            vec![trap(27, 35, 10, 13), trap(10, 20, 10, 13)]
        );
    }

    #[test]
    fn bridge_merges_the_next_active_trap() {
        assert_eq!(
            merge_hits(
                &[b'A'; 100],
                &[hit(20, 40, 10), hit(21, 45, 20), hit(22, 35, 15)],
                &params(),
            )
            .unwrap(),
            vec![trap(20, 45, 10, 23)]
        );
    }

    #[test]
    fn a_bridge_merges_only_one_neighbor() {
        let mut params = params();
        params.tube_offset = 20;
        assert_eq!(
            merge_hits(
                &[b'A'; 200],
                &[
                    hit(50, 90, 10),
                    hit(51, 91, 40),
                    hit(52, 92, 70),
                    hit(53, 93, 31),
                ],
                &params,
            )
            .unwrap(),
            vec![trap(52, 92, 70, 90), trap(50, 93, 10, 60)]
        );
    }

    #[test]
    fn expiration_follows_scan_order_not_a_global_sweep() {
        assert_eq!(
            merge_hits(
                &[b'A'; 100],
                &[hit(10, 50, 10), hit(11, 20, 30), hit(30, 40, 10)],
                &params(),
            )
            .unwrap(),
            vec![trap(11, 20, 30, 33), trap(10, 50, 10, 13)]
        );
    }

    #[test]
    fn bridge_merges_an_expired_successor_without_scanning_it() {
        assert_eq!(
            merge_hits(
                &[b'A'; 200],
                &[hit(50, 100, 10), hit(51, 60, 20), hit(80, 100, 15)],
                &params(),
            )
            .unwrap(),
            vec![trap(50, 100, 10, 23)]
        );
    }

    #[test]
    fn expired_traps_are_prepended_in_scan_order() {
        assert_eq!(
            merge_hits(
                &[b'A'; 200],
                &[hit(50, 60, 10), hit(51, 61, 30), hit(80, 100, 50)],
                &params(),
            )
            .unwrap(),
            vec![
                trap(80, 100, 50, 53),
                trap(51, 61, 30, 33),
                trap(50, 60, 10, 13)
            ]
        );
    }

    #[test]
    fn self_diagonals_are_excluded_and_hits_are_sorted() {
        assert_eq!(
            merge_hits(
                &[b'A'; 100],
                &[
                    hit(35, 45, 20),
                    hit(10, 20, 1),
                    hit(11, 21, 0),
                    hit(12, 22, -10),
                    hit(20, 30, 10),
                ],
                &params(),
            )
            .unwrap(),
            vec![trap(35, 45, 20, 23), trap(20, 30, 10, 13)]
        );
    }

    #[test]
    fn all_bytes_preserve_traps_without_ambiguity_splitting_or_clipping() {
        let hits = [hit(0, 2, 10), hit(3, 100, 30)];
        let expected = vec![trap(3, 100, 30, 33), trap(0, 2, 10, 13)];
        for byte in 0..=255 {
            assert_eq!(
                merge_hits(&[byte; 100], &hits, &params()).unwrap(),
                expected
            );
        }
    }

    #[test]
    fn reversed_hits_preserve_reference_geometry_and_expiration() {
        assert_eq!(
            merge_hits(&[b'A'; 100], &[hit(30, 20, 10)], &params()).unwrap(),
            vec![trap(30, 20, 10, 13)]
        );
        assert_eq!(
            merge_hits(&[b'A'; 100], &[hit(30, 20, 10), hit(31, 40, 15)], &params(),).unwrap(),
            vec![trap(31, 40, 15, 18), trap(30, 20, 10, 13)]
        );
    }

    #[test]
    fn empty_and_invalid_inputs() {
        assert!(merge_hits(&[], &[], &params()).unwrap().is_empty());
        assert!(merge_hits(&[b'A'; 10], &[hit(11, 4, 3)], &params()).is_err());
        assert!(merge_hits(&[b'A'; 10], &[hit(0, 11, 3)], &params()).is_err());
        let mut params = params();
        params.tube_offset = 0;
        assert!(merge_hits(&[b'A'; 10], &[hit(0, 10, 3)], &params).is_err());
    }

    #[test]
    #[ignore = "requires g++ and the reference source tree"]
    fn differential_against_compiled_reference() {
        use std::io::{BufRead, BufReader, Write};
        use std::process::{Command, Stdio};

        let reference = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../ref/pilercr1.06-src")
            .canonicalize()
            .unwrap();
        let executable = std::env::temp_dir().join(format!(
            "pilercr-merge-reference-{}{}",
            std::process::id(),
            std::env::consts::EXE_SUFFIX
        ));
        let mut compiler = Command::new("g++")
            .args([
                "-std=c++11",
                "-O2",
                "-ffunction-sections",
                "-fdata-sections",
                "-Wl,--gc-sections",
                "-I",
            ])
            .arg(&reference)
            .arg(reference.join("utils.cpp"))
            .arg(reference.join("mergefilterhits.cpp"))
            .args(["-x", "c++", "-", "-o"])
            .arg(&executable)
            .stdin(Stdio::piped())
            .spawn()
            .unwrap();
        compiler
            .stdin
            .take()
            .unwrap()
            .write_all(
                br#"#include "pilercr.h"
#include <iostream>
#include <sstream>
int k;
void Log(const char *, ...) {}
void Quit(const char *format, ...) {
    va_list args;
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    abort();
}
int main() {
    std::string line;
    while (std::getline(std::cin, line)) {
        std::istringstream input(line);
        int length, count;
        FilterParams p;
        input >> length >> p.WordSize >> p.SeedLength >> p.SeedDiffs >> p.TubeOffset >> count;
        k = p.WordSize;
        std::string sequence(length, 'A');
        for (int i = 0; i < length; ++i) sequence[i] = (char) i;
        std::vector<FilterHit> hits(count);
        for (int i = 0; i < count; ++i)
            input >> hits[i].QFrom >> hits[i].QTo >> hits[i].DiagIndex;
        int trap_count = -1;
        Trapezoid *traps = MergeFilterHits(sequence.data(), length, sequence.data(), length,
            true, hits.data(), count, p, &trap_count);
        std::cout << trap_count;
        for (int i = 0; i < trap_count; ++i) {
            if (!traps) abort();
            std::cout << ' ' << traps->bot << ' ' << traps->top
                << ' ' << traps->lft << ' ' << traps->rgt;
            traps = traps->next;
        }
        std::cout << '\n';
    }
}
"#,
            )
            .unwrap();
        assert!(compiler.wait().unwrap().success());

        let mut state = 0x6d65_7267_6568_6974_u64;
        let mut random = |bound: usize| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state as usize % bound
        };
        let mut input = String::new();
        let mut expected = Vec::new();
        for case in 0..50_007 {
            let word_size = 1 + random(12);
            let seed_diffs = random(5);
            let tube_offset = seed_diffs.max(1) + random(20);
            let generated_params = FilterParams {
                word_size,
                seed_length: word_size * (seed_diffs + 1) + random(30),
                seed_diffs,
                tube_offset,
                diameter: 512,
            };
            let sequence: Vec<_> = (0..512).map(|i| i as u8).collect();
            let count = random(81);
            let hits: Vec<_> = (0..count)
                .map(|_| {
                    let q_from = match case % 3 {
                        0 => 100 + random(20),
                        1 => 20 + random(460),
                        _ => 100 + 10 * random(30),
                    };
                    let q_to = q_from + word_size + random(513 - q_from - word_size);
                    let diagonal = random(q_from / tube_offset + 1) * tube_offset;
                    hit(q_from, q_to, diagonal as isize)
                })
                .collect();
            let (params, hits) = match case {
                50_000 => (
                    params(),
                    vec![hit(20, 40, 10), hit(21, 45, 20), hit(22, 35, 15)],
                ),
                50_001 => {
                    let mut p = params();
                    p.tube_offset = 20;
                    (
                        p,
                        vec![
                            hit(50, 90, 10),
                            hit(51, 91, 40),
                            hit(52, 92, 70),
                            hit(53, 93, 31),
                        ],
                    )
                }
                50_002 => (
                    params(),
                    vec![hit(10, 50, 10), hit(11, 20, 30), hit(30, 40, 10)],
                ),
                50_003 => (
                    params(),
                    vec![hit(50, 100, 10), hit(51, 60, 20), hit(80, 100, 15)],
                ),
                50_004 => (
                    params(),
                    vec![hit(50, 60, 10), hit(51, 61, 30), hit(80, 100, 50)],
                ),
                50_005 => (
                    params(),
                    vec![hit(20, 30, 10), hit(20, 35, 20), hit(20, 40, 15)],
                ),
                50_006 => (params(), vec![hit(30, 20, 10)]),
                _ => (generated_params, hits),
            };
            use std::fmt::Write;
            write!(
                input,
                "512 {} {} {} {} {}",
                params.word_size,
                params.seed_length,
                params.seed_diffs,
                params.tube_offset,
                hits.len()
            )
            .unwrap();
            for hit in &hits {
                write!(input, " {} {} {}", hit.q_from, hit.q_to, hit.diagonal).unwrap();
            }
            input.push('\n');
            expected.push(merge_hits(&sequence, &hits, &params).unwrap());
        }
        let mut child = Command::new(&executable)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .spawn()
            .unwrap();
        let mut stdin = child.stdin.take().unwrap();
        let writer = std::thread::spawn(move || stdin.write_all(input.as_bytes()).unwrap());
        let mut compared = 0;
        let mut mismatches = 0;
        let mut first = None;
        for (case, line) in BufReader::new(child.stdout.take().unwrap())
            .lines()
            .enumerate()
        {
            let values: Vec<i128> = line
                .unwrap()
                .split_whitespace()
                .map(|value| value.parse().unwrap())
                .collect();
            assert!(values[0] >= 0);
            assert_eq!(values.len(), 1 + values[0] as usize * 4);
            assert_eq!(values[0] as usize, expected[case].len(), "case {case}");
            let actual: Vec<_> = values[1..]
                .as_chunks::<4>()
                .0
                .iter()
                .map(|v| trap(v[0] as usize, v[1] as usize, v[2] as isize, v[3] as isize))
                .collect();
            if actual != expected[case] {
                mismatches += 1;
                first.get_or_insert_with(|| (case, actual, expected[case].clone()));
            }
            compared += 1;
        }
        writer.join().unwrap();
        assert!(child.wait().unwrap().success());
        std::fs::remove_file(executable).unwrap();
        eprintln!(
            "merge differential: {compared} cases, {mismatches} mismatches; first: {first:?}"
        );
        assert_eq!(compared, expected.len());
        assert_eq!(mismatches, 0);
    }

    #[test]
    fn unrepresentable_diagonals_return_an_error() {
        let hits = [FilterHit {
            q_from: 0,
            q_to: 10,
            diagonal: isize::MIN,
        }];
        assert!(merge_hits(&[b'A'; 10], &hits, &params()).is_err());
    }
}

use crate::{graph::ComponentMember, msa, pile::Pile, seq::Sequence};

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Array {
    pub id: usize,
    pub pos: usize,
    pub consensus: Vec<u8>,
    pub aligned_consensus: Vec<u8>,
    pub left_flanks: Vec<Vec<u8>>,
    pub repeats: Vec<Vec<u8>>,
    pub spacers: Vec<Vec<u8>>,
    pub cluster_rev_comp: bool,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ArrayParams {
    pub min_repeat: usize,
    pub max_repeat: usize,
    pub min_spacer: usize,
    pub max_spacer: usize,
    pub min_array: usize,
    pub min_repeat_ratio: f64,
    pub min_spacer_ratio: f64,
    pub min_cons: f64,
}

impl Default for ArrayParams {
    fn default() -> Self {
        Self {
            min_repeat: 16,
            max_repeat: 64,
            min_spacer: 8,
            max_spacer: 64,
            min_array: 3,
            min_repeat_ratio: 0.9,
            min_spacer_ratio: 0.75,
            min_cons: 0.9,
        }
    }
}

fn ratio(x: i64, y: i64) -> f64 {
    if x == 0 && y == 0 {
        0.0
    } else if x < y {
        x as f64 / y as f64
    } else {
        y as f64 / x as f64
    }
}

fn substitution(a: u8, b: u8) -> i32 {
    let a = a.to_ascii_uppercase();
    let b = b.to_ascii_uppercase();
    if !b"ACGT".contains(&a) || !b"ACGT".contains(&b) {
        0
    } else if a == b {
        50
    } else {
        -100
    }
}

#[derive(Debug, PartialEq, Eq)]
struct LocalAlignment {
    score: i32,
    start_a: usize,
    start_b: usize,
    path: Vec<u8>,
}

fn sw_simple(a: &[u8], b: &[u8]) -> Result<LocalAlignment, String> {
    if a.is_empty() || b.is_empty() {
        return Err("Local alignment requires nonempty sequences".into());
    }
    let width = b.len() + 1;
    let cells = (a.len() + 1)
        .checked_mul(width)
        .ok_or("Local alignment matrix is too large")?;
    let mut m = vec![0i32; cells];
    let mut d = vec![0i32; cells];
    let mut ins = vec![0i32; cells];
    let (mut top_score, mut top_i, mut top_j) = (-1, 0, 0);
    for i in 1..=a.len() {
        for j in 1..=b.len() {
            let p = i * width + j;
            let diag = p - width - 1;
            let s = m[diag].max(d[diag]).max(ins[diag]) + substitution(a[i - 1], b[j - 1]);
            if s > 0 {
                m[p] = s;
                if s > top_score {
                    top_score = s;
                    top_i = i;
                    top_j = j;
                }
            }
            d[p] = (m[p - width] - 400).max(d[p - width] - 25).max(0);
            ins[p] = (m[p - 1] - 400).max(ins[p - 1] - 25).max(0);
        }
    }
    if top_score < 0 {
        return Ok(LocalAlignment {
            score: 0,
            start_a: 0,
            start_b: 0,
            path: Vec::new(),
        });
    }
    let (mut i, mut j, mut edge) = (top_i, top_j, b'M');
    let mut path = Vec::new();
    while i != 0 && j != 0 {
        let p = i * width + j;
        let previous = edge;
        match edge {
            b'M' => {
                let diag = p - width - 1;
                let s = substitution(a[i - 1], b[j - 1]);
                if m[diag] + s == m[p] {
                    edge = b'M';
                } else if i > 1 && d[diag] + s == m[p] {
                    edge = b'D';
                } else if j > 1 && ins[diag] + s == m[p] {
                    edge = b'I';
                } else if m[p] == 0 {
                    break;
                } else {
                    return Err("Local alignment traceback failed at M".into());
                }
                i -= 1;
                j -= 1;
            }
            b'D' => {
                if i > 1 && m[p - width] - 400 == d[p] {
                    edge = b'M';
                } else if i > 1 && d[p - width] - 25 == d[p] {
                    edge = b'D';
                } else {
                    return Err("Local alignment traceback failed at D".into());
                }
                i -= 1;
            }
            b'I' => {
                if j > 1 && m[p - 1] - 400 == ins[p] {
                    edge = b'M';
                } else if j > 1 && ins[p - 1] - 25 == ins[p] {
                    edge = b'I';
                } else {
                    return Err("Local alignment traceback failed at I".into());
                }
                j -= 1;
            }
            _ => unreachable!(),
        }
        path.push(previous);
    }
    path.reverse();
    Ok(LocalAlignment {
        score: top_score,
        start_a: i,
        start_b: j,
        path,
    })
}

fn pile_to_seq<'a>(sequence: &'a Sequence, pile: &Pile) -> Result<(usize, &'a [u8]), String> {
    let (ci, _) = sequence
        .global_to_local(pile.lo)
        .ok_or("Pile is outside a contig")?;
    let contig = &sequence.contigs[ci];
    let end = contig
        .start
        .checked_add(contig.len)
        .ok_or("Contig coordinate overflow")?;
    let lo = pile.lo.saturating_sub(8).max(contig.start);
    let hi = pile.hi.saturating_add(8).min(end - 1);
    let bases = sequence
        .bases
        .get(lo..=hi)
        .ok_or("Invalid padded pile interval")?;
    Ok((lo, bases))
}

fn spacer(indexes: &[usize], index: usize, piles: &[Pile]) -> i64 {
    indexes[index + 1..]
        .iter()
        .find(|&&p| !piles[p].deleted)
        .map_or(-1, |&p| {
            piles[p].lo as i64 - piles[indexes[index]].hi as i64
        })
}

fn get_array(indexes: &[usize], index: usize, piles: &mut [Pile], min_array: usize) -> Vec<usize> {
    let root = indexes[index];
    let root_spacer = spacer(indexes, index, piles);
    if piles[root].deleted || root_spacer < 8 {
        return Vec::new();
    }
    let root_length = (piles[root].hi - piles[root].lo + 1) as i64;
    let mut selected = vec![root];
    for prev in (0..index).rev() {
        let p = indexes[prev];
        if piles[p].deleted {
            continue;
        }
        let s = spacer(indexes, prev, piles);
        if s < 8
            || ratio(root_length, (piles[p].hi - piles[p].lo + 1) as i64) < 0.6
            || ratio(root_spacer + 16, s + 16) < 0.6
        {
            break;
        }
        selected.push(p);
    }
    for next in index + 1..indexes.len() {
        let p = indexes[next];
        if next == indexes.len() - 1 {
            selected.push(p);
            break;
        }
        if piles[p].deleted {
            continue;
        }
        if ratio(root_length, (piles[p].hi - piles[p].lo + 1) as i64) < 0.6
            || ratio(root_spacer + 16, spacer(indexes, next, piles) + 16) < 0.6
        {
            break;
        }
        selected.push(p);
    }
    if selected.len() < min_array {
        piles[root].deleted = true;
        selected.clear();
    }
    selected
}

fn recompute_consensus(array: &mut Array) -> Result<(), String> {
    let cols = array.repeats.first().map_or(0, Vec::len);
    if array.repeats.iter().any(|r| r.len() != cols) {
        return Err("Array alignment rows have unequal lengths".into());
    }
    array.aligned_consensus = (0..cols)
        .map(|col| {
            let mut counts = [0usize; 5];
            for row in &array.repeats {
                if let Some(letter) = b"ACGT-"
                    .iter()
                    .position(|&c| c == row[col].to_ascii_uppercase())
                {
                    counts[letter] += 1;
                }
            }
            let mut best = 0;
            for i in 1..5 {
                if counts[i] > counts[best] {
                    best = i;
                }
            }
            b"ACGT-"[best]
        })
        .collect();
    array.consensus = array
        .aligned_consensus
        .iter()
        .copied()
        .filter(|&c| c != b'-' && !c.is_ascii_whitespace())
        .collect();
    Ok(())
}

fn align_array(
    sequence: &Sequence,
    piles: &[Pile],
    indexes: &[usize],
    consensus: &[u8],
    id: usize,
    min_array: usize,
) -> Result<Array, String> {
    let mut sorted = indexes.to_vec();
    sorted.sort_by_key(|&i| piles[i].lo);
    let mut locations = Vec::with_capacity(sorted.len());
    let mut hits = Vec::with_capacity(sorted.len() + 1);
    for index in sorted {
        let (lo, bases) = pile_to_seq(sequence, &piles[index])?;
        let alignment = sw_simple(bases, consensus)?;
        let length = alignment
            .path
            .iter()
            .filter(|&&c| c == b'M' || c == b'D')
            .count();
        locations.push((lo + alignment.start_a, length));
        hits.push(bases[alignment.start_a..alignment.start_a + length].to_vec());
    }
    hits.push(consensus.to_vec());
    let alignment = msa::align(&hits)?;
    let mut array = Array {
        id,
        ..Array::default()
    };
    for (row, &(lo, length)) in locations.iter().enumerate() {
        if let Some(&(next, _)) = locations.get(row + 1)
            && next as i64 - (lo + length) as i64 + 1 < 4
        {
            return Ok(Array {
                id,
                ..Array::default()
            });
        }
        let (ci, _) = sequence
            .global_to_local(lo)
            .ok_or("Localized repeat is outside a contig")?;
        let contig = &sequence.contigs[ci];
        let flank_start = lo.saturating_sub(10).max(contig.start);
        let mut flank = vec![b' '; 10 - (lo - flank_start)];
        flank.extend_from_slice(&sequence.bases[flank_start..lo]);
        array.left_flanks.push(flank);
        array.repeats.push(alignment.rows[row].clone());
        let start = lo + length;
        let end = if row + 1 == locations.len() {
            start.saturating_add(10).min(contig.start + contig.len)
        } else {
            locations[row + 1].0
        };
        array.spacers.push(
            sequence
                .bases
                .get(start..end)
                .ok_or("Invalid array spacer interval")?
                .to_vec(),
        );
    }
    array.pos = locations[0].0;
    crate::bounds::trace_rows(&array, "before");
    let bounds_result = crate::bounds::fix_bounds(&mut array, sequence, min_array);
    crate::bounds::trace_rows(&array, "after");
    bounds_result?;
    recompute_consensus(&mut array)?;
    Ok(array)
}

fn ungapped_length(row: &[u8]) -> usize {
    row.iter().filter(|&&c| c != b'-').count()
}

fn average_repeat(array: &Array) -> usize {
    array.repeats.iter().map(Vec::len).sum::<usize>() / array.repeats.len()
}

fn average_spacer(array: &Array) -> usize {
    array.spacers[..array.repeats.len() - 1]
        .iter()
        .map(Vec::len)
        .sum::<usize>()
        / (array.repeats.len() - 1)
}

fn array_length(array: &Array) -> usize {
    array
        .repeats
        .iter()
        .map(|r| ungapped_length(r))
        .sum::<usize>()
        + array.spacers[..array.repeats.len() - 1]
            .iter()
            .map(Vec::len)
            .sum::<usize>()
}

fn accepted(array: &Array, params: &ArrayParams) -> bool {
    acceptance_reason(array, params) == "accepted"
}

fn acceptance_reason(array: &Array, params: &ArrayParams) -> &'static str {
    if array.repeats.len() < params.min_array || array.repeats.len() < 2 {
        return "repeat-count";
    }
    let spacer_length = average_spacer(array);
    let repeat_length = average_repeat(array);
    if spacer_length < params.min_spacer || spacer_length > params.max_spacer {
        return "spacer-length";
    }
    if repeat_length < params.min_repeat || repeat_length > params.max_repeat {
        return "repeat-length";
    }
    let mut run = 0;
    for repeat in &array.repeats {
        let matches = repeat
            .iter()
            .zip(&array.aligned_consensus)
            .filter(|(a, b)| a.eq_ignore_ascii_case(b))
            .count();
        if matches as f64 / repeat.len() as f64 >= params.min_cons {
            run += 1;
        } else {
            run = 0;
        }
    }
    if run < params.min_array {
        return "conservation";
    }
    if 1.0 >= params.min_spacer_ratio {
        "accepted"
    } else {
        "spacer-ratio"
    }
}

fn merge_arrays(sequence: &Sequence, first: &mut Array, second: &Array) -> Result<bool, String> {
    let repeat_length = average_repeat(first);
    if repeat_length != average_repeat(second) {
        return Ok(false);
    }
    let spacer_length = average_spacer(first);
    if ratio(spacer_length as i64, average_spacer(second) as i64) < 0.95 {
        return Ok(false);
    }
    let last_spacer_start = first.pos + array_length(first);
    let end = last_spacer_start + spacer_length - 1;
    let end_contig = sequence.global_to_local(end).map(|p| p.0);
    if end_contig.is_none() || end_contig != sequence.global_to_local(second.pos).map(|p| p.0) {
        return Ok(false);
    }
    let distance = second.pos as i64 - end as i64 - 1;
    let unit = (repeat_length + spacer_length) as i64;
    let missed = match (1..=4).find(|&n| ratio(distance, n * unit) >= 0.95) {
        Some(n) => n as usize,
        None => return Ok(false),
    };
    let alignment = sw_simple(&first.consensus, &second.consensus)?;
    if alignment.score <= 0
        || ratio(alignment.path.len() as i64, first.consensus.len() as i64) < 0.95
    {
        return Ok(false);
    }
    let mut starts = Vec::with_capacity(missed);
    let mut lengths = Vec::with_capacity(missed);
    let mut repeats = Vec::with_capacity(missed);
    let mut previous_end = last_spacer_start - 1;
    for _ in 0..missed {
        let guess = previous_end + spacer_length;
        let start = match guess.checked_sub(8) {
            Some(start) => start,
            None => return Ok(false),
        };
        let end = guess + repeat_length + 8;
        let bases = match sequence.bases.get(start..=end) {
            Some(bases) => bases,
            None => return Ok(false),
        };
        let alignment = sw_simple(bases, &first.consensus)?;
        if ratio(alignment.path.len() as i64, first.consensus.len() as i64) < 0.95 {
            return Ok(false);
        }
        let repeat_start = start + alignment.start_a;
        let mut pos = repeat_start;
        let mut repeat = vec![b'-'; alignment.start_b];
        for edge in alignment.path {
            if edge == b'M' || edge == b'D' {
                repeat.push(sequence.bases[pos]);
                pos += 1;
            } else {
                repeat.push(b'-');
            }
        }
        let length = pos - repeat_start;
        if lengths.first().is_some_and(|&l| l != length) {
            return Ok(false);
        }
        lengths.push(length);
        starts.push(repeat_start);
        previous_end = pos - 1;
        repeat.resize(repeat.len().max(first.consensus.len()), b'-');
        repeats.push(repeat);
    }
    let mut flanks = Vec::with_capacity(missed);
    let mut spacers = Vec::with_capacity(missed);
    for (i, &start) in starts.iter().enumerate() {
        let flank_start = match start.checked_sub(10) {
            Some(start) => start,
            None => return Ok(false),
        };
        flanks.push(sequence.bases[flank_start..start].to_vec());
        let spacer_start = start + ungapped_length(&repeats[i]);
        let spacer_end = starts.get(i + 1).copied().unwrap_or(second.pos);
        if spacer_end <= spacer_start + 1 {
            return Ok(false);
        }
        spacers.push(sequence.bases[spacer_start..spacer_end].to_vec());
    }
    if starts[0] <= last_spacer_start + 1 {
        return Ok(false);
    }
    *first.spacers.last_mut().ok_or("Array has no spacers")? =
        sequence.bases[last_spacer_start..starts[0]].to_vec();
    first.repeats.extend(repeats);
    first.left_flanks.extend(flanks);
    first.spacers.extend(spacers);
    first.repeats.extend(second.repeats.iter().cloned());
    first.left_flanks.extend(second.left_flanks.iter().cloned());
    first.spacers.extend(second.spacers.iter().cloned());
    Ok(true)
}

pub fn find_arrays(
    sequence: &Sequence,
    piles: &mut [Pile],
    components: &[Vec<ComponentMember>],
    params: &ArrayParams,
) -> Result<Vec<Array>, String> {
    for component in components {
        for member in component {
            let pile = piles
                .get(member.pile_index)
                .ok_or("Component contains an invalid pile index")?;
            if pile.lo > pile.hi || sequence.global_to_local(pile.lo).is_none() {
                return Err("Invalid pile interval".into());
            }
        }
    }
    let draft_min_array = params.min_array.saturating_sub(2).max(3);
    let mut drafts = Vec::new();
    let mut id = 0;
    for component in components {
        if component.len() < draft_min_array {
            continue;
        }
        let mut indexes: Vec<_> = component.iter().map(|m| m.pile_index).collect();
        crate::sort::reference_sort_by_key(&mut indexes, |&p| piles[p].lo);
        let mut by_spacer: Vec<_> = (0..indexes.len()).collect();
        // Equal spacer ordering selects the draft root and MSA input order.
        crate::sort::reference_sort_by_key(&mut by_spacer, |&i| spacer(&indexes, i, piles));
        if std::env::var_os("TRACE_ARRAY").is_some() {
            eprintln!(
                "component {:?}",
                indexes
                    .iter()
                    .enumerate()
                    .map(|(i, &p)| (piles[p].lo, piles[p].hi, spacer(&indexes, i, piles)))
                    .collect::<Vec<_>>()
            );
        }

        for index in by_spacer {
            let selected = get_array(&indexes, index, piles, draft_min_array);
            if selected.is_empty() {
                continue;
            }
            id += 1;
            let sequences: Vec<_> = selected
                .iter()
                .map(|&p| pile_to_seq(sequence, &piles[p]).map(|(_, s)| s.to_vec()))
                .collect::<Result<_, _>>()?;
            let alignment = msa::align(&sequences)?;
            let min_cons = if sequences.len() >= 3 {
                0.95f64.min((sequences.len() - 1) as f64 / sequences.len() as f64 - 0.1)
            } else {
                0.95
            };
            let (consensus, _) = alignment.conserved_consensus(min_cons)?;
            if std::env::var_os("TRACE_ARRAY").is_some() {
                eprintln!(
                    "draft {id} piles={:?} consensus={}",
                    selected
                        .iter()
                        .map(|&p| (piles[p].lo, piles[p].hi))
                        .collect::<Vec<_>>(),
                    String::from_utf8_lossy(&consensus)
                );
            }
            if consensus.len() < 16 {
                continue;
            }
            for &p in &selected {
                piles[p].deleted = true;
            }
            if crate::trace::enabled() {
                let mut pile_indexes = selected.clone();
                pile_indexes.sort_by_key(|&p| piles[p].lo);
                crate::trace::emit(
                    "12-draft-arrays",
                    &[
                        id.to_string(),
                        pile_indexes
                            .iter()
                            .map(usize::to_string)
                            .collect::<Vec<_>>()
                            .join(","),
                        String::from_utf8_lossy(&consensus).into_owned(),
                    ],
                );
            }
            drafts.push((id, selected, consensus));
        }
    }
    let mut arrays = Vec::new();
    for (id, indexes, consensus) in drafts {
        let array = align_array(sequence, piles, &indexes, &consensus, id, params.min_array)?;
        let is_accepted = !array.repeats.is_empty() && accepted(&array, params);
        if crate::trace::enabled() {
            crate::trace::emit(
                "15-decisions",
                &[
                    array.id.to_string(),
                    if is_accepted { "accept" } else { "reject" }.into(),
                    acceptance_reason(&array, params).into(),
                ],
            );
        }
        if is_accepted {
            arrays.push(array);
        }
    }
    arrays.sort_by_key(|a| a.pos);
    let mut index = 0;
    while index + 1 < arrays.len() {
        let (left, right) = arrays.split_at_mut(index + 1);
        let merged = merge_arrays(sequence, &mut left[index], &right[0])?;
        if crate::trace::enabled() {
            crate::trace::emit(
                "15-decisions",
                &[
                    left[index].id.to_string(),
                    if merged { "merge" } else { "keep" }.into(),
                    right[0].id.to_string(),
                ],
            );
        }
        if merged {
            arrays.remove(index + 1);
        }
        index += 1;
    }
    for (index, array) in arrays.iter_mut().enumerate() {
        array.id = index + 1;
    }
    Ok(arrays)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::Contig;

    fn array(rows: &[&[u8]], spacer: usize) -> Array {
        let mut array = Array {
            repeats: rows.iter().map(|r| r.to_vec()).collect(),
            left_flanks: vec![vec![b' '; 10]; rows.len()],
            spacers: vec![vec![b'N'; spacer]; rows.len()],
            ..Array::default()
        };
        recompute_consensus(&mut array).unwrap();
        array
    }

    #[test]
    fn local_alignment_reference_scores_and_ties() {
        let alignment = sw_simple(b"TTACGTGGACGT", b"ACGT").unwrap();
        assert_eq!(
            alignment,
            LocalAlignment {
                score: 200,
                start_a: 2,
                start_b: 0,
                path: b"MMMM".to_vec()
            }
        );
        assert_eq!(sw_simple(b"NN", b"NN").unwrap().score, 0);
        assert!(sw_simple(b"A", b"C").unwrap().path.is_empty());
        let mut a = vec![b'A'; 16];
        a.extend_from_slice(b"CCCC");
        a.extend_from_slice(&[b'T'; 16]);
        let mut b = vec![b'A'; 16];
        b.extend_from_slice(&[b'T'; 16]);
        let alignment = sw_simple(&a, &b).unwrap();
        let mut expected = vec![b'M'; 16];
        expected.extend_from_slice(b"DDDD");
        expected.extend_from_slice(&[b'M'; 16]);
        assert_eq!(alignment.path, expected);
        assert_eq!(alignment.score, 1125);
        let reverse = sw_simple(&b, &a).unwrap();
        for edge in &mut expected {
            if *edge == b'D' {
                *edge = b'I';
            }
        }
        assert_eq!(reverse.path, expected);
    }

    #[test]
    fn consensus_uses_existing_columns_and_reference_ties() {
        let a = array(&[b"A-CN", b"C-AN", b"G-TN"], 8);
        assert_eq!(a.aligned_consensus, b"A-AA");
        assert_eq!(a.consensus, b"AAA");
    }

    #[test]
    fn acceptance_uses_trailing_run_gapped_lengths_and_ignores_repeat_ratio() {
        let mut params = ArrayParams {
            min_repeat: 4,
            max_repeat: 4,
            min_spacer: 8,
            max_spacer: 8,
            min_cons: 1.0,
            min_repeat_ratio: 1.0,
            ..ArrayParams::default()
        };
        let mut a = array(&[b"AC-G", b"AC-G", b"AC-G", b"TTTT"], 8);
        assert!(!accepted(&a, &params));
        a.repeats.reverse();
        a.spacers.last_mut().unwrap().resize(100, b'N');
        assert!(accepted(&a, &params));
        params.min_repeat_ratio = 2.0;
        assert!(accepted(&a, &params));
        a.spacers[0].resize(1, b'N');
        a.spacers[1].resize(15, b'N');
        assert!(accepted(&a, &params));
    }

    #[test]
    fn acceptance_reasons_preserve_filter_results() {
        let mut a = array(&[b"ACGT", b"ACGT", b"ACGT"], 8);
        let mut params = ArrayParams {
            min_repeat: 4,
            max_repeat: 4,
            min_spacer: 8,
            max_spacer: 8,
            min_cons: 1.0,
            ..ArrayParams::default()
        };
        assert_eq!(acceptance_reason(&a, &params), "accepted");
        params.min_spacer_ratio = 2.0;
        assert_eq!(acceptance_reason(&a, &params), "spacer-ratio");
        a.aligned_consensus = b"TTTT".to_vec();
        assert_eq!(acceptance_reason(&a, &params), "conservation");
        params.min_repeat = 5;
        assert_eq!(acceptance_reason(&a, &params), "repeat-length");
        params.min_spacer = 9;
        assert_eq!(acceptance_reason(&a, &params), "spacer-length");
        a.repeats.pop();
        assert_eq!(acceptance_reason(&a, &params), "repeat-count");
    }

    #[test]
    fn padded_piles_stay_inside_contigs() {
        let sequence = Sequence {
            bases: b"ACGTNNNTGCA".to_vec(),
            contigs: vec![
                Contig {
                    label: "a".into(),
                    start: 0,
                    len: 4,
                },
                Contig {
                    label: "b".into(),
                    start: 7,
                    len: 4,
                },
            ],
        };
        let pile = Pile {
            lo: 8,
            hi: 9,
            deleted: false,
            hit_indexes: Vec::new(),
        };
        let (lo, bases) = pile_to_seq(&sequence, &pile).unwrap();
        assert_eq!(lo, 7);
        assert_eq!(bases, b"TGCA");
    }

    #[test]
    fn draft_final_pile_is_added_even_when_deleted() {
        let mut piles: Vec<_> = [0, 27, 54]
            .into_iter()
            .map(|lo| Pile {
                lo,
                hi: lo + 19,
                deleted: false,
                hit_indexes: Vec::new(),
            })
            .collect();
        piles[2].deleted = true;
        assert_eq!(get_array(&[0, 1, 2], 0, &mut piles, 3), vec![0, 1, 2]);
    }

    #[test]
    fn finds_arrays_through_alignment_and_bounds() {
        let repeat = b"ACGTTGCAAGTCCTAG";
        let mut bases = b"TTTTTTTTTTTTTTTTTTTT".to_vec();
        let mut piles = Vec::new();
        for spacer in [
            b"AAAAAAAAAAAAAAAA",
            b"CCCCCCCCCCCCCCCC",
            b"GGGGGGGGGGGGGGGG",
        ] {
            let lo = bases.len();
            bases.extend_from_slice(repeat);
            piles.push(Pile {
                lo,
                hi: bases.len() - 1,
                deleted: false,
                hit_indexes: Vec::new(),
            });
            bases.extend_from_slice(spacer);
        }
        let sequence = Sequence {
            contigs: vec![Contig {
                label: "a".into(),
                start: 0,
                len: bases.len(),
            }],
            bases,
        };
        let component = (0..3)
            .map(|pile_index| ComponentMember {
                pile_index,
                rev: false,
            })
            .collect();
        let directory = std::env::temp_dir().join(format!(
            "pilercr-array-trace-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let arrays = crate::trace::with_directory(&directory, || {
            find_arrays(&sequence, &mut piles, &[component], &ArrayParams::default())
        })
        .unwrap();
        let drafts = std::fs::read_to_string(directory.join("12-draft-arrays.tsv")).unwrap();
        assert_eq!(
            drafts.lines().skip(1).collect::<Vec<_>>(),
            ["1\t0,1,2\tACGTTGCAAGTCCTAGG"]
        );
        let bounds = std::fs::read_to_string(directory.join("14-bounds.tsv")).unwrap();
        let rows: Vec<_> = bounds.lines().skip(1).collect();
        assert_eq!(rows.len(), 6);
        assert_eq!(
            rows[0],
            "1\tbefore\t20\t0\tTTTTTTTTTT\tACGTTGCAAGTCCTAG-\tAAAAAAAAAAAAAAAA"
        );
        assert_eq!(
            rows[5],
            "1\tafter\t20\t2\tCCCCCCCCCC\tACGTTGCAAGTCCTAG\tGGGGGGGGGGG"
        );
        let decisions = std::fs::read_to_string(directory.join("15-decisions.tsv")).unwrap();
        assert_eq!(
            decisions.lines().skip(1).collect::<Vec<_>>(),
            ["1\taccept\taccepted"]
        );
        std::fs::remove_dir_all(&directory).unwrap();
        assert_eq!(arrays.len(), 1);
        assert_eq!(arrays[0].id, 1);
        assert_eq!(arrays[0].pos, 20);
        assert_eq!(arrays[0].consensus, repeat);
        assert_eq!(arrays[0].repeats.len(), 3);
        assert_eq!(arrays[0].spacers.last().unwrap(), b"GGGGGGGGGGG");
        assert!(piles.iter().all(|p| p.deleted));
    }

    #[test]
    fn merging_localizes_missing_repeats_without_recomputing_consensus() {
        let repeat = b"ACGTTGCAAGTCCTAG";
        let spacer = b"NNNNNNNNNNNNNNNN";
        let mut bases = vec![b'N'; 20];
        for _ in 0..7 {
            bases.extend_from_slice(repeat);
            bases.extend_from_slice(spacer);
        }
        let sequence = Sequence {
            contigs: vec![Contig {
                label: "a".into(),
                start: 0,
                len: bases.len(),
            }],
            bases,
        };
        let mut first = array(&[repeat, repeat, repeat], 16);
        first.pos = 20;
        let mut second = first.clone();
        second.pos = 20 + 4 * 32;
        let consensus = first.consensus.clone();
        assert!(merge_arrays(&sequence, &mut first, &second).unwrap());
        assert_eq!(first.repeats.len(), 7);
        assert_eq!(first.spacers.len(), 7);
        assert_eq!(first.left_flanks.len(), 7);
        assert_eq!(first.repeats[3], repeat);
        assert_eq!(first.spacers[2], spacer);
        assert_eq!(first.spacers[3], spacer);
        assert_eq!(first.consensus, consensus);
    }
}

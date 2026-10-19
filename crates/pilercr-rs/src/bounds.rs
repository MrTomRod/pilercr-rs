use crate::{array::Array, seq::Sequence};

pub(crate) fn trace_rows(array: &Array, phase: &str) {
    if !crate::trace::enabled() {
        return;
    }
    for (row_index, ((left_flank, repeat), spacer)) in array
        .left_flanks
        .iter()
        .zip(&array.repeats)
        .zip(&array.spacers)
        .enumerate()
    {
        crate::trace::emit(
            "14-bounds",
            &[
                array.id.to_string(),
                phase.into(),
                array.pos.to_string(),
                row_index.to_string(),
                String::from_utf8_lossy(left_flank).into_owned(),
                String::from_utf8_lossy(repeat).into_owned(),
                String::from_utf8_lossy(spacer).into_owned(),
            ],
        );
    }
}

pub fn fix_bounds(array: &mut Array, sequence: &Sequence, min_array: usize) -> Result<(), String> {
    let count = array.repeats.len();
    if count == 0 || array.left_flanks.len() != count || array.spacers.len() != count {
        return Err("array must have matching, nonempty repeat, flank, and spacer rows".into());
    }
    let width = array.repeats[0].len();
    if width == 0 || array.repeats.iter().any(|repeat| repeat.len() != width) {
        return Err("array repeats must have the same nonzero aligned width".into());
    }

    let draft_min = min_array.saturating_sub(2).max(3);
    if partial_repeat(array, 0) {
        if count <= draft_min {
            return Ok(());
        }
        let offset =
            array.repeats[0].iter().filter(|&&c| c != b'-').count() + array.spacers[0].len();
        array.pos = array
            .pos
            .checked_add(offset)
            .ok_or("array position overflow")?;
        array.left_flanks.remove(0);
        array.repeats.remove(0);
        array.spacers.remove(0);
    }
    let count = array.repeats.len();
    if partial_repeat(array, count - 1) {
        if count <= draft_min {
            return Ok(());
        }
        array.left_flanks.pop();
        array.repeats.pop();
        array.spacers.pop();
    }

    fix_terminal_gaps(array)?;
    while left_conserved(array) {
        let (_, local) = sequence
            .global_to_local(array.pos)
            .ok_or("array position is outside a contig")?;
        if local == 0 {
            break;
        }
        extend_left(array, sequence)?;
    }
    while right_conserved(array) {
        for (repeat, spacer) in array.repeats.iter_mut().zip(&mut array.spacers) {
            repeat.push(spacer.remove(0));
        }
    }
    while !column_conserved(array, 0)? {
        trim_left(array)?;
    }
    loop {
        let width = array.repeats[0].len();
        if width == 0 {
            return Err("no conserved repeat columns remain".into());
        }
        if column_conserved(array, width - 1)? {
            break;
        }
        for (repeat, spacer) in array.repeats.iter_mut().zip(&mut array.spacers) {
            let c = repeat.pop().unwrap();
            if c != b'-' {
                spacer.insert(0, c);
            }
        }
    }
    fix_terminal_gaps(array)
}

fn gap_count(repeat: &[u8]) -> usize {
    repeat.iter().filter(|&&c| c == b'-').count()
}

fn partial_repeat(array: &Array, row: usize) -> bool {
    let sum: usize = array.repeats.iter().map(|repeat| gap_count(repeat)).sum();
    let average = (sum as f64 / array.repeats.len() as f64 + 0.5) as usize;
    let count = gap_count(&array.repeats[row]);
    count >= 3 && count > average * 2
}

fn column_conserved(array: &Array, column: usize) -> Result<bool, String> {
    let mut counts = [0usize; 4];
    for repeat in &array.repeats {
        let c = repeat
            .get(column)
            .ok_or("no conserved repeat columns remain")?;
        let letter = match c.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => continue,
        };
        counts[letter] += 1;
    }
    let maximum = *counts.iter().max().unwrap();
    let count = array.repeats.len();
    Ok(count - maximum <= 1 || maximum as f64 / count as f64 >= 0.90)
}

fn left_conserved(array: &Array) -> bool {
    match array.left_flanks[0].last() {
        Some(&c) if c != 0 && c != b'-' => array
            .left_flanks
            .iter()
            .all(|flank| flank.last() == Some(&c)),
        _ => false,
    }
}

fn right_conserved(array: &Array) -> bool {
    match array.spacers[0].first() {
        Some(&c) if c != 0 => array
            .spacers
            .iter()
            .all(|spacer| spacer.first() == Some(&c)),
        _ => false,
    }
}

fn extend_left(array: &mut Array, sequence: &Sequence) -> Result<(), String> {
    let (_, local) = sequence
        .global_to_local(array.pos)
        .ok_or("array position is outside a contig")?;
    if local == 0 {
        return Err("left extension precedes the contig".into());
    }
    let pos = array.pos - 1;
    let context = if local > 10 {
        sequence.bases[pos - 10]
    } else {
        b' '
    };
    array.pos = pos;
    let count = array.repeats.len();
    for row in 0..count {
        let flank = &array.left_flanks[row];
        let mut c = *flank.last().unwrap();
        if c == b' ' {
            c = b'-';
        }
        array.repeats[row].insert(0, c);
        if row != count - 1 {
            let spacer = &mut array.spacers[row];
            spacer.truncate(spacer.len().saturating_sub(1));
        }
        if row > 0 {
            let previous = &array.spacers[row - 1];
            let mut flank = vec![b' '; 10usize.saturating_sub(previous.len())];
            flank.extend_from_slice(&previous[previous.len().saturating_sub(10)..]);
            array.left_flanks[row] = flank;
        } else {
            let flank = &mut array.left_flanks[0];
            flank.truncate(9);
            flank.insert(0, context);
        }
    }
    Ok(())
}

fn trim_left(array: &mut Array) -> Result<(), String> {
    if array
        .repeats
        .iter()
        .zip(&array.left_flanks)
        .any(|(repeat, flank)| repeat[0] != b'-' && flank.is_empty())
    {
        return Err("cannot move a repeat base into an empty left flank".into());
    }
    if array.repeats[0][0] != b'-' {
        array.pos = array.pos.checked_add(1).ok_or("array position overflow")?;
    }
    for row in 0..array.repeats.len() {
        let c = array.repeats[row].remove(0);
        if c != b'-' {
            array.left_flanks[row].remove(0);
            array.left_flanks[row].push(c);
            if row > 0 {
                array.spacers[row - 1].push(c);
            }
        }
    }
    Ok(())
}

fn fix_terminal_gaps(array: &mut Array) -> Result<(), String> {
    let count = array.repeats.len();
    for row in 0..count {
        let gap = array.repeats[row]
            .iter()
            .take_while(|&&c| c == b'-')
            .count();
        if (1..=3).contains(&gap)
            && row > 0
            && array.spacers[row - 1].len() >= 10 + gap
            && array.left_flanks[row].len() > gap
        {
            let flank = &array.left_flanks[row];
            let start = flank.len() - gap;
            array.repeats[row][..gap].copy_from_slice(&flank[start..]);
            let previous = &mut array.spacers[row - 1];
            previous.truncate(previous.len() - gap);
            array.left_flanks[row] = previous[previous.len() - 10..].to_vec();
        }
        let gap = array.repeats[row]
            .iter()
            .rev()
            .take_while(|&&c| c == b'-')
            .count();
        if gap == 0 || gap > 2 {
            continue;
        }
        let spacer = &mut array.spacers[row];
        if spacer.len() < gap {
            if row != count - 1 {
                return Err("SpacerLength < GapLength".into());
            }
            continue;
        }
        let repeat = &mut array.repeats[row];
        let start = repeat.len() - gap;
        repeat[start..].copy_from_slice(&spacer[..gap]);
        spacer.drain(..gap);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sequence() -> Sequence {
        Sequence::from_reader(format!(">test\n{}\n", "N".repeat(256)).as_bytes()).unwrap()
    }

    fn array(repeats: &[&[u8]]) -> Array {
        let count = repeats.len();
        Array {
            id: 7,
            pos: 100,
            consensus: b"old".to_vec(),
            aligned_consensus: b"o-ld".to_vec(),
            left_flanks: (0..count).map(|row| vec![b"ACGT"[row % 4]; 10]).collect(),
            repeats: repeats.iter().map(|repeat| repeat.to_vec()).collect(),
            spacers: (0..count).map(|row| vec![b"CGTA"[row % 4]; 12]).collect(),
            cluster_rev_comp: true,
        }
    }

    #[test]
    fn preserves_metadata_and_consensus() {
        let mut a = array(&[b"ACGT", b"acgt", b"ACGT"]);
        let before = a.repeats.clone();
        fix_bounds(&mut a, &sequence(), 3).unwrap();
        assert_eq!(a.repeats, before);
        assert_eq!(a.pos, 100);
        assert_eq!(a.id, 7);
        assert!(a.cluster_rev_comp);
        assert_eq!(a.consensus, b"old");
        assert_eq!(a.aligned_consensus, b"o-ld");
    }

    #[test]
    fn partial_repeat_at_draft_minimum_returns_before_any_fix() {
        let mut a = array(&[b"----AC", b"TTTTAC", b"TTTTAC", b"TTTTAC"]);
        let before = a.repeats.clone();
        fix_bounds(&mut a, &sequence(), 6).unwrap();
        assert_eq!(a.repeats, before);
        assert_eq!(a.pos, 100);
    }

    #[test]
    fn initial_deletion_precedes_final_minimum_check() {
        let mut a = array(&[
            b"-----ACGT",
            b"TTTTTACGT",
            b"TTTTTACGT",
            b"TTTTTACGT",
            b"TTTT-----",
        ]);
        let next_flank = a.left_flanks[1].clone();
        fix_bounds(&mut a, &sequence(), 6).unwrap();
        assert_eq!(a.pos, 116);
        assert_eq!(a.repeats.len(), 4);
        assert_eq!(a.repeats.last().unwrap(), b"TTTT-----");
        assert_eq!(a.left_flanks[0], next_flank);
    }

    #[test]
    fn removes_final_partial_without_merging_its_right_flank() {
        let mut a = array(&[b"ACGTTTTT", b"ACGTTTTT", b"ACGTTTTT", b"ACGT----"]);
        let spacer = a.spacers[2].clone();
        fix_bounds(&mut a, &sequence(), 3).unwrap();
        assert_eq!(a.repeats.len(), 3);
        assert_eq!(a.pos, 100);
        assert_eq!(a.spacers[2], spacer);
    }

    #[test]
    fn terminal_gap_pass_is_left_then_right_for_each_row() {
        let mut a = array(&[b"AACGTT", b"--CGTT", b"AACG--"]);
        a.spacers[0] = b"CCCCCCCCCCAA".to_vec();
        a.left_flanks[1] = b"CCCCCCCCAA".to_vec();
        a.spacers[2] = b"TTCGG".to_vec();
        fix_bounds(&mut a, &sequence(), 3).unwrap();
        assert_eq!(a.repeats, vec![b"AACGTT".to_vec(); 3]);
        assert_eq!(a.spacers[0], b"CCCCCCCCCC");
        assert_eq!(a.left_flanks[1], b"CCCCCCCCCC");
        assert_eq!(a.spacers[2], b"CGG");
    }

    #[test]
    fn extends_left_before_right_and_rebuilds_flanks() {
        let mut a = array(&[b"AC", b"AC", b"AC"]);
        a.pos = 5;
        a.left_flanks = vec![b"         G".to_vec(); 3];
        a.spacers = vec![b"TACG".to_vec(), b"TCAG".to_vec(), b"TG".to_vec()];
        fix_bounds(&mut a, &sequence(), 3).unwrap();
        assert_eq!(a.pos, 4);
        assert_eq!(a.repeats, vec![b"GACT".to_vec(); 3]);
        assert_eq!(a.left_flanks[0], b"          ");
        assert_eq!(a.left_flanks[1], b"       TAC");
        assert_eq!(a.left_flanks[2], b"       TCA");
        assert_eq!(
            a.spacers,
            vec![b"AC".to_vec(), b"CA".to_vec(), b"G".to_vec()]
        );
    }

    #[test]
    fn trims_left_then_right_without_moving_gaps_as_bases() {
        let mut a = array(&[b"-ACGA", b"CACGC", b"GACGG"]);
        a.spacers[0] = b"C".to_vec();
        fix_bounds(&mut a, &sequence(), 3).unwrap();
        assert_eq!(a.pos, 100);
        assert_eq!(a.repeats, vec![b"ACG".to_vec(); 3]);
        assert_eq!(a.left_flanks[0], b"AAAAAAAAAA");
        assert_eq!(a.left_flanks[1], b"CCCCCCCCCC");
        assert_eq!(a.spacers[0], b"ACC");
        assert_eq!(a.spacers[1], b"CGGGGGGGGGGGGG");
        assert_eq!(a.spacers[2], b"GTTTTTTTTTTTT");
    }

    #[test]
    fn second_gap_pass_fills_gaps_exposed_by_trimming() {
        let mut a = array(&[b"AACGTA", b"CACG-C", b"GACGTG"]);
        a.spacers[1] = b"TCCCCCCCCCCC".to_vec();
        fix_bounds(&mut a, &sequence(), 3).unwrap();
        assert_eq!(a.pos, 101);
        assert_eq!(
            a.repeats,
            vec![b"ACGT".to_vec(), b"ACGC".to_vec(), b"ACGT".to_vec()]
        );
        assert_eq!(a.spacers[1], b"TCCCCCCCCCCCG");
    }

    #[test]
    fn right_gap_shortage_is_fatal_only_for_internal_rows() {
        let mut a = array(&[b"AC-", b"ACT", b"ACT"]);
        a.spacers[0].clear();
        assert_eq!(
            fix_bounds(&mut a, &sequence(), 3).unwrap_err(),
            "SpacerLength < GapLength"
        );
        let mut a = array(&[b"ACT", b"ACT", b"AC-"]);
        a.spacers[2].clear();
        fix_bounds(&mut a, &sequence(), 3).unwrap();
        assert_eq!(a.repeats[2], b"AC-");
    }

    #[test]
    fn conservation_uses_one_difference_or_ninety_percent() {
        let mut a = array(&[b"AC", b"AC", b"GC"]);
        assert!(column_conserved(&a, 0).unwrap());
        a.repeats = vec![b"AC".to_vec(); 20];
        a.repeats[18][0] = b'N';
        a.repeats[19][0] = b'-';
        assert!(column_conserved(&a, 0).unwrap());
        a.repeats[17][0] = b'G';
        assert!(!column_conserved(&a, 0).unwrap());
    }

    #[test]
    fn extension_is_case_sensitive_and_right_allows_gaps() {
        let mut a = array(&[b"AC", b"AC", b"AC"]);
        a.left_flanks = vec![
            b"AAAAAAAAAG".to_vec(),
            b"CCCCCCCCCG".to_vec(),
            b"GGGGGGGGGg".to_vec(),
        ];
        a.spacers = vec![b"-C".to_vec(), b"-G".to_vec(), b"-T".to_vec()];
        assert!(!left_conserved(&a));
        assert!(right_conserved(&a));
        fix_bounds(&mut a, &sequence(), 3).unwrap();
        assert_eq!(a.repeats, vec![b"AC".to_vec(); 3]);
        assert_eq!(a.spacers, vec![b"C".to_vec(), b"G".to_vec(), b"T".to_vec()]);
    }

    #[test]
    fn invalid_alignment_is_an_error() {
        let mut a = array(&[b"AC", b"AC", b"AC"]);
        a.repeats[1].pop();
        assert!(
            fix_bounds(&mut a, &sequence(), 3)
                .unwrap_err()
                .contains("aligned width")
        );
    }

    #[test]
    fn extends_beyond_flank_width_preserving_genomic_case_and_ambiguity() {
        let sequence =
            Sequence::from_reader(&b">test\nCaaGGnTCCAAGGACcaaGGnTCCAAGGACcaaGGnTCCAAGGACT\n"[..])
                .unwrap();
        let mut a = array(&[b"AC", b"AC", b"AC"]);
        a.pos = 13;
        a.left_flanks = vec![b"GGnTCCAAGG".to_vec(); 3];
        a.spacers = vec![
            b"caaGGnTCCAAGG".to_vec(),
            b"caaGGnTCCAAGG".to_vec(),
            b"T".to_vec(),
        ];
        fix_bounds(&mut a, &sequence, 3).unwrap();
        assert_eq!(a.pos, 1);
        assert_eq!(a.repeats, vec![b"aaGGnTCCAAGGAC".to_vec(); 3]);
        assert_eq!(
            a.left_flanks,
            vec![
                b"         C".to_vec(),
                b"         c".to_vec(),
                b"         c".to_vec()
            ]
        );
        assert_eq!(a.spacers, vec![b"c".to_vec(), b"c".to_vec(), b"T".to_vec()]);
    }

    #[test]
    fn contig_start_uses_blanks_without_reading_previous_contig_or_padding() {
        for fasta in [
            &b">test\nGACGACGACT\n"[..],
            &b">previous\nCCCCCCCCCCCCCCCCCCCC\n>test\nGACGACGACT\n"[..],
        ] {
            let sequence = Sequence::from_reader(fasta).unwrap();
            let start = sequence.contigs.last().unwrap().start;
            let mut a = array(&[b"AC", b"AC", b"AC"]);
            a.pos = start + 1;
            a.left_flanks = vec![b"         G".to_vec(); 3];
            a.spacers = vec![b"G".to_vec(), b"G".to_vec(), b"T".to_vec()];
            fix_bounds(&mut a, &sequence, 3).unwrap();
            assert_eq!(a.pos, start);
            assert_eq!(a.repeats, vec![b"GAC".to_vec(); 3]);
            assert_eq!(a.left_flanks, vec![b"          ".to_vec(); 3]);
            assert_eq!(a.spacers, vec![Vec::new(), Vec::new(), b"T".to_vec()]);
            fix_bounds(&mut a, &sequence, 3).unwrap();
            assert_eq!(a.pos, start);
            assert_eq!(a.repeats, vec![b"GAC".to_vec(); 3]);
        }
    }
}

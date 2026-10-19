use crate::{array::Array, seq::Sequence};

pub(crate) fn trace_rows(array: &Array, phase: &str) {
    if !crate::trace::enabled() {
        return;
    }
    for (row_index, element) in array.elements.iter().enumerate() {
        crate::trace::emit(
            "14-bounds",
            &[
                array.id.to_string(),
                phase.into(),
                array.pos.to_string(),
                row_index.to_string(),
                String::from_utf8_lossy(&element.left_flank).into_owned(),
                String::from_utf8_lossy(&element.repeat).into_owned(),
                String::from_utf8_lossy(&element.spacer).into_owned(),
            ],
        );
    }
}

pub fn fix_bounds(array: &mut Array, sequence: &Sequence, min_array: usize) -> Result<(), String> {
    let count = array.elements.len();
    if count == 0 {
        return Err("array must have matching, nonempty repeat, flank, and spacer rows".into());
    }
    let width = array.elements[0].repeat.len();
    if width == 0 || array.elements.iter().any(|e| e.repeat.len() != width) {
        return Err("array repeats must have the same nonzero aligned width".into());
    }

    let draft_min = min_array.saturating_sub(2).max(3);
    if partial_repeat(array, 0) {
        if count <= draft_min {
            return Ok(());
        }
        let offset = array.elements[0]
            .repeat
            .iter()
            .filter(|&&c| c != b'-')
            .count()
            + array.elements[0].spacer.len();
        array.pos = array
            .pos
            .checked_add(offset)
            .ok_or("array position overflow")?;
        array.elements.remove(0);
    }
    let count = array.elements.len();
    if partial_repeat(array, count - 1) {
        if count <= draft_min {
            return Ok(());
        }
        array.elements.pop();
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
        for element in &mut array.elements {
            element.repeat.push(element.spacer.remove(0));
        }
    }
    while !column_conserved(array, 0)? {
        trim_left(array)?;
    }
    loop {
        let width = array.elements[0].repeat.len();
        if width == 0 {
            return Err("no conserved repeat columns remain".into());
        }
        if column_conserved(array, width - 1)? {
            break;
        }
        for element in &mut array.elements {
            let c = element.repeat.pop().unwrap();
            if c != b'-' {
                element.spacer.insert(0, c);
            }
        }
    }
    fix_terminal_gaps(array)
}

fn gap_count(repeat: &[u8]) -> usize {
    repeat.iter().filter(|&&c| c == b'-').count()
}

fn partial_repeat(array: &Array, row: usize) -> bool {
    let sum: usize = array.elements.iter().map(|e| gap_count(&e.repeat)).sum();
    let average = (sum as f64 / array.elements.len() as f64 + 0.5) as usize;
    let count = gap_count(&array.elements[row].repeat);
    count >= 3 && count > average * 2
}

fn column_conserved(array: &Array, column: usize) -> Result<bool, String> {
    let mut counts = [0usize; 4];
    for repeat in array.elements.iter().map(|e| &e.repeat) {
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
    let count = array.elements.len();
    Ok(count - maximum <= 1 || maximum as f64 / count as f64 >= 0.90)
}

fn left_conserved(array: &Array) -> bool {
    match array.elements[0].left_flank.last() {
        Some(&c) if c != 0 && c != b'-' => array
            .elements
            .iter()
            .all(|e| e.left_flank.last() == Some(&c)),
        _ => false,
    }
}

fn right_conserved(array: &Array) -> bool {
    match array.elements[0].spacer.first() {
        Some(&c) if c != 0 => array.elements.iter().all(|e| e.spacer.first() == Some(&c)),
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
    let count = array.elements.len();
    for row in 0..count {
        let mut c = *array.elements[row].left_flank.last().unwrap();
        if c == b' ' {
            c = b'-';
        }
        array.elements[row].repeat.insert(0, c);
        if row != count - 1 {
            let spacer = &mut array.elements[row].spacer;
            spacer.truncate(spacer.len().saturating_sub(1));
        }
        if row > 0 {
            let previous = &array.elements[row - 1].spacer;
            let mut flank = vec![b' '; 10usize.saturating_sub(previous.len())];
            flank.extend_from_slice(&previous[previous.len().saturating_sub(10)..]);
            array.elements[row].left_flank = flank;
        } else {
            let flank = &mut array.elements[0].left_flank;
            flank.truncate(9);
            flank.insert(0, context);
        }
    }
    Ok(())
}

fn trim_left(array: &mut Array) -> Result<(), String> {
    if array
        .elements
        .iter()
        .any(|e| e.repeat[0] != b'-' && e.left_flank.is_empty())
    {
        return Err("cannot move a repeat base into an empty left flank".into());
    }
    if array.elements[0].repeat[0] != b'-' {
        array.pos = array.pos.checked_add(1).ok_or("array position overflow")?;
    }
    for row in 0..array.elements.len() {
        let c = array.elements[row].repeat.remove(0);
        if c != b'-' {
            array.elements[row].left_flank.remove(0);
            array.elements[row].left_flank.push(c);
            if row > 0 {
                array.elements[row - 1].spacer.push(c);
            }
        }
    }
    Ok(())
}

fn fix_terminal_gaps(array: &mut Array) -> Result<(), String> {
    let count = array.elements.len();
    for row in 0..count {
        let gap = array.elements[row]
            .repeat
            .iter()
            .take_while(|&&c| c == b'-')
            .count();
        if (1..=3).contains(&gap)
            && row > 0
            && array.elements[row - 1].spacer.len() >= 10 + gap
            && array.elements[row].left_flank.len() > gap
        {
            let start = array.elements[row].left_flank.len() - gap;
            let flank = array.elements[row].left_flank[start..].to_vec();
            array.elements[row].repeat[..gap].copy_from_slice(&flank);
            let previous = &mut array.elements[row - 1].spacer;
            previous.truncate(previous.len() - gap);
            array.elements[row].left_flank = previous[previous.len() - 10..].to_vec();
        }
        let gap = array.elements[row]
            .repeat
            .iter()
            .rev()
            .take_while(|&&c| c == b'-')
            .count();
        if gap == 0 || gap > 2 {
            continue;
        }
        let element = &mut array.elements[row];
        if element.spacer.len() < gap {
            if row != count - 1 {
                return Err("SpacerLength < GapLength".into());
            }
            continue;
        }
        let start = element.repeat.len() - gap;
        element.repeat[start..].copy_from_slice(&element.spacer[..gap]);
        element.spacer.drain(..gap);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::array::ArrayElement;

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
            elements: (0..count)
                .map(|row| ArrayElement {
                    left_flank: vec![b"ACGT"[row % 4]; 10],
                    repeat: repeats[row].to_vec(),
                    spacer: vec![b"CGTA"[row % 4]; 12],
                })
                .collect(),
            cluster_rev_comp: true,
        }
    }

    fn repeat_rows(a: &Array) -> Vec<Vec<u8>> {
        a.elements.iter().map(|e| e.repeat.clone()).collect()
    }

    fn spacer_rows(a: &Array) -> Vec<Vec<u8>> {
        a.elements.iter().map(|e| e.spacer.clone()).collect()
    }

    fn flank_rows(a: &Array) -> Vec<Vec<u8>> {
        a.elements.iter().map(|e| e.left_flank.clone()).collect()
    }

    fn set_left_flanks(a: &mut Array, flanks: Vec<Vec<u8>>) {
        for (element, flank) in a.elements.iter_mut().zip(flanks) {
            element.left_flank = flank;
        }
    }

    fn set_spacers(a: &mut Array, spacers: Vec<Vec<u8>>) {
        for (element, spacer) in a.elements.iter_mut().zip(spacers) {
            element.spacer = spacer;
        }
    }

    #[test]
    fn preserves_metadata_and_consensus() {
        let mut a = array(&[b"ACGT", b"acgt", b"ACGT"]);
        let before = repeat_rows(&a);
        fix_bounds(&mut a, &sequence(), 3).unwrap();
        assert_eq!(repeat_rows(&a), before);
        assert_eq!(a.pos, 100);
        assert_eq!(a.id, 7);
        assert!(a.cluster_rev_comp);
        assert_eq!(a.consensus, b"old");
        assert_eq!(a.aligned_consensus, b"o-ld");
    }

    #[test]
    fn partial_repeat_at_draft_minimum_returns_before_any_fix() {
        let mut a = array(&[b"----AC", b"TTTTAC", b"TTTTAC", b"TTTTAC"]);
        let before = repeat_rows(&a);
        fix_bounds(&mut a, &sequence(), 6).unwrap();
        assert_eq!(repeat_rows(&a), before);
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
        let next_flank = a.elements[1].left_flank.clone();
        fix_bounds(&mut a, &sequence(), 6).unwrap();
        assert_eq!(a.pos, 116);
        assert_eq!(a.elements.len(), 4);
        assert_eq!(a.elements.last().unwrap().repeat, b"TTTT-----");
        assert_eq!(a.elements[0].left_flank, next_flank);
    }

    #[test]
    fn removes_final_partial_without_merging_its_right_flank() {
        let mut a = array(&[b"ACGTTTTT", b"ACGTTTTT", b"ACGTTTTT", b"ACGT----"]);
        let spacer = a.elements[2].spacer.clone();
        fix_bounds(&mut a, &sequence(), 3).unwrap();
        assert_eq!(a.elements.len(), 3);
        assert_eq!(a.pos, 100);
        assert_eq!(a.elements[2].spacer, spacer);
    }

    #[test]
    fn terminal_gap_pass_is_left_then_right_for_each_row() {
        let mut a = array(&[b"AACGTT", b"--CGTT", b"AACG--"]);
        a.elements[0].spacer = b"CCCCCCCCCCAA".to_vec();
        a.elements[1].left_flank = b"CCCCCCCCAA".to_vec();
        a.elements[2].spacer = b"TTCGG".to_vec();
        fix_bounds(&mut a, &sequence(), 3).unwrap();
        assert_eq!(repeat_rows(&a), vec![b"AACGTT".to_vec(); 3]);
        assert_eq!(a.elements[0].spacer, b"CCCCCCCCCC");
        assert_eq!(a.elements[1].left_flank, b"CCCCCCCCCC");
        assert_eq!(a.elements[2].spacer, b"CGG");
    }

    #[test]
    fn extends_left_before_right_and_rebuilds_flanks() {
        let mut a = array(&[b"AC", b"AC", b"AC"]);
        a.pos = 5;
        set_left_flanks(&mut a, vec![b"         G".to_vec(); 3]);
        set_spacers(
            &mut a,
            vec![b"TACG".to_vec(), b"TCAG".to_vec(), b"TG".to_vec()],
        );
        fix_bounds(&mut a, &sequence(), 3).unwrap();
        assert_eq!(a.pos, 4);
        assert_eq!(repeat_rows(&a), vec![b"GACT".to_vec(); 3]);
        assert_eq!(a.elements[0].left_flank, b"          ");
        assert_eq!(a.elements[1].left_flank, b"       TAC");
        assert_eq!(a.elements[2].left_flank, b"       TCA");
        assert_eq!(
            spacer_rows(&a),
            vec![b"AC".to_vec(), b"CA".to_vec(), b"G".to_vec()]
        );
    }

    #[test]
    fn trims_left_then_right_without_moving_gaps_as_bases() {
        let mut a = array(&[b"-ACGA", b"CACGC", b"GACGG"]);
        a.elements[0].spacer = b"C".to_vec();
        fix_bounds(&mut a, &sequence(), 3).unwrap();
        assert_eq!(a.pos, 100);
        assert_eq!(repeat_rows(&a), vec![b"ACG".to_vec(); 3]);
        assert_eq!(a.elements[0].left_flank, b"AAAAAAAAAA");
        assert_eq!(a.elements[1].left_flank, b"CCCCCCCCCC");
        assert_eq!(a.elements[0].spacer, b"ACC");
        assert_eq!(a.elements[1].spacer, b"CGGGGGGGGGGGGG");
        assert_eq!(a.elements[2].spacer, b"GTTTTTTTTTTTT");
    }

    #[test]
    fn second_gap_pass_fills_gaps_exposed_by_trimming() {
        let mut a = array(&[b"AACGTA", b"CACG-C", b"GACGTG"]);
        a.elements[1].spacer = b"TCCCCCCCCCCC".to_vec();
        fix_bounds(&mut a, &sequence(), 3).unwrap();
        assert_eq!(a.pos, 101);
        assert_eq!(
            repeat_rows(&a),
            vec![b"ACGT".to_vec(), b"ACGC".to_vec(), b"ACGT".to_vec()]
        );
        assert_eq!(a.elements[1].spacer, b"TCCCCCCCCCCCG");
    }

    #[test]
    fn right_gap_shortage_is_fatal_only_for_internal_rows() {
        let mut a = array(&[b"AC-", b"ACT", b"ACT"]);
        a.elements[0].spacer.clear();
        assert_eq!(
            fix_bounds(&mut a, &sequence(), 3).unwrap_err(),
            "SpacerLength < GapLength"
        );
        let mut a = array(&[b"ACT", b"ACT", b"AC-"]);
        a.elements[2].spacer.clear();
        fix_bounds(&mut a, &sequence(), 3).unwrap();
        assert_eq!(a.elements[2].repeat, b"AC-");
    }

    #[test]
    fn conservation_uses_one_difference_or_ninety_percent() {
        let mut a = array(&[b"AC", b"AC", b"GC"]);
        assert!(column_conserved(&a, 0).unwrap());
        a.elements = (0..20)
            .map(|_| ArrayElement {
                repeat: b"AC".to_vec(),
                ..ArrayElement::default()
            })
            .collect();
        a.elements[18].repeat[0] = b'N';
        a.elements[19].repeat[0] = b'-';
        assert!(column_conserved(&a, 0).unwrap());
        a.elements[17].repeat[0] = b'G';
        assert!(!column_conserved(&a, 0).unwrap());
    }

    #[test]
    fn extension_is_case_sensitive_and_right_allows_gaps() {
        let mut a = array(&[b"AC", b"AC", b"AC"]);
        set_left_flanks(
            &mut a,
            vec![
                b"AAAAAAAAAG".to_vec(),
                b"CCCCCCCCCG".to_vec(),
                b"GGGGGGGGGg".to_vec(),
            ],
        );
        set_spacers(&mut a, vec![b"-C".to_vec(), b"-G".to_vec(), b"-T".to_vec()]);
        assert!(!left_conserved(&a));
        assert!(right_conserved(&a));
        fix_bounds(&mut a, &sequence(), 3).unwrap();
        assert_eq!(repeat_rows(&a), vec![b"AC".to_vec(); 3]);
        assert_eq!(
            spacer_rows(&a),
            vec![b"C".to_vec(), b"G".to_vec(), b"T".to_vec()]
        );
    }

    #[test]
    fn invalid_alignment_is_an_error() {
        let mut a = array(&[b"AC", b"AC", b"AC"]);
        a.elements[1].repeat.pop();
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
        set_left_flanks(&mut a, vec![b"GGnTCCAAGG".to_vec(); 3]);
        set_spacers(
            &mut a,
            vec![
                b"caaGGnTCCAAGG".to_vec(),
                b"caaGGnTCCAAGG".to_vec(),
                b"T".to_vec(),
            ],
        );
        fix_bounds(&mut a, &sequence, 3).unwrap();
        assert_eq!(a.pos, 1);
        assert_eq!(repeat_rows(&a), vec![b"aaGGnTCCAAGGAC".to_vec(); 3]);
        assert_eq!(
            flank_rows(&a),
            vec![
                b"         C".to_vec(),
                b"         c".to_vec(),
                b"         c".to_vec()
            ]
        );
        assert_eq!(
            spacer_rows(&a),
            vec![b"c".to_vec(), b"c".to_vec(), b"T".to_vec()]
        );
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
            set_left_flanks(&mut a, vec![b"         G".to_vec(); 3]);
            set_spacers(&mut a, vec![b"G".to_vec(), b"G".to_vec(), b"T".to_vec()]);
            fix_bounds(&mut a, &sequence, 3).unwrap();
            assert_eq!(a.pos, start);
            assert_eq!(repeat_rows(&a), vec![b"GAC".to_vec(); 3]);
            assert_eq!(flank_rows(&a), vec![b"          ".to_vec(); 3]);
            assert_eq!(spacer_rows(&a), vec![Vec::new(), Vec::new(), b"T".to_vec()]);
            fix_bounds(&mut a, &sequence, 3).unwrap();
            assert_eq!(a.pos, start);
            assert_eq!(repeat_rows(&a), vec![b"GAC".to_vec(); 3]);
        }
    }
}

use crate::{array::Array, report, seq::Sequence};
use std::fmt::Write;

fn similar(a: &[u8], b: &[u8]) -> Result<bool, String> {
    for candidate in [a.to_vec(), report::reverse_complement(a)] {
        let aligned = crate::msa::align(&[candidate, b.to_vec()])?;
        let same = aligned.rows[0]
            .iter()
            .zip(&aligned.rows[1])
            .filter(|(a, b)| **a != b'-' && a.eq_ignore_ascii_case(b))
            .count();
        if same as f64 / a.len().max(b.len()) as f64 > 0.899 {
            return Ok(true);
        }
    }
    Ok(false)
}

pub fn consensus_fasta(
    sequence: &Sequence,
    arrays: &[Array],
    trim: bool,
) -> Result<String, String> {
    for a in arrays {
        report::location(sequence, a)?;
    }
    let mut sorted: Vec<_> = arrays.iter().collect();
    sorted.sort_by_key(|a| a.pos);
    let mut out = String::new();
    for (index, a) in sorted.iter().enumerate() {
        if trim {
            let mut discard = false;
            for previous in &sorted[..index] {
                if similar(&a.consensus, &previous.consensus)? {
                    discard = true;
                    break;
                }
            }
            if discard {
                continue;
            }
        }
        let (ci, local, _) = report::location(sequence, a)?;
        let label = &sequence.contigs[ci].label;
        if label.contains(['\r', '\n']) {
            return Err("FASTA label contains a newline".into());
        }
        write!(out, ">{label}[Array{};Pos={}]", a.id, local + 1).unwrap();
        for chunk in a.consensus.chunks(60) {
            out.push('\n');
            out.push_str(std::str::from_utf8(chunk).map_err(|e| e.to_string())?);
        }
        out.push('\n');
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::Contig;

    #[test]
    fn headers_wrapping_sorting_and_trim() {
        let sequence = Sequence {
            bases: vec![b'A'; 300],
            contigs: vec![Contig {
                label: "scf1 description".into(),
                start: 0,
                len: 300,
            }],
        };
        let make = |id, pos| Array {
            id,
            pos,
            consensus: vec![b'A'; 61],
            aligned_consensus: vec![b'A'; 61],
            elements: vec![crate::array::ArrayElement {
                left_flank: vec![],
                repeat: vec![b'A'; 61],
                spacer: b"TT".to_vec(),
            }],
            cluster_rev_comp: false,
        };
        let arrays = [make(2, 100), make(1, 0)];
        let first = format!(">scf1 description[Array1;Pos=1]\n{}\nA\n", "A".repeat(60));
        let second = format!(">scf1 description[Array2;Pos=101]\n{}\nA\n", "A".repeat(60));
        assert_eq!(
            consensus_fasta(&sequence, &arrays, false).unwrap(),
            format!("{first}{second}")
        );
        assert_eq!(consensus_fasta(&sequence, &arrays, true).unwrap(), first);
        assert!(similar(b"AAAACCCCGGGG", b"CCCCGGGGTTTT").unwrap());
        assert_eq!(consensus_fasta(&sequence, &[], true).unwrap(), "");
    }
}

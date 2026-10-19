use crate::{array::Array, seq::Sequence};
use std::fmt::Write;

const INFO: &str = concat!(
    "Help on reading this report\n",
    "===========================\n",
    "\n",
    "This report has three sections: Detailed, Summary by Similarity\n",
    "and Summary by Position.\n",
    "\n",
    "The detailed section shows each repeat in each putative\n",
    "CRISPR array.\n",
    "\n",
    "The summary sections give one line for each array.\n",
    "\n",
    "An 'array' is a contiguous sequence of CRISPR repeats\n",
    "looking like this:\n",
    "\n",
    "    REPEAT Spacer REPEAT Spacer REPEAT ... Spacer REPEAT\n",
    "\n",
    "Within one array, repeats have high similarity and spacers\n",
    "are, roughly speaking, unique within a window around the array.\n",
    "In a given array, each repeat has a similar length, and each\n",
    "spacer has a similar length. With default parameters, the\n",
    "algorithm allows a fair amount of variability in order to\n",
    "maximize sensitivity. This may allow identification of\n",
    "inactive (\"fossil\") arrays, and may in rare cases also\n",
    "induce false positives due to other classes of repeats\n",
    "such as microsatellites, LTRs and arrays of RNA genes.\n",
    "\n",
    "\n",
    "Columns in the detailed section are:\n",
    "\n",
    "  Pos               Sequence position, starting at 1 for the first base.\n",
    "  Repeat            Length of the repeat.\n",
    "  %id               Identity with the consensus sequence.\n",
    "  Spacer            Length of spacer to the right of this repeat.\n",
    "  Left flank        10 bases to the left of this repeat.\n",
    "  Repeat            Sequence of this repeat.\n",
    "                      Dots indicate positions where this repeat\n",
    "                      agrees with the consensus sequence below.\n",
    "  Spacer            Sequence of spacer to the right of this repeat,\n",
    "                      or 10 bases if this is the last repeat.\n",
    "\n",
    "The left flank sequence duplicates the end of the spacer for the preceding\n",
    "repeat; it is provided to facilitate visual identification of cases\n",
    "where the algorithm does not correctly identify repeat endpoints.\n",
    "\n",
    "At the end of each array there is a sub-heading that gives the average\n",
    "repeat length, average spacer length and consensus sequence.\n",
    "\n",
    "Columns in the summary sections are:\n",
    "\n",
    "  Array             Number 1, 2 ... referring back to the detailed report.\n",
    "  Sequence          FASTA label of the sequence. May be truncated.\n",
    "  From              Start position of array.\n",
    "  To                End position of array.\n",
    "  # copies          Number of repeats in the array.\n",
    "  Repeat            Average repeat length.\n",
    "  Spacer            Average spacer length.\n",
    "  +                 +/-, indicating orientation relative to first array in group.\n",
    "  Distance          Distance from previous array.\n",
    "  Consensus         Consensus sequence.\n",
    "\n",
    "In the Summary by Similarity section, arrays are grouped by similarity of their\n",
    "consensus sequences. If consensus sequences are sufficiently similar, they are\n",
    "aligned to each other to indicate probable relationships between arrays.\n",
    "\n",
    "In the Summary by Position section, arrays are sorted by position within the\n",
    "input sequence file.\n",
    "\n",
    "The Distance column facilitates identification of cases where a single\n",
    "array has been reported as two adjacent arrays. In such a case, (a) the\n",
    "consensus sequences will be similar or identical, and (b) the distance\n",
    "will be approximately a small multiple of the repeat length + spacer length.\n",
    "\n",
    "Use the -noinfo option to turn off this help.\n",
    "Use the -help option to get a list of command line options.\n",
    "\n",
);

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum CoordinateMode {
    #[default]
    Corrected,
    Legacy,
}

pub(crate) fn ungapped(s: &[u8]) -> usize {
    s.iter().filter(|&&c| c != b'-').count()
}

fn spacer_average(a: &Array) -> usize {
    a.elements[..a.elements.len() - 1]
        .iter()
        .map(|e| e.spacer.len())
        .sum::<usize>()
        / (a.elements.len() - 1).max(1)
}

fn summary_repeat_average(a: &Array, mode: CoordinateMode) -> usize {
    a.elements
        .iter()
        .map(|e| {
            if mode == CoordinateMode::Legacy {
                e.repeat.len()
            } else {
                ungapped(&e.repeat)
            }
        })
        .sum::<usize>()
        / a.elements.len()
}

/// Location and summary metrics for one array within its sequence.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ArrayLocation {
    /// Index of the contig containing the array.
    pub contig: usize,
    /// Zero-based offset of the array start within that contig.
    pub start: usize,
    /// Total array length in bases (ungapped repeats plus internal spacers).
    pub len: usize,
    /// Number of repeats (copies) in the array.
    pub copies: usize,
    /// Average ungapped repeat length, matching the corrected summary column.
    pub repeat_len: usize,
    /// Average spacer length between repeats, matching the summary column.
    pub spacer_len: usize,
}

/// Locates an array and derives the same coordinates and averages used by the
/// corrected report summary. Validation is delegated to `location` so the
/// helper cannot drift from the report writers.
pub fn locate(sequence: &Sequence, array: &Array) -> Result<ArrayLocation, String> {
    let (contig, start, len) = location(sequence, array)?;
    Ok(ArrayLocation {
        contig,
        start,
        len,
        copies: array.elements.len(),
        repeat_len: summary_repeat_average(array, CoordinateMode::Corrected),
        spacer_len: spacer_average(array),
    })
}

pub(crate) fn location(sequence: &Sequence, a: &Array) -> Result<(usize, usize, usize), String> {
    let (contig, local) = sequence
        .global_to_local(a.pos)
        .ok_or_else(|| format!("array {} starts outside a contig", a.id))?;
    let n = a.elements.len();
    if n == 0
        || a.aligned_consensus.is_empty()
        || a.consensus.is_empty()
        || a.elements
            .iter()
            .any(|e| e.repeat.len() != a.aligned_consensus.len() || ungapped(&e.repeat) == 0)
    {
        return Err(format!("array {} has inconsistent alignment rows", a.id));
    }
    for row in a
        .elements
        .iter()
        .flat_map(|e| [&e.repeat, &e.spacer, &e.left_flank])
        .chain(std::iter::once(&a.consensus))
        .chain(std::iter::once(&a.aligned_consensus))
    {
        if !row.iter().all(|c| c.is_ascii_graphic() || *c == b' ') {
            return Err(format!("array {} contains invalid sequence text", a.id));
        }
    }
    let length = a
        .elements
        .iter()
        .map(|e| ungapped(&e.repeat))
        .chain(a.elements[..n - 1].iter().map(|e| e.spacer.len()))
        .try_fold(0usize, |sum, len| sum.checked_add(len))
        .ok_or_else(|| format!("array {} length overflow", a.id))?;
    if local
        .checked_add(length)
        .is_none_or(|end| end > sequence.contigs[contig].len)
    {
        return Err(format!("array {} crosses a contig boundary", a.id));
    }
    Ok((contig, local, length))
}

fn letters(s: &[u8]) -> String {
    String::from_utf8_lossy(s).into_owned()
}

pub(crate) fn reverse_complement(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|c| match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' | b'U' => b'A',
            b'a' => b't',
            b'c' => b'g',
            b'g' => b'c',
            b't' | b'u' => b'a',
            b'R' => b'Y',
            b'Y' => b'R',
            b'K' => b'M',
            b'M' => b'K',
            b'B' => b'V',
            b'V' => b'B',
            b'D' => b'H',
            b'H' => b'D',
            b'r' => b'y',
            b'y' => b'r',
            b'k' => b'm',
            b'm' => b'k',
            b'b' => b'v',
            b'v' => b'b',
            b'd' => b'h',
            b'h' => b'd',
            _ => *c,
        })
        .collect()
}

fn tuple_counts(s: &[u8]) -> [u8; 256] {
    let mut counts = [0u8; 256];
    for i in 0..s.len().saturating_sub(5) {
        let mut tuple = 0;
        for c in &s[i..i + 4] {
            tuple = tuple * 4
                + match c.to_ascii_uppercase() {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    _ => 3,
                };
        }
        counts[tuple] = counts[tuple].wrapping_add(1);
    }
    counts
}

fn clusters(arrays: &[Array]) -> Vec<Vec<(usize, bool)>> {
    let counts: Vec<_> = arrays
        .iter()
        .flat_map(|a| {
            [
                tuple_counts(&a.consensus),
                tuple_counts(&reverse_complement(&a.consensus)),
            ]
        })
        .collect();
    let totals: Vec<u32> = counts
        .iter()
        .map(|c| c.iter().map(|&x| u32::from(x)).sum::<u32>().max(1))
        .collect();
    let near = |i: usize, j: usize| {
        let common: u32 = counts[i]
            .iter()
            .zip(&counts[j])
            .map(|(&a, &b)| u32::from(a.min(b)))
            .sum();
        i == j || (1.0 - common as f64 / totals[i].min(totals[j]) as f64) as f32 <= 0.5
    };
    let mut done = vec![false; arrays.len()];
    let mut result = Vec::new();
    loop {
        let mut best = None;
        let mut best_count = 0;
        for i in 0..counts.len() {
            if done[i / 2] {
                continue;
            }
            let count = (0..counts.len())
                .filter(|&j| i != j && !done[j / 2] && near(i, j))
                .count();
            if count > best_count {
                best = Some(i);
                best_count = count;
            }
        }
        let Some(best) = best else { break };
        let mut cluster = Vec::new();
        for i in 0..counts.len() {
            if !done[i / 2] && near(best, i) {
                done[i / 2] = true;
                cluster.push((i / 2, i % 2 == 1));
            }
        }
        let flip = cluster[0].1;
        for (_, reverse) in &mut cluster {
            *reverse ^= flip;
        }
        result.push(cluster);
    }
    for (i, &is_done) in done.iter().enumerate() {
        if !is_done {
            result.push(vec![(i, false)]);
        }
    }
    result
}

fn summary_prefix(out: &mut String, sequence: &Sequence, a: &Array, mode: CoordinateMode) {
    let (ci, local, length) = location(sequence, a).unwrap();
    let label: String = sequence.contigs[ci].label.chars().take(16).collect();
    let average = summary_repeat_average(a, mode);
    write!(
        out,
        "{:>5}  {:>16}  {:>10}  {:>10}  {:>8}  {:>6}  {:>6}  ",
        a.id,
        label,
        local + 1,
        length,
        a.elements.len(),
        average,
        spacer_average(a)
    )
    .unwrap();
}

pub fn text(
    sequence: &Sequence,
    arrays: &[Array],
    input_name: &str,
    mode: CoordinateMode,
    no_info: bool,
) -> Result<String, String> {
    for a in arrays {
        location(sequence, a)?;
    }
    let mut out = String::new();
    if !no_info {
        out.push_str(INFO);
    }
    write!(
        out,
        "pilercr v1.06\nBy Robert C. Edgar\n\n{input_name}: {} putative CRISPR arrays found.\n\n",
        arrays.len()
    )
    .unwrap();
    if arrays.is_empty() {
        return Ok(out);
    }
    out.push_str("\n\nDETAIL REPORT\n\n");
    for a in arrays {
        let (ci, local, _) = location(sequence, a)?;
        let width = a.aligned_consensus.len();
        write!(
            out,
            "\n\nArray {}\n>{}\n\n",
            a.id, sequence.contigs[ci].label
        )
        .unwrap();
        writeln!(
            out,
            "       Pos  Repeat     %id  Spacer  Left flank    Repeat{}    Spacer",
            " ".repeat(width.saturating_sub(6))
        )
        .unwrap();
        let rule = format!(
            "==========  ======  ======  ======  ==========    {}",
            "=".repeat(width)
        );
        writeln!(out, "{rule}    ======").unwrap();
        let mut pos = local;
        for (i, element) in a.elements.iter().enumerate() {
            let repeat = &element.repeat;
            let length = if mode == CoordinateMode::Legacy {
                repeat.len()
            } else {
                ungapped(repeat)
            };
            let same = repeat
                .iter()
                .zip(&a.aligned_consensus)
                .filter(|(a, b)| a.eq_ignore_ascii_case(b))
                .count();
            write!(
                out,
                "{:>10}  {:>6}  {:>6.1}  ",
                pos + 1,
                length,
                100.0 * same as f64 / width as f64
            )
            .unwrap();
            if i + 1 == a.elements.len() {
                out.push_str("        ");
            } else {
                write!(out, "{:>6}  ", element.spacer.len()).unwrap();
            }
            let flank = &element.left_flank;
            write!(
                out,
                "{:>10}    ",
                letters(&flank[flank.len().saturating_sub(10)..])
            )
            .unwrap();
            for (&c, &cons) in repeat.iter().zip(&a.aligned_consensus) {
                out.push(if c == cons && c != b'-' {
                    '.'
                } else {
                    c as char
                });
            }
            writeln!(out, "    {}", letters(&element.spacer)).unwrap();
            pos = pos
                .checked_add(length)
                .and_then(|p| p.checked_add(element.spacer.len()))
                .ok_or_else(|| "coordinate overflow".to_string())?;
        }
        let repeat_length = if mode == CoordinateMode::Legacy {
            width
        } else {
            a.elements
                .iter()
                .map(|e| ungapped(&e.repeat))
                .sum::<usize>()
                / a.elements.len()
        };
        writeln!(out, "{rule}").unwrap();
        writeln!(
            out,
            "{:>10}  {:>6}          {:>6}                {}",
            a.elements.len(),
            repeat_length,
            spacer_average(a),
            letters(&a.aligned_consensus)
        )
        .unwrap();
    }
    out.push_str("\n\nSUMMARY BY SIMILARITY\n\n\n\nArray          Sequence    Position      Length  # Copies  Repeat  Spacer  +  Consensus\n=====  ================  ==========  ==========  ========  ======  ======  =  =========\n");
    for cluster in clusters(arrays) {
        let rows: Vec<_> = cluster
            .iter()
            .map(|&(i, rev)| {
                if rev {
                    reverse_complement(&arrays[i].consensus)
                } else {
                    arrays[i].consensus.clone()
                }
            })
            .collect();
        let aligned = if rows.len() == 1 {
            rows
        } else {
            crate::msa::align(&rows)?.rows
        };
        for (&(i, reverse), row) in cluster.iter().zip(&aligned) {
            summary_prefix(&mut out, sequence, &arrays[i], mode);
            writeln!(out, "{}  {}", if reverse { '-' } else { '+' }, letters(row)).unwrap();
        }
        if cluster.len() > 1 {
            out.push_str(&" ".repeat(78));
            for col in 0..aligned[0].len() {
                let mut counts = [0usize; 5];
                for row in &aligned {
                    counts[match row[col].to_ascii_uppercase() {
                        b'A' => 0,
                        b'C' => 1,
                        b'G' => 2,
                        b'T' | b'U' => 3,
                        _ => 4,
                    }] += 1;
                }
                let mut best = 0;
                for i in 1..5 {
                    if counts[i] > counts[best] {
                        best = i;
                    }
                }
                let count = counts[best];
                out.push(if best == 4 {
                    ' '
                } else if count == aligned.len() {
                    '*'
                } else if count as f64 / aligned.len() as f64 >= 0.95
                    || (aligned.len() > 5 && count == aligned.len() - 1)
                {
                    ':'
                } else {
                    ' '
                });
            }
            out.push('\n');
        }
        out.push('\n');
    }
    out.push_str("\n\nSUMMARY BY POSITION\n\n");
    let mut sorted: Vec<_> = arrays.iter().collect();
    sorted.sort_by_key(|a| a.pos);
    let mut previous: Option<&Array> = None;
    let mut last_contig = None;
    for a in sorted {
        let (ci, _, _) = location(sequence, a)?;
        if last_contig != Some(ci) {
            write!(out, "\n\n>{}\n\nArray          Sequence    Position      Length  # Copies  Repeat  Spacer    Distance  Consensus\n=====  ================  ==========  ==========  ========  ======  ======  ==========  =========\n", sequence.contigs[ci].label).unwrap();
            last_contig = Some(ci);
        }
        summary_prefix(&mut out, sequence, a, mode);
        let distance = previous.and_then(|p| {
            let (pci, _, len) = location(sequence, p).ok()?;
            let end = p.pos.checked_add(len)?.checked_add(spacer_average(p))?;
            (pci == ci && sequence.global_to_local(end.checked_sub(1)?)?.0 == ci)
                .then_some(a.pos as i128 - end as i128)
        });
        match distance {
            Some(d) if d != -1 => write!(out, "{d:>10}  ").unwrap(),
            _ => out.push_str("            "),
        }
        writeln!(out, "{}", letters(&a.consensus)).unwrap();
        previous = Some(a);
    }
    Ok(out)
}

#[derive(serde::Serialize)]
struct JsonReport<'a> {
    arrays: Vec<JsonArray<'a>>,
}

#[derive(serde::Serialize)]
struct JsonArray<'a> {
    id: usize,
    sequence: &'a str,
    start: usize,
    end: usize,
    consensus: String,
    aligned_consensus: String,
    cluster_rev_comp: bool,
    repeats: Vec<JsonRepeat>,
    spacers: Vec<JsonSpacer>,
    right_flank: String,
}

#[derive(serde::Serialize)]
struct JsonRepeat {
    start: usize,
    end: usize,
    sequence: String,
    aligned_sequence: String,
}

#[derive(serde::Serialize)]
struct JsonSpacer {
    start: usize,
    end: usize,
    sequence: String,
}

pub fn json(sequence: &Sequence, arrays: &[Array]) -> Result<String, String> {
    let mut report = JsonReport {
        arrays: Vec::with_capacity(arrays.len()),
    };
    for a in arrays {
        let (ci, local, length) = location(sequence, a)?;
        let mut repeats = Vec::with_capacity(a.elements.len());
        let mut pos = local + 1;
        for element in &a.elements {
            let len = ungapped(&element.repeat);
            let sequence: Vec<_> = element
                .repeat
                .iter()
                .copied()
                .filter(|&c| c != b'-')
                .collect();
            repeats.push(JsonRepeat {
                start: pos,
                end: pos + len - 1,
                sequence: letters(&sequence),
                aligned_sequence: letters(&element.repeat),
            });
            pos += len + element.spacer.len();
        }
        let mut spacers = Vec::with_capacity(a.elements.len().saturating_sub(1));
        pos = local + 1;
        for i in 0..a.elements.len() - 1 {
            pos += ungapped(&a.elements[i].repeat);
            spacers.push(JsonSpacer {
                start: pos,
                end: pos + a.elements[i].spacer.len().saturating_sub(1),
                sequence: letters(&a.elements[i].spacer),
            });
            pos += a.elements[i].spacer.len();
        }
        report.arrays.push(JsonArray {
            id: a.id,
            sequence: &sequence.contigs[ci].label,
            start: local + 1,
            end: local + length,
            consensus: letters(&a.consensus),
            aligned_consensus: letters(&a.aligned_consensus),
            cluster_rev_comp: a.cluster_rev_comp,
            repeats,
            spacers,
            right_flank: letters(&a.elements.last().unwrap().spacer),
        });
    }
    let mut out = serde_json::to_string(&report).map_err(|e| e.to_string())?;
    out.push('\n');
    Ok(out)
}

fn percent_encoded(s: &str) -> String {
    let mut out = String::new();
    for b in s.bytes() {
        if b.is_ascii_alphanumeric()
            || matches!(
                b,
                b'.' | b':' | b'^' | b'*' | b'$' | b'@' | b'!' | b'+' | b'_' | b'?' | b'-' | b'|'
            )
        {
            out.push(b as char);
        } else {
            write!(out, "%{b:02X}").unwrap();
        }
    }
    out
}

pub fn gff(sequence: &Sequence, arrays: &[Array]) -> Result<String, String> {
    let mut out = String::from("##gff-version 3\n");
    for a in arrays {
        let (ci, local, length) = location(sequence, a)?;
        let label = percent_encoded(&sequence.contigs[ci].label);
        let parent = format!("array{}", a.id);
        writeln!(out, "{label}\tpilercr\trepeat_region\t{}\t{}\t.\t.\t.\tID={parent};Note=CRISPR;consensus={}", local + 1, local + length, percent_encoded(&letters(&a.consensus))).unwrap();
        let mut pos = local + 1;
        for (i, element) in a.elements.iter().enumerate() {
            let len = ungapped(&element.repeat);
            writeln!(out, "{label}\tpilercr\tdirect_repeat\t{pos}\t{}\t.\t.\t.\tID={parent}.repeat{};Parent={parent}", pos + len - 1, i + 1).unwrap();
            pos += len;
            let len = element.spacer.len();
            if i + 1 < a.elements.len() && len != 0 {
                writeln!(out, "{label}\tpilercr\tspacer\t{pos}\t{}\t.\t.\t.\tID={parent}.spacer{};Parent={parent}", pos + len - 1, i + 1).unwrap();
                pos += len;
            }
        }
    }
    Ok(out)
}

pub fn tsv(sequence: &Sequence, arrays: &[Array]) -> Result<String, String> {
    let mut out =
        String::from("array\tsequence\tstart\tend\tlength\tcopies\trepeat\tspacer\tconsensus\n");
    for a in arrays {
        let (ci, local, length) = location(sequence, a)?;
        let label = sequence.contigs[ci]
            .label
            .replace('\\', "\\\\")
            .replace('\t', "\\t")
            .replace('\r', "\\r")
            .replace('\n', "\\n");
        let repeat = a
            .elements
            .iter()
            .map(|e| ungapped(&e.repeat))
            .sum::<usize>()
            / a.elements.len();
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            a.id,
            label,
            local + 1,
            local + length,
            length,
            a.elements.len(),
            repeat,
            spacer_average(a),
            letters(&a.consensus)
        )
        .unwrap();
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::array::ArrayElement;
    use crate::seq::Contig;
    use std::{fs, path::Path};

    fn sample() -> (Sequence, Array) {
        let seq = Sequence {
            bases: vec![b'A'; 100],
            contigs: vec![Contig {
                start: 0,
                len: 100,
                label: "test".into(),
            }],
        };
        let a = Array {
            id: 1,
            pos: 2,
            consensus: b"ACGT".to_vec(),
            aligned_consensus: b"AC-GT".to_vec(),
            elements: vec![
                ArrayElement {
                    left_flank: b"AA".to_vec(),
                    repeat: b"AC-GT".to_vec(),
                    spacer: b"TT".to_vec(),
                },
                ArrayElement {
                    left_flank: b"AAAAAAAAAAAA".to_vec(),
                    repeat: b"Ac-GA".to_vec(),
                    spacer: b"GG".to_vec(),
                },
            ],
            cluster_rev_comp: false,
        };
        (seq, a)
    }

    #[test]
    fn gap_coordinates_and_mismatch_dots() {
        let (seq, a) = sample();
        let fixed = text(
            &seq,
            std::slice::from_ref(&a),
            "in",
            CoordinateMode::Corrected,
            true,
        )
        .unwrap();
        let legacy = text(
            &seq,
            std::slice::from_ref(&a),
            "in",
            CoordinateMode::Legacy,
            true,
        )
        .unwrap();
        assert!(fixed.contains("         9       4"));
        assert!(legacy.contains("        10       5"));
        assert!(fixed.contains("..-.."));
        assert!(fixed.contains(".c-.A"));
        assert!(fixed.contains("        AA    "));
        assert!(!fixed.contains("AAAAAAAAAAA"));
        let json = json(&seq, std::slice::from_ref(&a)).unwrap();
        assert!(json.contains("\"start\":9,\"end\":12"));
        assert!(json.contains("\"right_flank\":\"GG\""));
        let gff = gff(&seq, std::slice::from_ref(&a)).unwrap();
        assert!(gff.contains("direct_repeat\t9\t12"));
        assert_eq!(gff.lines().count(), 5);
        assert!(
            tsv(&seq, &[a])
                .unwrap()
                .contains("1\ttest\t3\t12\t10\t2\t4\t2\tACGT")
        );
    }

    #[test]
    fn no_array_reference_report() {
        let seq = Sequence {
            bases: vec![],
            contigs: vec![],
        };
        let expected =
            "pilercr v1.06\nBy Robert C. Edgar\n\nempty.fa: 0 putative CRISPR arrays found.\n\n";
        assert_eq!(
            text(&seq, &[], "empty.fa", CoordinateMode::Corrected, true).unwrap(),
            expected
        );
        assert_eq!(
            text(&seq, &[], "empty.fa", CoordinateMode::Legacy, false).unwrap(),
            format!("{INFO}{expected}")
        );
        assert_eq!(json(&seq, &[]).unwrap(), "{\"arrays\":[]}\n");
        assert_eq!(gff(&seq, &[]).unwrap(), "##gff-version 3\n");
    }

    #[test]
    fn safe_escaping_and_invalid_arrays() {
        assert_eq!(percent_encoded("a;b= c\té"), "a%3Bb%3D%20c%09%C3%A9");
        let (seq, mut a) = sample();
        a.pos = 99;
        assert!(text(&seq, &[a], "in", CoordinateMode::Corrected, true).is_err());
    }

    #[test]
    fn json_uses_serde_escaping_and_stable_field_order() {
        let (mut seq, a) = sample();
        let label = "chr\"\\\n\t\u{01}é";
        seq.contigs[0].label = label.into();
        let out = json(&seq, std::slice::from_ref(&a)).unwrap();
        let expected = "{\"arrays\":[{\"id\":1,\"sequence\":\"chr\\\"\\\\\\n\\t\\u0001é\",\"start\":3,\"end\":12,\"consensus\":\"ACGT\",\"aligned_consensus\":\"AC-GT\",\"cluster_rev_comp\":false,\"repeats\":[{\"start\":3,\"end\":6,\"sequence\":\"ACGT\",\"aligned_sequence\":\"AC-GT\"},{\"start\":9,\"end\":12,\"sequence\":\"AcGA\",\"aligned_sequence\":\"Ac-GA\"}],\"spacers\":[{\"start\":7,\"end\":8,\"sequence\":\"TT\"}],\"right_flank\":\"GG\"}]}\n";
        assert_eq!(out, expected);

        let parsed: serde_json::Value = serde_json::from_str(&out).unwrap();
        assert_eq!(parsed["arrays"].as_array().unwrap().len(), 1);
        assert_eq!(parsed["arrays"][0]["sequence"], label);
        assert_eq!(
            parsed["arrays"][0]["repeats"][1]["aligned_sequence"],
            "Ac-GA"
        );
        assert_eq!(parsed["arrays"][0]["spacers"][0]["sequence"], "TT");
        assert_eq!(parsed["arrays"][0]["right_flank"], "GG");
    }

    #[test]
    fn locate_uses_local_coordinates_and_summary_averages() {
        let sequence =
            Sequence::from_reader(&b">first\nAAAA\n>second\nACGTACGTACGT\n"[..]).unwrap();
        assert_eq!(sequence.contigs[0].start, 0);
        assert_eq!(sequence.contigs[1].start, 64);
        let a = Array {
            id: 1,
            pos: sequence.contigs[1].start + 2,
            consensus: b"ACGT".to_vec(),
            aligned_consensus: b"ACGT".to_vec(),
            elements: vec![
                ArrayElement {
                    left_flank: b"AAA".to_vec(),
                    repeat: b"ACGT".to_vec(),
                    spacer: b"TT".to_vec(),
                },
                ArrayElement {
                    left_flank: b"CCC".to_vec(),
                    repeat: b"AC-T".to_vec(),
                    spacer: b"GG".to_vec(),
                },
            ],
            cluster_rev_comp: false,
        };
        let expected = ArrayLocation {
            contig: 1,
            start: 2,
            len: 9,
            copies: 2,
            repeat_len: 3,
            spacer_len: 2,
        };
        assert_eq!(locate(&sequence, &a), Ok(expected));
        assert_eq!(crate::locate(&sequence, &a), Ok(expected));
        assert_eq!(location(&sequence, &a), Ok((1, 2, 9)));
    }

    #[test]
    fn locate_rejects_invalid_arrays() {
        let (sequence, mut a) = sample();
        a.pos = 99;
        assert!(locate(&sequence, &a).is_err());
        let (sequence, mut a) = sample();
        a.elements.clear();
        assert!(locate(&sequence, &a).is_err());
    }

    #[test]
    fn canonical_golden_reports() {
        let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/excerpts");
        let mut ids: Vec<_> = fs::read_dir(&dir)
            .unwrap()
            .filter_map(|entry| {
                let name = entry.unwrap().file_name().into_string().unwrap();
                name.strip_prefix("excerpt_")?
                    .strip_suffix(".fa")
                    .filter(|id| !id.contains('.'))
                    .map(str::to_owned)
            })
            .collect();
        ids.sort();
        assert!(!ids.is_empty(), "{}", dir.display());
        for id in ids {
            let reference =
                fs::read_to_string(dir.join(format!("excerpt_{id}.legacy.txt"))).unwrap();
            let lines: Vec<_> = reference.lines().collect();
            let mut arrays = Vec::new();
            let mut i = lines
                .iter()
                .position(|&l| l == "DETAIL REPORT")
                .unwrap_or(lines.len())
                + 1;
            while i < lines.len() && lines[i] != "SUMMARY BY SIMILARITY" {
                if let Some(id) = lines[i].strip_prefix("Array ") {
                    let id = id.parse().unwrap();
                    i += 5;
                    let mut repeats = Vec::new();
                    let mut spacers = Vec::new();
                    let mut flanks = Vec::new();
                    let mut pos = 0;
                    while !lines[i].starts_with("==========") {
                        let l = lines[i];
                        if repeats.is_empty() {
                            pos = l[..10].trim().parse::<usize>().unwrap() - 1;
                        }
                        let fields: Vec<_> = l[34..].split_whitespace().collect();
                        flanks.push(fields[0].as_bytes().to_vec());
                        repeats.push(fields[1].as_bytes().to_vec());
                        spacers.push(fields.get(2).unwrap_or(&"").as_bytes().to_vec());
                        i += 1;
                    }
                    let consensus = lines[i + 1]
                        .split_whitespace()
                        .last()
                        .unwrap()
                        .as_bytes()
                        .to_vec();
                    for repeat in &mut repeats {
                        for (c, &cons) in repeat.iter_mut().zip(&consensus) {
                            if *c == b'.' {
                                *c = cons;
                            }
                        }
                    }
                    let elements = flanks
                        .into_iter()
                        .zip(repeats)
                        .zip(spacers)
                        .map(|((left_flank, repeat), spacer)| ArrayElement {
                            left_flank,
                            repeat,
                            spacer,
                        })
                        .collect();
                    arrays.push(Array {
                        id,
                        pos,
                        consensus: consensus.clone(),
                        aligned_consensus: consensus,
                        elements,
                        cluster_rev_comp: false,
                    });
                }
                i += 1;
            }
            let seq = Sequence::from_path(dir.join(format!("excerpt_{id}.fa"))).unwrap();
            let input = lines
                .iter()
                .find_map(|l| {
                    l.split_once(": ")
                        .filter(|(_, r)| r.ends_with("putative CRISPR arrays found."))
                        .map(|(l, _)| l)
                })
                .unwrap();
            let actual = text(&seq, &arrays, input, CoordinateMode::Legacy, false).unwrap();
            if actual != reference {
                let mismatch = actual
                    .lines()
                    .zip(reference.lines())
                    .position(|(a, b)| a != b)
                    .unwrap_or(actual.lines().count().min(reference.lines().count()));
                panic!(
                    "golden mismatch at line {}: actual {:?}, expected {:?}",
                    mismatch + 1,
                    actual.lines().nth(mismatch),
                    reference.lines().nth(mismatch)
                );
            }
        }
    }
}

use pilercr_rs::{
    Array, Config, Sequence, detect,
    report::{self, CoordinateMode},
};
use serde_json::Value;
use std::{
    fs,
    path::{Path, PathBuf},
};

fn verify_report(sequence: &Sequence, arrays: &[Array], input: &str) -> usize {
    assert_eq!(CoordinateMode::default(), CoordinateMode::Corrected);
    let legacy = report::text(sequence, arrays, input, CoordinateMode::Legacy, true).unwrap();
    let structured: Value = serde_json::from_str(&report::json(sequence, arrays).unwrap()).unwrap();
    let json_arrays = structured["arrays"].as_array().unwrap();
    assert_eq!(json_arrays.len(), arrays.len(), "{input}");
    if arrays.is_empty() {
        assert!(pilercr_parser::parse(&legacy).is_err(), "{input}");
        assert!(legacy.contains("0 putative CRISPR arrays found."));
        return 0;
    }
    let parsed =
        pilercr_parser::parse(&legacy).unwrap_or_else(|error| panic!("{input}: {error:?}"));
    assert_eq!(parsed.len(), arrays.len(), "{input}");
    let corrected = report::text(sequence, arrays, input, CoordinateMode::default(), true).unwrap();
    let mut corrected_rows = corrected.lines().filter_map(|line| {
        let fields: Vec<_> = line.split_whitespace().collect();
        (fields.len() >= 5 && fields[2].contains('.') && fields[2].parse::<f64>().is_ok())
            .then_some(fields)
    });
    let mut gaps = 0;
    for ((native, parsed), json) in arrays.iter().zip(&parsed).zip(json_arrays) {
        let (contig_index, local) = sequence.global_to_local(native.pos).unwrap();
        let contig = &sequence.contigs[contig_index];
        assert_eq!(parsed.accession, contig.label, "{input}");
        assert_eq!(parsed.order + 1, native.id, "{input}");
        assert_eq!(parsed.start, local, "{input}");
        assert_eq!(
            parsed.consensus_repeat_sequence.as_bytes(),
            native.aligned_consensus
        );
        assert_eq!(parsed.repeat_spacers.len(), native.elements.len());
        assert_eq!(json["id"], native.id);
        assert_eq!(json["sequence"], contig.label);
        assert_eq!(json["start"], local + 1);
        assert_eq!(
            json["repeats"].as_array().unwrap().len(),
            native.elements.len()
        );
        assert_eq!(
            json["spacers"].as_array().unwrap().len(),
            native.elements.len() - 1
        );
        let mut pos = local;
        let mut array_end = local;
        for (index, rs) in parsed.repeat_spacers.iter().enumerate() {
            let element = &native.elements[index];
            let repeat: Vec<_> = element
                .repeat
                .iter()
                .copied()
                .filter(|&base| base != b'-')
                .collect();
            gaps += element.repeat.len() - repeat.len();
            let repeat_end = pos + repeat.len();
            let spacer_end = repeat_end + element.spacer.len();
            assert_eq!((rs.start, rs.repeat_start), (pos, pos), "{input}");
            assert_eq!(rs.repeat_end, repeat_end, "{input}");
            assert_eq!(
                (rs.spacer_start, rs.spacer_end, rs.end),
                (repeat_end, spacer_end, spacer_end)
            );
            assert_eq!(rs.repeat.as_bytes(), repeat, "{input}");
            assert_eq!(rs.spacer.as_bytes(), element.spacer, "{input}");
            assert!(
                sequence.bases[contig.start + pos..contig.start + repeat_end]
                    .eq_ignore_ascii_case(rs.repeat.as_bytes()),
                "{input}: repeat {index}"
            );
            assert!(
                sequence.bases[contig.start + repeat_end..contig.start + spacer_end]
                    .eq_ignore_ascii_case(rs.spacer.as_bytes()),
                "{input}: spacer {index}"
            );
            let row = corrected_rows.next().unwrap();
            assert_eq!(row[0].parse::<usize>().unwrap(), pos + 1, "{input}");
            assert_eq!(row[1].parse::<usize>().unwrap(), repeat.len(), "{input}");
            let json_repeat = &json["repeats"][index];
            assert_eq!(json_repeat["start"], rs.repeat_start + 1);
            assert_eq!(json_repeat["end"], rs.repeat_end);
            assert_eq!(json_repeat["sequence"], rs.repeat);
            if index + 1 < native.elements.len() {
                let spacer = &json["spacers"][index];
                assert_eq!(spacer["start"], rs.spacer_start + 1);
                assert_eq!(spacer["end"], rs.spacer_end);
                assert_eq!(spacer["sequence"], rs.spacer);
            } else {
                assert_eq!(json["right_flank"], rs.spacer);
            }
            array_end = repeat_end;
            pos = spacer_end;
        }
        assert_eq!(parsed.end, pos, "{input}");
        assert_eq!(json["end"], array_end);
    }
    assert!(corrected_rows.next().is_none());
    gaps
}

/// Verify the native reports of `inputs` (label, path) against the upstream
/// parser. Returns (positive, empty, corrected gaps, failures).
fn verify_inputs(inputs: &[(String, PathBuf)]) -> (usize, usize, usize, Vec<String>) {
    let mut failures = Vec::new();
    let mut gaps = 0;
    let mut positive = 0;
    let mut empty = 0;
    for (input, path) in inputs {
        let result = std::panic::catch_unwind(|| {
            let sequence = Sequence::from_path(path).unwrap();
            let arrays = detect(&sequence, &Config::default()).unwrap();
            (verify_report(&sequence, &arrays, input), arrays.is_empty())
        });
        match result {
            Ok((count, is_empty)) => {
                gaps += count;
                if is_empty {
                    empty += 1;
                } else {
                    positive += 1;
                }
            }
            Err(error) => failures.push(format!(
                "{input}: {}",
                error
                    .downcast_ref::<String>()
                    .map(String::as_str)
                    .or_else(|| error.downcast_ref::<&str>().copied())
                    .unwrap_or("non-string panic")
            )),
        }
    }
    (positive, empty, gaps, failures)
}

#[test]
fn committed_excerpt_legacy_reports_recover_genomic_coordinates() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/excerpts");
    let inputs: Vec<_> = fs::read_dir(&dir)
        .unwrap()
        .filter_map(|entry| {
            let name = entry.ok()?.file_name().into_string().ok()?;
            let id = name.strip_prefix("excerpt_")?.strip_suffix(".fa")?;
            if id.contains('.') {
                return None;
            }
            Some((format!("excerpt_{id}.fa"), dir.join(&name)))
        })
        .collect();
    assert!(!inputs.is_empty(), "no excerpts in {}", dir.display());
    let (positive, empty, gaps, failures) = verify_inputs(&inputs);
    eprintln!(
        "committed excerpt parser integration: {positive} positive, {empty} empty reports, {gaps} corrected gaps"
    );
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    assert_eq!(positive + empty, inputs.len());
    assert!(positive > 0);
}

#[test]
#[ignore = "requires test-data/; run with --ignored"]
fn all_52_native_legacy_reports_recover_genomic_coordinates() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let inputs: Vec<_> = (1..=52)
        .map(|id| {
            (
                format!("{id}.fasta"),
                root.join("test-data/assemblies")
                    .join(format!("{id}.fasta")),
            )
        })
        .collect();
    let (positive, empty, gaps, failures) = verify_inputs(&inputs);
    eprintln!(
        "parser integration: {positive} positive, {empty} empty reports, {gaps} corrected gaps"
    );
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    assert_eq!(positive + empty, 52);
    assert!(positive > 0);
    assert!(gaps > 0);
}

#[test]
fn zero_array_report_is_valid_but_unsupported_by_upstream_parser() {
    let sequence = Sequence::from_reader(&b">empty\nACGT\n"[..]).unwrap();
    assert_eq!(verify_report(&sequence, &[], "empty.fasta"), 0);
}

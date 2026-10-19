//! Independent GFF3 parsing and structural validation for `report::gff`.
//!
//! The default test uses `noodles-gff` as an independent GFF3 parser. It checks
//! directives, percent-encoded sequence IDs and attribute values, 1-based
//! inclusive coordinates, and parent/child referential integrity.
//!
//! A second, ignored test invokes the canonical GenomeTools validator
//! (`gt gff3validator`). It is opt-in because `gt` is not part of this
//! repository's build environment.

use std::collections::HashMap;
use std::fs;
use std::process::Command;

use noodles_gff::feature::record::Strand;
use noodles_gff::{self as gff, LineBuf};
use pilercr_rs::report;
use pilercr_rs::seq::Contig;
use pilercr_rs::{Array, Sequence};

const TRICKY_LABEL: &str = "chr1;alpha\tbeta%\u{03b1}";
const TRICKY_CONSENSUS: &[u8] = b"AC;GT=%";

fn sequence() -> Sequence {
    let mut bases = vec![b'A'; 400];
    // Keep the generated bases deterministic and sequence-like.
    for (index, base) in bases.iter_mut().enumerate() {
        *base = b"ACGT"[index % 4];
    }
    Sequence {
        bases,
        contigs: vec![
            Contig {
                start: 0,
                len: 150,
                label: TRICKY_LABEL.into(),
            },
            Contig {
                start: 150,
                len: 250,
                label: "chr2".into(),
            },
        ],
    }
}

fn array(id: usize, pos: usize, consensus: &[u8], spacer_len: usize) -> Array {
    let repeats = vec![b"ACGT".to_vec(); 3];
    let spacers = vec![vec![b'T'; spacer_len]; 3];
    Array {
        id,
        pos,
        consensus: consensus.to_vec(),
        aligned_consensus: b"ACGT".to_vec(),
        left_flanks: vec![Vec::new(); 3],
        repeats,
        spacers,
        cluster_rev_comp: false,
    }
}

fn representative_output() -> String {
    let sequence = sequence();
    let arrays = [
        // Contig 0: local 10, so the 1-based region is 11..=28.
        array(3, 10, TRICKY_CONSENSUS, 3),
        // Contig 1: global 180, local 30, so the 1-based region is 31..=46.
        array(7, 180, b"TTAA", 2),
    ];
    report::gff(&sequence, &arrays).unwrap()
}

#[derive(Debug)]
struct Feature {
    seqid: Vec<u8>,
    ty: String,
    start: usize,
    end: usize,
    id: String,
    parent: Option<String>,
}

fn as_bytes<T: AsRef<[u8]> + ?Sized>(value: &T) -> &[u8] {
    value.as_ref()
}

fn attribute<'a>(record: &'a gff::feature::RecordBuf, tag: &[u8]) -> Option<&'a [u8]> {
    let value = record.attributes().get(tag)?;
    let value = value.as_string()?;
    Some(as_bytes(value))
}

fn parse(gff_text: &str) -> (Vec<gff::DirectiveBuf>, Vec<Feature>) {
    let mut reader = gff::io::Reader::new(gff_text.as_bytes());
    let mut directives = Vec::new();
    let mut features = Vec::new();
    let mut saw_record = false;

    for line in reader.line_bufs() {
        match line.unwrap() {
            LineBuf::Directive(directive) => {
                assert!(
                    !saw_record,
                    "directive {:?} appears after the first record",
                    directive.key()
                );
                directives.push(directive);
            }
            LineBuf::Comment(comment) => panic!("unexpected comment: {comment:?}"),
            LineBuf::Record(record) => {
                saw_record = true;
                features.push(Feature {
                    seqid: record.reference_sequence_name().to_vec(),
                    ty: String::from_utf8(record.ty().to_vec()).unwrap(),
                    start: usize::from(record.start()),
                    end: usize::from(record.end()),
                    id: String::from_utf8(attribute(&record, b"ID").unwrap().to_vec()).unwrap(),
                    parent: attribute(&record, b"Parent")
                        .map(|value| String::from_utf8(value.to_vec()).unwrap()),
                });
            }
        }
    }

    (directives, features)
}

#[test]
fn parses_as_valid_gff3_with_consistent_parent_child_coordinates() {
    let output = representative_output();

    let (directives, features) = parse(&output);

    assert_eq!(directives.len(), 1);
    let version = &directives[0];
    assert_eq!(as_bytes(version.key()), b"gff-version");
    match version.value().expect("missing gff-version value") {
        gff::directive_buf::value::Value::String(value) => {
            assert_eq!(as_bytes(value), b"3");
        }
        other => panic!("unexpected gff-version value: {other:?}"),
    }
    assert_eq!(output.lines().next(), Some("##gff-version 3"));
    assert_eq!(
        output.lines().filter(|line| line.starts_with("##")).count(),
        1
    );

    // Two arrays, each with one region, three repeats, and two spacers.
    assert_eq!(features.len(), 12);

    let mut by_id: HashMap<&str, &Feature> = HashMap::new();
    for feature in &features {
        assert!(feature.start >= 1, "coordinates are 1-based: {feature:?}");
        assert!(
            feature.start <= feature.end,
            "inclusive coordinates: {feature:?}"
        );
        assert!(
            by_id.insert(&feature.id, feature).is_none(),
            "duplicate ID {}",
            feature.id
        );
    }

    for feature in &features {
        if let Some(parent_id) = &feature.parent {
            let parent = by_id
                .get(parent_id.as_str())
                .unwrap_or_else(|| panic!("missing parent {parent_id}"));
            assert_eq!(parent.ty, "repeat_region", "parent type for {feature:?}");
            assert_eq!(
                feature.seqid, parent.seqid,
                "cross-contig parent: {feature:?}"
            );
            assert!(
                parent.start <= feature.start && feature.end <= parent.end,
                "child outside parent: {feature:?} in {parent:?}"
            );
        } else {
            assert_eq!(feature.ty, "repeat_region", "orphan feature: {feature:?}");
        }
    }

    let mut children: HashMap<&str, Vec<&Feature>> = HashMap::new();
    for feature in &features {
        if let Some(parent_id) = &feature.parent {
            children
                .entry(parent_id.as_str())
                .or_default()
                .push(feature);
        }
    }

    for region in features.iter().filter(|f| f.ty == "repeat_region") {
        let mut kids = children
            .get(region.id.as_str())
            .cloned()
            .unwrap_or_else(|| panic!("region {} has no children", region.id));
        kids.sort_by_key(|feature| feature.start);

        assert_eq!(
            kids.iter().map(|f| f.ty.as_str()).collect::<Vec<_>>(),
            [
                "direct_repeat",
                "spacer",
                "direct_repeat",
                "spacer",
                "direct_repeat"
            ]
        );
        assert_eq!(kids.first().unwrap().start, region.start);
        assert_eq!(kids.last().unwrap().end, region.end);
        for pair in kids.windows(2) {
            assert_eq!(
                pair[0].end + 1,
                pair[1].start,
                "children must tile the region without gaps or overlaps: {pair:?}"
            );
        }
    }

    let chr1 = by_id.get("array3").expect("missing array3");
    assert_eq!(chr1.seqid, TRICKY_LABEL.as_bytes());
    assert_eq!((chr1.start, chr1.end), (11, 28));
    let chr2 = by_id.get("array7").expect("missing array7");
    assert_eq!(chr2.seqid, b"chr2");
    assert_eq!((chr2.start, chr2.end), (31, 46));

    // The raw output must percent-encode everything GFF3 reserves in a seqid.
    assert!(
        output.contains("chr1%3Balpha%09beta%25%CE%B1\tpilercr\trepeat_region\t"),
        "seqid was not percent-encoded:\n{output}"
    );
    // ... and in an attribute value (`;`, `=`, and `%` are reserved).
    assert!(
        output.contains("consensus=AC%3BGT%3D%25"),
        "consensus was not percent-encoded:\n{output}"
    );
}

#[test]
fn parses_records_after_percent_decoding() {
    let output = representative_output();
    let (_, features) = parse(&output);

    // `noodles-gff` exposes decoded values, so these assertions prove the
    // exporter both escaped the raw bytes and preserved their meaning.
    let region = features
        .iter()
        .find(|feature| feature.id == "array3")
        .expect("missing array3");
    assert_eq!(region.seqid, TRICKY_LABEL.as_bytes());

    let mut reader = gff::io::Reader::new(output.as_bytes());
    let mut saw_consensus = false;
    for line in reader.line_bufs() {
        if let LineBuf::Record(record) = line.unwrap()
            && attribute(&record, b"ID") == Some(b"array3")
        {
            assert!(attribute(&record, b"Note") == Some(b"CRISPR"));
            assert_eq!(attribute(&record, b"consensus"), Some(TRICKY_CONSENSUS));
            assert_eq!(record.score(), None);
            assert_eq!(record.strand(), Strand::None);
            assert_eq!(record.phase(), None);
            assert_eq!(as_bytes(record.source()), b"pilercr");
            assert_eq!(as_bytes(record.ty()), b"repeat_region");
            saw_consensus = true;
        }
    }
    assert!(saw_consensus, "array3 region was not parsed");
}

#[test]
fn empty_output_is_a_valid_stream_with_only_the_version_directive() {
    let output = report::gff(
        &Sequence {
            bases: Vec::new(),
            contigs: Vec::new(),
        },
        &[],
    )
    .unwrap();
    let (directives, features) = parse(&output);
    assert_eq!(directives.len(), 1);
    assert_eq!(as_bytes(directives[0].key()), b"gff-version");
    assert!(features.is_empty());
}

/// Canonical validation with GenomeTools. This is intentionally ignored: `gt`
/// is not installed in the standard development environment. Install
/// GenomeTools and run:
///
/// ```text
/// cargo test -p pilercr-rs --test gff3_validation -- --ignored
/// ```
#[test]
#[ignore = "requires GenomeTools (gt) on PATH"]
fn passes_genometools_gff3validator() {
    let output = representative_output();
    let path = std::env::temp_dir().join(format!(
        "pilercr-rs-gff3-validation-{}.gff3",
        std::process::id()
    ));
    fs::write(&path, &output).unwrap();

    let result = Command::new("gt")
        .args(["gff3validator", path.to_str().unwrap()])
        .output();
    let result = result.unwrap_or_else(|error| {
        panic!(
            "failed to run 'gt gff3validator' ({error}); install GenomeTools or run the default parser test instead"
        )
    });
    let _ = fs::remove_file(&path);

    assert!(
        result.status.success(),
        "gt gff3validator failed with {}\nstdout:\n{}\nstderr:\n{}",
        result.status,
        String::from_utf8_lossy(&result.stdout),
        String::from_utf8_lossy(&result.stderr)
    );
}

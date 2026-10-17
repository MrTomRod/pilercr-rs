use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sequence {
    pub bases: Vec<u8>,
    pub contigs: Vec<Contig>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Contig {
    pub label: String,
    pub start: usize,
    pub len: usize,
}

impl Sequence {
    pub fn from_reader(mut reader: impl BufRead) -> io::Result<Self> {
        let mut sequence = Self {
            bases: Vec::new(),
            contigs: Vec::new(),
        };
        let mut line = Vec::new();
        loop {
            line.clear();
            if reader.read_until(b'\n', &mut line)? == 0 {
                break;
            }
            if line[0] == b'>' {
                if line.last() != Some(&b'\n') {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "unfinished FASTA header",
                    ));
                }
                let label: Vec<u8> = line[1..line.len() - 1]
                    .iter()
                    .copied()
                    .filter(|&byte| byte != b'\r')
                    .collect();
                let label = String::from_utf8(label)
                    .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
                if !sequence.contigs.is_empty() {
                    let mut padding = 32 - sequence.bases.len() % 32;
                    if padding < 50 {
                        padding += 32;
                    }
                    sequence.bases.resize(sequence.bases.len() + padding, b'!');
                }
                sequence.contigs.push(Contig {
                    label,
                    start: sequence.bases.len(),
                    len: 0,
                });
            } else {
                for &byte in &line {
                    if matches!(byte, b' ' | b'\t'..=b'\r') {
                        continue;
                    }
                    let contig = sequence.contigs.last_mut().ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            "sequence data before FASTA header",
                        )
                    })?;
                    sequence.bases.push(byte);
                    contig.len += 1;
                }
            }
        }
        Ok(sequence)
    }

    pub fn from_path(path: impl AsRef<Path>) -> io::Result<Self> {
        Self::from_reader(BufReader::new(File::open(path)?))
    }

    pub fn global_to_local(&self, pos: usize) -> Option<(usize, usize)> {
        if pos >= self.bases.len() {
            return None;
        }
        let index = self.contigs.partition_point(|contig| contig.start <= pos);
        let index = index.checked_sub(1)?;
        let contig = &self.contigs[index];
        let local = pos - contig.start;
        (local < contig.len).then_some((index, local))
    }

    pub fn local_to_global(&self, contig: usize, local: usize) -> Option<usize> {
        let contig = self.contigs.get(contig)?;
        if local >= contig.len {
            return None;
        }
        contig
            .start
            .checked_add(local)
            .filter(|&pos| pos < self.bases.len())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{fs, io::Cursor};

    fn parse(input: &[u8]) -> Sequence {
        Sequence::from_reader(Cursor::new(input)).unwrap()
    }

    #[test]
    fn padding_matches_reference_across_all_residues() {
        for len in 0..=96 {
            let mut input = b">first\n".to_vec();
            input.extend(vec![b'a'; len]);
            input.extend_from_slice(b"\n>second\nCg\n>third\nT");
            let sequence = parse(&input);
            let second_start = len + 64 - len % 32;
            assert_eq!(sequence.contigs[0].len, len);
            assert_eq!(sequence.contigs[1].start, second_start);
            assert_eq!(sequence.contigs[2].start, second_start + 64);
            assert!(sequence.bases[len..second_start].iter().all(|&b| b == b'!'));
            assert!(
                sequence.bases[second_start + 2..second_start + 64]
                    .iter()
                    .all(|&b| b == b'!')
            );
            assert_eq!(sequence.bases.len(), second_start + 65);
        }
    }

    #[test]
    fn coordinates_round_trip_and_reject_every_padding_position() {
        let sequence = parse(b">one\nAcG\n>two\nTn\n>three\na");
        for (index, contig) in sequence.contigs.iter().enumerate() {
            for local in 0..contig.len {
                let global = sequence.local_to_global(index, local).unwrap();
                assert_eq!(global, contig.start + local);
                assert_eq!(sequence.global_to_local(global), Some((index, local)));
            }
            assert_eq!(sequence.local_to_global(index, contig.len), None);
            assert_eq!(sequence.local_to_global(index, usize::MAX), None);
        }
        for global in 0..sequence.bases.len() {
            if sequence.bases[global] == b'!' {
                assert_eq!(sequence.global_to_local(global), None);
            }
        }
        assert_eq!(sequence.global_to_local(sequence.bases.len()), None);
        assert_eq!(sequence.global_to_local(usize::MAX), None);
        assert_eq!(sequence.local_to_global(sequence.contigs.len(), 0), None);
        assert_eq!(sequence.local_to_global(usize::MAX, 0), None);
    }

    #[test]
    fn preserves_full_utf8_headers_case_and_crlf() {
        let sequence =
            parse(">  first café\tannotation  \r\nAc g\tT\r\n>second\r label\r\nnN\r\n".as_bytes());
        assert_eq!(sequence.contigs[0].label, "  first café\tannotation  ");
        assert_eq!(sequence.contigs[1].label, "second label");
        assert_eq!(&sequence.bases[..4], b"AcgT");
        assert_eq!(&sequence.bases[64..], b"nN");
    }

    #[test]
    fn preserves_ambiguity_punctuation_and_non_utf8_sequence_bytes() {
        let sequence = parse(b">ambiguity\nACGTNnRrYySWKMBDHVUu-.!?\xff\nA>B\x0bC\x0cD");
        assert_eq!(sequence.bases, b"ACGTNnRrYySWKMBDHVUu-.!?\xffA>BCD");
        assert_eq!(sequence.contigs[0].len, sequence.bases.len());
        let bang = sequence.bases.iter().position(|&b| b == b'!').unwrap();
        assert_eq!(sequence.global_to_local(bang), Some((0, bang)));
    }

    #[test]
    fn empty_and_whitespace_only_input_succeed() {
        for input in [b"".as_slice(), b" \t\r\n\x0b\x0c\n"] {
            let sequence = parse(input);
            assert!(sequence.bases.is_empty());
            assert!(sequence.contigs.is_empty());
            assert_eq!(sequence.global_to_local(0), None);
            assert_eq!(sequence.local_to_global(0, 0), None);
        }
    }

    #[test]
    fn empty_contigs_have_no_valid_coordinates() {
        let sequence = parse(b">\n>empty\n>base\nA\n>last\n");
        assert_eq!(sequence.contigs.len(), 4);
        assert_eq!(sequence.contigs[0].label, "");
        assert_eq!(sequence.contigs[0].start, 0);
        assert_eq!(sequence.contigs[1].start, 64);
        assert_eq!(sequence.contigs[2].start, 128);
        assert_eq!(sequence.contigs[3].start, 192);
        for index in [0, 1, 3] {
            let contig = &sequence.contigs[index];
            assert_eq!(contig.len, 0);
            assert_eq!(sequence.global_to_local(contig.start), None);
            assert_eq!(sequence.local_to_global(index, 0), None);
        }
        assert_eq!(sequence.global_to_local(128), Some((2, 0)));
        let single = parse(b">only\n");
        assert!(single.bases.is_empty());
        assert_eq!(single.contigs[0].len, 0);
        assert_eq!(single.global_to_local(0), None);
    }

    #[test]
    fn malformed_input_is_invalid_data() {
        for input in [
            b"ACGT\n".as_slice(),
            b"\n A\n>header\n",
            b" >indented\nA",
            b">unfinished",
            b">unfinished\r",
            b">valid\nA\n>unfinished",
            b">\xff\nA",
        ] {
            let error = Sequence::from_reader(Cursor::new(input)).unwrap_err();
            assert_eq!(error.kind(), io::ErrorKind::InvalidData, "{input:?}");
        }
    }

    #[test]
    fn reader_errors_are_propagated() {
        struct FailingReader;
        impl io::Read for FailingReader {
            fn read(&mut self, _: &mut [u8]) -> io::Result<usize> {
                Err(io::Error::other("read failed"))
            }
        }
        let reader = BufReader::new(FailingReader);
        let error = Sequence::from_reader(reader).unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::Other);
    }

    #[test]
    fn loads_committed_excerpts_from_paths() {
        let dir = match option_env!("CARGO_MANIFEST_DIR") {
            Some(manifest) => Path::new(manifest).join("tests/fixtures/excerpts"),
            None => Path::new("crates/pilercr-core/tests/fixtures/excerpts").to_path_buf(),
        };
        let mut paths: Vec<_> = fs::read_dir(&dir)
            .unwrap()
            .map(|entry| entry.unwrap().path())
            .filter(|path| {
                path.extension().is_some_and(|ext| ext == "fa")
                    && path
                        .file_stem()
                        .is_some_and(|stem| !stem.to_string_lossy().contains('.'))
            })
            .collect();
        paths.sort();
        assert!(!paths.is_empty(), "{}", dir.display());
        for path in paths {
            let sequence = Sequence::from_path(&path).unwrap();
            assert!(!sequence.contigs.is_empty(), "{}", path.display());
            for (index, contig) in sequence.contigs.iter().enumerate() {
                assert_eq!(contig.start % 32, 0);
                if contig.len != 0 {
                    assert_eq!(sequence.global_to_local(contig.start), Some((index, 0)));
                    assert_eq!(
                        sequence.global_to_local(contig.start + contig.len - 1),
                        Some((index, contig.len - 1))
                    );
                }
            }
        }
    }
}

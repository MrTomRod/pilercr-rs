pub mod array;
pub mod bounds;
pub mod dp;
pub mod fasta_out;
pub mod filter;
pub mod graph;
pub mod kmer;
pub mod msa;
pub mod params;
pub mod pile;
pub mod report;
pub mod seq;
mod sort;
pub mod trace;

pub use array::{Array, ArrayElement, ArrayParams};
pub use dp::DpParams;
pub use report::{ArrayLocation, locate};
pub use seq::Sequence;
pub use sort::SortMode;

#[derive(Clone, Debug, PartialEq)]
pub struct Config {
    pub array: ArrayParams,
    pub dp: DpParams,
    pub diameter: usize,
    pub memory_budget: usize,
    /// Tie-breaking behaviour for the unstable sorts in detection.
    ///
    /// Defaults to [`SortMode::Modern`]; use [`SortMode::Legacy`] to target
    /// the shipped 2011 executable's GCC 4.4.3 sort order.
    pub sort_mode: SortMode,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            array: ArrayParams::default(),
            dp: DpParams::default(),
            diameter: 512,
            memory_budget: params::DEFAULT_MEMORY_BUDGET,
            sort_mode: SortMode::default(),
        }
    }
}

impl Config {
    pub fn validate(&self) -> Result<(), String> {
        let array = &self.array;
        if array.min_repeat == 0
            || array.min_repeat > array.max_repeat
            || array.min_spacer == 0
            || array.min_spacer > array.max_spacer
            || array.min_array < 3
        {
            return Err(
                "invalid array length bounds or minimum array size (must be at least 3)".into(),
            );
        }
        for ratio in [
            array.min_repeat_ratio,
            array.min_spacer_ratio,
            array.min_cons,
        ] {
            if !ratio.is_finite() || !(0.0..=1.0).contains(&ratio) {
                return Err(
                    "array ratios and conservation must be finite and between 0 and 1".into(),
                );
            }
        }
        params::get_params(
            0,
            self.diameter,
            self.dp.min_hit_length,
            self.dp.min_identity,
            self.memory_budget,
        )?;
        dp::align_traps(
            &Sequence {
                bases: Vec::new(),
                contigs: Vec::new(),
            },
            &[],
            &self.dp,
        )?;
        Ok(())
    }
}

pub fn detect(sequence: &Sequence, config: &Config) -> Result<Vec<Array>, String> {
    config.validate()?;
    dp::align_traps(sequence, &[], &config.dp)?;
    let diameter = config.diameter.min(sequence.bases.len()).max(1);
    let (filter_params, _) = params::get_params(
        sequence.bases.len(),
        diameter,
        config.dp.min_hit_length,
        config.dp.min_identity,
        config.memory_budget,
    )?;
    if trace::enabled() {
        trace::emit(
            "01-parameters",
            &[
                filter_params.word_size.to_string(),
                filter_params.seed_length.to_string(),
                filter_params.seed_diffs.to_string(),
                filter_params.tube_offset.to_string(),
                filter_params.diameter.to_string(),
                config.dp.min_hit_length.to_string(),
                format!("{:.17}", config.dp.min_identity),
                config.dp.max_hit_length.to_string(),
            ],
        );
        for (index, contig) in sequence.contigs.iter().enumerate() {
            trace::emit(
                "02-contigs",
                &[
                    index.to_string(),
                    contig.start.to_string(),
                    contig.len.to_string(),
                    contig.label.clone(),
                ],
            );
        }
    }
    if diameter < filter_params.word_size || sequence.bases.len() <= diameter {
        return Ok(Vec::new());
    }
    let filter_hits = kmer::filter(&sequence.bases, &filter_params)?;
    if trace::enabled() {
        for (index, hit) in filter_hits.iter().enumerate() {
            trace::emit(
                "03-filter-hits",
                &[
                    index.to_string(),
                    hit.q_from.to_string(),
                    hit.q_to.to_string(),
                    hit.diagonal.to_string(),
                ],
            );
        }
    }
    let mut traps = filter::merge_hits(&sequence.bases, &filter_hits, &filter_params)?;
    if trace::enabled() {
        for (index, trap) in traps.iter().enumerate() {
            trace::emit(
                "04-trapezoids",
                &[
                    index.to_string(),
                    trap.bottom.to_string(),
                    trap.top.to_string(),
                    trap.left.to_string(),
                    trap.right.to_string(),
                ],
            );
        }
    }
    traps.retain(|trap| {
        trap.top
            .checked_sub(trap.bottom)
            .is_some_and(|height| height >= filter_params.word_size)
    });
    let mut hits = dp::align_traps(sequence, &traps, &config.dp)?;
    dp::sort_hits_with_mode(&mut hits, config.sort_mode);
    if trace::enabled() {
        for (index, hit) in hits.iter().enumerate() {
            trace::emit(
                "05-dp-hits",
                &[
                    index.to_string(),
                    hit.a_lo.to_string(),
                    hit.a_hi.to_string(),
                    hit.b_lo.to_string(),
                    hit.b_hi.to_string(),
                    hit.low_diagonal.to_string(),
                    hit.high_diagonal.to_string(),
                    hit.score.to_string(),
                    format!("{:.9}", hit.error),
                ],
            );
        }
    }
    let images = pile::hits_to_images_with_mode(&hits, config.sort_mode);
    if trace::enabled() {
        for (index, image) in images.iter().enumerate() {
            trace::emit(
                "06-images",
                &[
                    index.to_string(),
                    image.lo.to_string(),
                    image.hi.to_string(),
                    image.hit_index.to_string(),
                    u8::from(image.is_a).to_string(),
                ],
            );
        }
    }
    let mut piles = pile::find_piles(&images)?;
    emit_piles("07-piles", &piles.piles);
    pile::expand_piles(&mut piles.piles, sequence, 8);
    emit_piles("08-expanded-piles", &piles.piles);
    if trace::enabled() {
        for (index, (a, b)) in piles
            .hit_to_pile_a
            .iter()
            .zip(&piles.hit_to_pile_b)
            .enumerate()
        {
            trace::emit(
                "09-hit-piles",
                &[index.to_string(), a.to_string(), b.to_string()],
            );
        }
    }
    let edges = graph::make_graph(&piles)?;
    if trace::enabled() {
        for (index, edge) in edges.iter().enumerate() {
            trace::emit(
                "10-edges",
                &[
                    index.to_string(),
                    edge.node1.to_string(),
                    edge.node2.to_string(),
                    u8::from(edge.rev).to_string(),
                ],
            );
        }
    }
    let components =
        graph::connected_components(&edges, config.array.min_array.saturating_sub(2).max(3));
    if trace::enabled() {
        for (index, component) in components.iter().enumerate() {
            for (member_index, member) in component.iter().enumerate() {
                trace::emit(
                    "11-components",
                    &[
                        index.to_string(),
                        member_index.to_string(),
                        member.pile_index.to_string(),
                        u8::from(member.rev).to_string(),
                    ],
                );
            }
        }
    }
    let arrays = array::find_arrays(sequence, &mut piles.piles, &components, &config.array)?;
    if trace::enabled() {
        for array in &arrays {
            trace::emit(
                "16-final-arrays",
                &[
                    array.id.to_string(),
                    array.pos.to_string(),
                    String::from_utf8_lossy(&array.consensus).into_owned(),
                    String::from_utf8_lossy(&array.aligned_consensus).into_owned(),
                    array.elements.len().to_string(),
                ],
            );
        }
    }
    Ok(arrays)
}

fn emit_piles(stage: &str, piles: &[pile::Pile]) {
    if trace::enabled() {
        for (index, pile) in piles.iter().enumerate() {
            trace::emit(
                stage,
                &[
                    index.to_string(),
                    pile.lo.to_string(),
                    pile.hi.to_string(),
                    u8::from(pile.deleted).to_string(),
                    pile.hit_indexes
                        .iter()
                        .map(usize::to_string)
                        .collect::<Vec<_>>()
                        .join(","),
                ],
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_and_short_inputs_have_no_arrays() {
        for fasta in ["", ">empty\n", ">short\nACGTACGT\n"] {
            let sequence = Sequence::from_reader(fasta.as_bytes()).unwrap();
            assert!(detect(&sequence, &Config::default()).unwrap().is_empty());
        }
    }

    #[test]
    fn detects_a_synthetic_array_through_the_complete_pipeline() {
        let repeat = b"ACGTTGCAAGTCCTAGGATCGTACCAAG";
        let mut fasta = b">synthetic\n".to_vec();
        fasta.extend(vec![b'N'; 100]);
        for spacer in [
            b"GATACCGTCTAGTGCATACGATGCCTATCGATACCTGAGT",
            b"CTTGATGACCTAGACTACGGATCATGCTAGTCGATGCTAC",
            b"TAGCTACGATGCTGATCGTAGCTAGCATCGATGCTAGTCA",
        ] {
            fasta.extend_from_slice(repeat);
            fasta.extend_from_slice(spacer);
        }
        fasta.extend(vec![b'N'; 500]);
        let sequence = Sequence::from_reader(fasta.as_slice()).unwrap();
        let arrays = detect(&sequence, &Config::default()).unwrap();
        assert_eq!(arrays.len(), 1);
        assert_eq!(arrays[0].pos, 100);
        assert_eq!(arrays[0].consensus, repeat);
        assert_eq!(arrays[0].elements.len(), 3);
        let directory = std::env::temp_dir().join(format!(
            "pilercr-determinism-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        assert!(!trace::enabled());
        assert!(!directory.exists());
        std::fs::create_dir(&directory).unwrap();
        let render = |arrays: &[Array]| {
            vec![
                report::text(
                    &sequence,
                    arrays,
                    "synthetic",
                    report::CoordinateMode::Corrected,
                    false,
                )
                .unwrap(),
                report::json(&sequence, arrays).unwrap(),
                report::gff(&sequence, arrays).unwrap(),
                report::tsv(&sequence, arrays).unwrap(),
                fasta_out::consensus_fasta(&sequence, arrays, false).unwrap(),
            ]
        };
        let expected = render(&arrays);
        for name in ["first", "second"] {
            let traced = trace::with_directory(&directory.join(name), || {
                detect(&sequence, &Config::default())
            })
            .unwrap();
            assert_eq!(traced, arrays);
            assert_eq!(render(&traced), expected);
        }
        for stage in trace::STAGES {
            let filename = format!("{stage}.tsv");
            let first = std::fs::read(directory.join("first").join(&filename)).unwrap();
            assert_eq!(
                first,
                std::fs::read(directory.join("second").join(&filename)).unwrap()
            );
            assert!(first.starts_with(
                format!("# pilercr-trace-v1\t{stage}\tzero-based;native-endpoints\n").as_bytes()
            ));
        }
        assert!(!trace::enabled());
        std::fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn invalid_configuration_is_rejected_even_for_empty_input() {
        let sequence = Sequence::from_reader(&b""[..]).unwrap();
        let mut config = Config::default();
        config.array.min_array = 0;
        assert!(detect(&sequence, &config).is_err());
        config = Config::default();
        config.dp.min_identity = f64::NAN;
        assert!(detect(&sequence, &config).is_err());
    }
}

use crate::dp::DpHit;
use crate::seq::Sequence;
use crate::sort;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Image {
    pub lo: usize,
    pub hi: usize,
    pub hit_index: usize,
    pub is_a: bool,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Pile {
    pub lo: usize,
    pub hi: usize,
    pub deleted: bool,
    pub hit_indexes: Vec<usize>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PileSet {
    pub piles: Vec<Pile>,
    pub hit_to_pile_a: Vec<usize>,
    pub hit_to_pile_b: Vec<usize>,
}

/// Converts hits to images and sorts image endpoints by `lo`.
///
/// The image order is observable because [`find_piles`] assigns each pile's
/// `hi` from the current image rather than taking a maximum.
pub fn hits_to_images(hits: &[DpHit]) -> Vec<Image> {
    let mut images = Vec::with_capacity(hits.len() * 2);
    for (hit_index, hit) in hits.iter().enumerate() {
        images.push(Image {
            lo: hit.a_lo,
            hi: hit.a_hi,
            hit_index,
            is_a: true,
        });
        images.push(Image {
            lo: hit.b_lo,
            hi: hit.b_hi,
            hit_index,
            is_a: false,
        });
    }
    sort::reference_sort_by_key(&mut images, |image| image.lo);
    images
}

pub fn find_piles(images: &[Image]) -> Result<PileSet, String> {
    if !images.len().is_multiple_of(2) {
        return Err("each hit must have one A image and one B image".into());
    }
    let hit_count = images.len() / 2;
    let mut seen = vec![[false; 2]; hit_count];
    for (index, image) in images.iter().enumerate() {
        if image.lo > image.hi || (index > 0 && image.lo < images[index - 1].lo) {
            return Err("images must have inclusive bounds and be sorted by lo".into());
        }
        let sides = seen
            .get_mut(image.hit_index)
            .ok_or("invalid image hit index")?;
        let side = usize::from(!image.is_a);
        if sides[side] {
            return Err("duplicate hit image".into());
        }
        sides[side] = true;
    }
    let mut result = PileSet {
        piles: Vec::new(),
        hit_to_pile_a: vec![0; hit_count],
        hit_to_pile_b: vec![0; hit_count],
    };
    for image in images {
        if result.piles.last().is_none_or(|pile| image.lo > pile.hi) {
            result.piles.push(Pile {
                lo: image.lo,
                hi: image.hi,
                deleted: false,
                hit_indexes: Vec::new(),
            });
        }
        let pile_index = result.piles.len() - 1;
        let pile = &mut result.piles[pile_index];
        pile.hi = image.hi;
        if image.is_a {
            result.hit_to_pile_a[image.hit_index] = pile_index;
            pile.hit_indexes.push(image.hit_index);
        } else {
            result.hit_to_pile_b[image.hit_index] = pile_index;
        }
    }
    Ok(result)
}

pub fn expand_piles(piles: &mut [Pile], sequence: &Sequence, margin: usize) {
    if margin == 0 {
        return;
    }
    for pile in piles {
        let Some((contig_index, _)) = sequence.global_to_local(pile.lo) else {
            continue;
        };
        let contig = &sequence.contigs[contig_index];
        let contig_hi = contig.start + contig.len - 1;
        pile.lo = pile.lo.saturating_sub(margin).max(contig.start);
        pile.hi = pile.hi.saturating_add(margin).min(contig_hi);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::Contig;

    fn hit(a_lo: usize, a_hi: usize, b_lo: usize, b_hi: usize) -> DpHit {
        DpHit {
            a_lo,
            a_hi,
            b_lo,
            b_hi,
            low_diagonal: 0,
            high_diagonal: 0,
            score: 0,
            error: 0.0,
        }
    }

    #[test]
    fn reference_sort_matches_cpp_tie_order() {
        let mut values: Vec<_> = (0..32).map(|i| (i / 4, i)).collect();
        sort::reference_sort_by_key(&mut values, |value| value.0);
        assert_eq!(
            values.iter().map(|value| value.1).collect::<Vec<_>>(),
            vec![
                1, 2, 3, 0, 4, 5, 6, 7, 9, 11, 10, 8, 12, 13, 14, 15, 16, 19, 18, 17, 20, 21, 22,
                23, 24, 25, 26, 27, 28, 29, 30, 31,
            ]
        );
        let mut equal: Vec<_> = (0..32).map(|i| (0, i)).collect();
        sort::reference_sort_by_key(&mut equal, |value| value.0);
        assert_ne!(
            equal.iter().map(|value| value.1).collect::<Vec<_>>(),
            (0..32).collect::<Vec<_>>()
        );
    }

    #[test]
    fn image_ties_below_partition_threshold_keep_both_sides() {
        let images = hits_to_images(&[hit(5, 20, 5, 6), hit(5, 10, 0, 1)]);
        assert_eq!(
            images
                .iter()
                .map(|i| (i.lo, i.hi, i.hit_index, i.is_a))
                .collect::<Vec<_>>(),
            vec![
                (0, 1, 1, false),
                (5, 20, 0, true),
                (5, 6, 0, false),
                (5, 10, 1, true)
            ]
        );
        let piles = find_piles(&images).unwrap();
        assert_eq!(
            piles.piles.iter().map(|p| (p.lo, p.hi)).collect::<Vec<_>>(),
            vec![(0, 1), (5, 10)]
        );
        assert_eq!(piles.piles[1].hit_indexes, vec![0, 1]);
        assert_eq!(piles.hit_to_pile_a, vec![1, 1]);
        assert_eq!(piles.hit_to_pile_b, vec![1, 0]);
    }

    #[test]
    fn assembly_11_image_ties_preserve_reference_pile_geometry() {
        let coordinates = [
            (22910, 22940, 22971, 23001),
            (23121, 23161, 23187, 23227),
            (23190, 23282, 23246, 23338),
            (23124, 23161, 23246, 23283),
            (23190, 23226, 23302, 23338),
            (23124, 23160, 23302, 23338),
            (23190, 23227, 23368, 23405),
            (23190, 23226, 23434, 23470),
            (23302, 23338, 23434, 23470),
            (23368, 23404, 23434, 23470),
            (23124, 23160, 23434, 23470),
            (24984, 25004, 25262, 25282),
            (29637, 29664, 29666, 29694),
            (42836, 42876, 42879, 42919),
            (77737, 77762, 77821, 77846),
            (84778, 84796, 85158, 85176),
        ];
        let mut hits: Vec<_> = coordinates
            .into_iter()
            .map(|(a, b, c, d)| hit(a, b, c, d))
            .collect();
        crate::dp::sort_hits(&mut hits);
        let piles = find_piles(&hits_to_images(&hits)).unwrap();
        assert_eq!(piles.piles.len(), 17);
        assert_eq!((piles.piles[3].lo, piles.piles[3].hi), (23187, 23283));
        assert_eq!((piles.piles[4].lo, piles.piles[4].hi), (23302, 23338));
    }

    #[test]
    fn inclusive_overlap_and_reference_high_endpoint_replacement() {
        let images = hits_to_images(&[hit(0, 20, 5, 7), hit(7, 8, 9, 9)]);
        let piles = find_piles(&images).unwrap();
        assert_eq!(
            piles.piles.iter().map(|p| (p.lo, p.hi)).collect::<Vec<_>>(),
            vec![(0, 8), (9, 9)]
        );
        assert_eq!(piles.piles[0].hit_indexes, vec![0, 1]);
        assert_eq!(piles.hit_to_pile_a, vec![0, 0]);
        assert_eq!(piles.hit_to_pile_b, vec![0, 1]);
        assert!(piles.piles.iter().all(|p| !p.deleted));
        assert_eq!(find_piles(&[]).unwrap(), PileSet::default());
    }

    #[test]
    fn invalid_images_are_rejected() {
        let images = hits_to_images(&[hit(1, 2, 3, 4)]);
        assert!(find_piles(&images[..1]).is_err());
        let mut bad = images.clone();
        bad.reverse();
        assert!(find_piles(&bad).is_err());
        bad = images.clone();
        bad[0].hi = 0;
        assert!(find_piles(&bad).is_err());
        bad = images.clone();
        bad[0].hit_index = 1;
        assert!(find_piles(&bad).is_err());
        bad = images;
        bad[1].is_a = true;
        assert!(find_piles(&bad).is_err());
    }

    #[test]
    fn expansion_clamps_to_start_contig_and_keeps_membership() {
        let sequence = Sequence {
            bases: vec![b'A'; 20],
            contigs: vec![
                Contig {
                    start: 0,
                    len: 5,
                    label: "a".into(),
                },
                Contig {
                    start: 10,
                    len: 10,
                    label: "b".into(),
                },
            ],
        };
        let mut set = find_piles(&hits_to_images(&[hit(0, 0, 4, 4), hit(10, 10, 19, 19)])).unwrap();
        set.piles[1].deleted = true;
        let before = set.clone();
        expand_piles(&mut set.piles, &sequence, 0);
        assert_eq!(set, before);
        expand_piles(&mut set.piles, &sequence, usize::MAX);
        assert_eq!(
            set.piles.iter().map(|p| (p.lo, p.hi)).collect::<Vec<_>>(),
            vec![(0, 4), (0, 4), (10, 19), (10, 19)]
        );
        assert_eq!(set.hit_to_pile_a, before.hit_to_pile_a);
        assert_eq!(set.hit_to_pile_b, before.hit_to_pile_b);
        assert!(set.piles[1].deleted);
        for (pile, original) in set.piles.iter().zip(&before.piles) {
            assert_eq!(pile.hit_indexes, original.hit_indexes);
        }
        let mut gap = Pile {
            lo: 5,
            hi: 8,
            deleted: false,
            hit_indexes: vec![],
        };
        let original = gap.clone();
        expand_piles(std::slice::from_mut(&mut gap), &sequence, 2);
        assert_eq!(gap, original);
        let mut interior = Pile {
            lo: 13,
            hi: 15,
            deleted: false,
            hit_indexes: vec![],
        };
        expand_piles(std::slice::from_mut(&mut interior), &sequence, 2);
        assert_eq!((interior.lo, interior.hi), (11, 17));
    }
}

use crate::pile::PileSet;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Edge {
    pub node1: usize,
    pub node2: usize,
    pub rev: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ComponentMember {
    pub pile_index: usize,
    pub rev: bool,
}

pub fn make_graph(piles: &PileSet) -> Result<Vec<Edge>, String> {
    if piles.hit_to_pile_a.len() != piles.hit_to_pile_b.len() {
        return Err("A and B pile maps must have equal lengths".into());
    }
    piles
        .hit_to_pile_a
        .iter()
        .zip(&piles.hit_to_pile_b)
        .map(|(&node1, &node2)| {
            if node1 >= piles.piles.len() || node2 >= piles.piles.len() {
                return Err("hit maps to an invalid pile index".into());
            }
            Ok(Edge {
                node1,
                node2,
                rev: false,
            })
        })
        .collect()
}

pub fn connected_components(edges: &[Edge], min_size: usize) -> Vec<Vec<ComponentMember>> {
    let Some(max_index) = edges.iter().flat_map(|edge| [edge.node1, edge.node2]).max() else {
        return Vec::new();
    };
    let mut neighbors = vec![Vec::new(); max_index + 1];
    for edge in edges {
        neighbors[edge.node1].push((edge.node2, edge.rev));
        neighbors[edge.node2].push((edge.node1, edge.rev));
    }
    let mut orientations = vec![None; neighbors.len()];
    let mut components = Vec::new();
    for root in (0..neighbors.len()).rev() {
        if orientations[root].is_some() {
            continue;
        }
        orientations[root] = Some(false);
        let mut pending = vec![root];
        let mut component = Vec::new();
        while let Some(node) = pending.pop() {
            let rev = orientations[node].unwrap();
            for &(neighbor, edge_rev) in &neighbors[node] {
                if orientations[neighbor].is_none() {
                    orientations[neighbor] = Some(rev ^ edge_rev);
                    pending.push(neighbor);
                }
            }
            component.push(ComponentMember {
                pile_index: node,
                rev,
            });
        }
        if component.len() >= min_size {
            component.reverse();
            components.push(component);
        }
    }
    components
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pile::Pile;

    fn members(component: &[ComponentMember]) -> Vec<(usize, bool)> {
        component
            .iter()
            .map(|member| (member.pile_index, member.rev))
            .collect()
    }

    #[test]
    fn graph_preserves_hit_order_duplicates_and_self_edges() {
        let mut piles = PileSet {
            piles: vec![
                Pile {
                    lo: 0,
                    hi: 1,
                    deleted: true,
                    hit_indexes: vec![]
                };
                3
            ],
            hit_to_pile_a: vec![2, 0, 0, 1],
            hit_to_pile_b: vec![0, 2, 2, 1],
        };
        let edges = make_graph(&piles).unwrap();
        assert_eq!(
            edges,
            vec![
                Edge {
                    node1: 2,
                    node2: 0,
                    rev: false
                },
                Edge {
                    node1: 0,
                    node2: 2,
                    rev: false
                },
                Edge {
                    node1: 0,
                    node2: 2,
                    rev: false
                },
                Edge {
                    node1: 1,
                    node2: 1,
                    rev: false
                },
            ]
        );
        let components = connected_components(&edges, 1);
        assert_eq!(
            components.iter().map(|c| members(c)).collect::<Vec<_>>(),
            vec![vec![(0, false), (2, false)], vec![(1, false)]]
        );
        piles.hit_to_pile_b.pop();
        assert!(make_graph(&piles).is_err());
        piles.hit_to_pile_a.pop();
        piles.hit_to_pile_a[0] = 3;
        assert!(make_graph(&piles).is_err());
        assert!(make_graph(&PileSet::default()).unwrap().is_empty());
    }

    #[test]
    fn components_are_undirected_with_reference_order_and_reversal() {
        let edges = vec![
            Edge {
                node1: 3,
                node2: 0,
                rev: true,
            },
            Edge {
                node1: 3,
                node2: 2,
                rev: false,
            },
            Edge {
                node1: 2,
                node2: 1,
                rev: true,
            },
        ];
        let components = connected_components(&edges, 4);
        assert_eq!(
            members(&components[0]),
            vec![(0, true), (1, true), (2, false), (3, false)]
        );
        let reversed: Vec<_> = edges
            .iter()
            .map(|e| Edge {
                node1: e.node2,
                node2: e.node1,
                rev: e.rev,
            })
            .collect();
        assert_eq!(connected_components(&reversed, 4), components);
        assert!(connected_components(&edges, 5).is_empty());
    }

    #[test]
    fn sparse_indices_include_isolates_and_empty_edges_have_no_nodes() {
        let edges = [Edge {
            node1: 1,
            node2: 3,
            rev: false,
        }];
        let components = connected_components(&edges, 0);
        assert_eq!(
            components.iter().map(|c| members(c)).collect::<Vec<_>>(),
            vec![
                vec![(1, false), (3, false)],
                vec![(2, false)],
                vec![(0, false)]
            ]
        );
        assert_eq!(connected_components(&edges, 2).len(), 1);
        assert!(connected_components(&[], 0).is_empty());
    }

    #[test]
    fn conflicting_cycles_keep_first_discovered_orientation() {
        let edges = [
            Edge {
                node1: 2,
                node2: 0,
                rev: false,
            },
            Edge {
                node1: 2,
                node2: 1,
                rev: false,
            },
            Edge {
                node1: 1,
                node2: 0,
                rev: true,
            },
            Edge {
                node1: 2,
                node2: 2,
                rev: true,
            },
        ];
        assert_eq!(
            members(&connected_components(&edges, 1)[0]),
            vec![(0, false), (1, false), (2, false)]
        );
    }
}

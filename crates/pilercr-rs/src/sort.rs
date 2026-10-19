//! Tie-breaking-compatible port of modern libstdc++ `std::sort`.
//!
//! PILER-CR sorts DP hits and image endpoints with an unstable comparator
//! (`a.lo < b.lo`) and then assigns pile `Hi` from the image order, so the
//! permutation of equal keys is observable in the report. This implementation
//! reproduces the modern libstdc++ introsort used by the project's acceptance
//! target.

/// Sorts `values` like modern libstdc++ `std::sort`, using `key` for ordering.
pub(crate) fn reference_sort_by_key<T: Clone, K: Ord>(values: &mut [T], key: impl Fn(&T) -> K) {
    fn adjust<T: Clone, K: Ord>(
        values: &mut [T],
        mut hole: usize,
        value: T,
        key: &impl Fn(&T) -> K,
    ) {
        let top = hole;
        let len = values.len();
        let mut child = hole;
        while child < (len - 1) / 2 {
            child = 2 * (child + 1);
            if key(&values[child]) < key(&values[child - 1]) {
                child -= 1;
            }
            values[hole] = values[child].clone();
            hole = child;
        }
        if len.is_multiple_of(2) && child == (len - 2) / 2 {
            child = 2 * child + 1;
            values[hole] = values[child].clone();
            hole = child;
        }
        while hole > top {
            let parent = (hole - 1) / 2;
            if key(&values[parent]) >= key(&value) {
                break;
            }
            values[hole] = values[parent].clone();
            hole = parent;
        }
        values[hole] = value;
    }
    fn partition<T: Clone, K: Ord>(values: &mut [T], mut depth: u32, key: &impl Fn(&T) -> K) {
        let mut end = values.len();
        while end > 16 {
            if depth == 0 {
                for parent in (0..=(end - 2) / 2).rev() {
                    let value = values[parent].clone();
                    adjust(&mut values[..end], parent, value, key);
                }
                for last in (1..end).rev() {
                    let value = values[last].clone();
                    values[last] = values[0].clone();
                    adjust(&mut values[..last], 0, value, key);
                }
                return;
            }
            depth -= 1;
            let (a, b, c) = (1, end / 2, end - 1);
            let median = if key(&values[a]) < key(&values[b]) {
                if key(&values[b]) < key(&values[c]) {
                    b
                } else if key(&values[a]) < key(&values[c]) {
                    c
                } else {
                    a
                }
            } else if key(&values[a]) < key(&values[c]) {
                a
            } else if key(&values[b]) < key(&values[c]) {
                c
            } else {
                b
            };
            values.swap(0, median);
            let (mut left, mut right) = (1, end);
            loop {
                while key(&values[left]) < key(&values[0]) {
                    left += 1;
                }
                right -= 1;
                while key(&values[0]) < key(&values[right]) {
                    right -= 1;
                }
                if left >= right {
                    break;
                }
                values.swap(left, right);
                left += 1;
            }
            partition(&mut values[left..end], depth, key);
            end = left;
        }
    }
    if values.len() < 2 {
        return;
    }
    partition(values, 2 * values.len().ilog2(), &key);
    for i in 1..values.len() {
        let mut j = i;
        while j > 0 && key(&values[j]) < key(&values[j - 1]) {
            values.swap(j, j - 1);
            j -= 1;
        }
    }
}

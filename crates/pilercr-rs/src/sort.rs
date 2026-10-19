//! Tie-breaking-compatible ports of libstdc++ `std::sort` implementations.
//!
//! PILER-CR sorts DP hits and image endpoints with an unstable comparator
//! (`a.lo < b.lo`) and then assigns pile `Hi` from the image order, so the
//! permutation of equal keys is observable in the report. The default mode
//! reproduces the modern libstdc++ introsort used by the project's acceptance
//! target. [`SortMode::Legacy`] reproduces GCC 4.4.3's libstdc++ introsort,
//! which is what the shipped 2011 executable used.

/// Selects the `std::sort` tie-breaking behaviour used during detection.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SortMode {
    /// Modern libstdc++ introsort (the default acceptance target).
    #[default]
    Modern,
    /// GCC 4.4.3 libstdc++ introsort, matching the shipped 2011 executable.
    Legacy,
}

pub(crate) fn sort_by_key<T: Copy, K: Ord>(
    mode: SortMode,
    values: &mut [T],
    key: impl Fn(&T) -> K,
) {
    match mode {
        SortMode::Modern => reference_sort_by_key(values, key),
        SortMode::Legacy => legacy_sort_by_key(values, key),
    }
}

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

const S_THRESHOLD: usize = 16;

/// libstdc++ 4.4.3 `std::__lg`: index of the highest set bit.
fn lg_443(n: usize) -> usize {
    if n == 0 {
        return 0;
    }
    usize::BITS as usize - 1 - n.leading_zeros() as usize
}

/// Sorts `values` like GCC 4.4.3 libstdc++ `std::sort`, using `key`.
///
/// GCC 4.4.3 differs from modern libstdc++ in two observable ways:
/// the pivot is the median *value* of `first`, `first + n/2` and `last - 1`
/// (passed by value, with no median-to-first swap), and
/// `__unguarded_partition` scans `[first, last)`.
pub(crate) fn legacy_sort_by_key<T: Copy, K: Ord>(values: &mut [T], key: impl Fn(&T) -> K) {
    let n = values.len();
    if n == 0 {
        return;
    }
    introsort_loop(values, 0, n, lg_443(n) * 2, &key);
    final_insertion_sort(values, 0, n, &key);
}

fn introsort_loop<T: Copy, K: Ord>(
    v: &mut [T],
    first: usize,
    mut last: usize,
    mut depth_limit: usize,
    key: &impl Fn(&T) -> K,
) {
    while last - first > S_THRESHOLD {
        if depth_limit == 0 {
            partial_sort(v, first, last, key);
            return;
        }
        depth_limit -= 1;
        let pivot = median(&v[first], &v[first + (last - first) / 2], &v[last - 1], key);
        let cut = unguarded_partition(v, first, last, &pivot, key);
        introsort_loop(v, cut, last, depth_limit, key);
        last = cut;
    }
}

/// libstdc++ 4.4.3 `std::__median`: median of three values under `less`.
fn median<T: Copy, K: Ord>(a: &T, b: &T, c: &T, key: &impl Fn(&T) -> K) -> T {
    if key(a) < key(b) {
        if key(b) < key(c) {
            *b
        } else if key(a) < key(c) {
            *c
        } else {
            *a
        }
    } else if key(a) < key(c) {
        *a
    } else if key(b) < key(c) {
        *c
    } else {
        *b
    }
}

fn unguarded_partition<T: Copy, K: Ord>(
    v: &mut [T],
    mut first: usize,
    mut last: usize,
    pivot: &T,
    key: &impl Fn(&T) -> K,
) -> usize {
    // `first` advances in the scan loop below.
    loop {
        while key(&v[first]) < key(pivot) {
            first += 1;
        }
        last -= 1;
        while key(pivot) < key(&v[last]) {
            last -= 1;
        }
        if first >= last {
            return first;
        }
        v.swap(first, last);
        first += 1;
    }
}

fn final_insertion_sort<T: Copy, K: Ord>(
    v: &mut [T],
    first: usize,
    last: usize,
    key: &impl Fn(&T) -> K,
) {
    if last - first > S_THRESHOLD {
        insertion_sort(v, first, first + S_THRESHOLD, key);
        unguarded_insertion_sort(v, first + S_THRESHOLD, last, key);
    } else {
        insertion_sort(v, first, last, key);
    }
}

fn insertion_sort<T: Copy, K: Ord>(v: &mut [T], first: usize, last: usize, key: &impl Fn(&T) -> K) {
    if first == last {
        return;
    }
    for i in (first + 1)..last {
        let value = v[i];
        if key(&value) < key(&v[first]) {
            v.copy_within(first..i, first + 1);
            v[first] = value;
        } else {
            unguarded_linear_insert(v, i, value, key);
        }
    }
}

fn unguarded_insertion_sort<T: Copy, K: Ord>(
    v: &mut [T],
    first: usize,
    last: usize,
    key: &impl Fn(&T) -> K,
) {
    for i in first..last {
        let value = v[i];
        unguarded_linear_insert(v, i, value, key);
    }
}

fn unguarded_linear_insert<T: Copy, K: Ord>(
    v: &mut [T],
    mut last: usize,
    value: T,
    key: &impl Fn(&T) -> K,
) {
    let mut next = last - 1;
    while key(&value) < key(&v[next]) {
        v[last] = v[next];
        last = next;
        if next == 0 {
            break;
        }
        next -= 1;
    }
    v[last] = value;
}

// --- heap fallback (`std::__partial_sort`) ---

fn partial_sort<T: Copy, K: Ord>(v: &mut [T], first: usize, last: usize, key: &impl Fn(&T) -> K) {
    // With `middle == last`, `__heap_select` reduces to `__make_heap`.
    make_heap(v, first, last, key);
    sort_heap(v, first, last, key);
}

fn make_heap<T: Copy, K: Ord>(v: &mut [T], first: usize, last: usize, key: &impl Fn(&T) -> K) {
    let len = last - first;
    if len < 2 {
        return;
    }
    let mut parent = (len - 2) / 2;
    loop {
        let value = v[first + parent];
        adjust_heap(v, first, parent, len, value, key);
        if parent == 0 {
            break;
        }
        parent -= 1;
    }
}

/// libstdc++ `std::__adjust_heap`. `value` is the displaced element that is
/// sifted back into the heap; it must be captured before `hole` is overwritten.
fn adjust_heap<T: Copy, K: Ord>(
    v: &mut [T],
    first: usize,
    mut hole: usize,
    len: usize,
    value: T,
    key: &impl Fn(&T) -> K,
) {
    let top = hole;
    let mut second_child = 2 * (hole + 1);
    while second_child < len {
        if key(&v[first + second_child]) < key(&v[first + second_child - 1]) {
            second_child -= 1;
        }
        v[first + hole] = v[first + second_child];
        hole = second_child;
        second_child = 2 * (second_child + 1);
    }
    if second_child == len {
        v[first + hole] = v[first + second_child - 1];
        hole = second_child - 1;
    }
    push_heap(v, first, hole, top, value, key);
}

fn push_heap<T: Copy, K: Ord>(
    v: &mut [T],
    first: usize,
    mut hole: usize,
    top: usize,
    value: T,
    key: &impl Fn(&T) -> K,
) {
    let mut parent = hole.saturating_sub(1) / 2;
    while hole > top && key(&v[first + parent]) < key(&value) {
        v[first + hole] = v[first + parent];
        hole = parent;
        parent = hole.saturating_sub(1) / 2;
    }
    v[first + hole] = value;
}

fn sort_heap<T: Copy, K: Ord>(v: &mut [T], first: usize, mut last: usize, key: &impl Fn(&T) -> K) {
    while last - first > 1 {
        last -= 1;
        v.swap(first, last);
        let value = v[first];
        adjust_heap(v, first, 0, last - first, value, key);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ids<T>(values: &[(u32, T)]) -> Vec<T>
    where
        T: Copy,
    {
        values.iter().map(|value| value.1).collect()
    }

    #[test]
    fn legacy_equal_keys_pin_gcc_443_permutation() {
        // For n > 16, GCC 4.4.3's first partition reverses equal keys around
        // the pivot; insertion sort preserves that order. This is distinct
        // from modern libstdc++, which swaps the median to the front first.
        let mut equal: Vec<_> = (0..32).map(|id| (0u32, id)).collect();
        legacy_sort_by_key(&mut equal, |value| value.0);
        assert_eq!(ids(&equal), (0..32).rev().collect::<Vec<_>>());

        let mut seventeen: Vec<_> = (0..17).map(|id| (0u32, id)).collect();
        legacy_sort_by_key(&mut seventeen, |value| value.0);
        assert_eq!(ids(&seventeen), (0..17).rev().collect::<Vec<_>>());

        // At or below the introsort threshold, insertion sort is stable.
        let mut small: Vec<_> = (0..16).map(|id| (0u32, id)).collect();
        legacy_sort_by_key(&mut small, |value| value.0);
        assert_eq!(ids(&small), (0..16).collect::<Vec<_>>());
    }

    #[test]
    fn legacy_mixed_duplicate_keys_pin_gcc_443_permutation() {
        // Reference permutation from the GCC 4.4.3 algorithm for this key
        // pattern. Ties are not stable and follow the old introsort exactly.
        let keys = [
            5, 1, 3, 5, 1, 1, 7, 2, 2, 1, 7, 5, 6, 5, 4, 1, 7, 0, 6, 7, 6, 7, 2, 7, 7, 0, 1, 1, 4,
            5, 6, 3,
        ];
        let mut values: Vec<_> = keys
            .into_iter()
            .enumerate()
            .map(|(id, key)| (key, id))
            .collect();
        legacy_sort_by_key(&mut values, |value| value.0);
        assert_eq!(
            ids(&values),
            vec![
                17, 25, 1, 15, 4, 5, 26, 27, 9, 22, 8, 7, 2, 31, 28, 14, 29, 13, 11, 3, 0, 18, 20,
                12, 30, 16, 19, 21, 23, 24, 10, 6,
            ]
        );
    }

    #[test]
    fn legacy_and_modern_differ_on_equal_keys() {
        let original: Vec<_> = (0..32).map(|id| (0u32, id)).collect();
        let mut modern = original.clone();
        let mut legacy = original;
        sort_by_key(SortMode::Modern, &mut modern, |value| value.0);
        sort_by_key(SortMode::Legacy, &mut legacy, |value| value.0);
        assert_ne!(ids(&modern), ids(&legacy));
        assert_eq!(ids(&legacy), (0..32).rev().collect::<Vec<_>>());
    }

    #[test]
    fn legacy_preserves_multiset_and_sorts_duplicate_keys() {
        let mut state = 0x1234_5678_9abc_def0u64;
        let mut values = Vec::new();
        for id in 0..500 {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            let key = (state >> 32) as u32 % 23;
            values.push((key, id));
        }
        let mut sorted = values;
        legacy_sort_by_key(&mut sorted, |value| value.0);

        let mut actual_ids: Vec<_> = ids(&sorted);
        actual_ids.sort_unstable();
        assert_eq!(actual_ids, (0..500).collect::<Vec<_>>());
        assert!(
            sorted.windows(2).all(|pair| pair[0].0 <= pair[1].0),
            "legacy sorter must order by key"
        );
    }

    #[test]
    fn legacy_sorts_distinct_keys_ascending() {
        let mut values: Vec<_> = (0..100).rev().map(|id| (id, id)).collect();
        legacy_sort_by_key(&mut values, |value| value.0);
        assert_eq!(ids(&values), (0..100).collect::<Vec<_>>());
    }

    #[test]
    fn legacy_heap_fallback_sorts_and_does_not_underflow() {
        // `__partial_sort` is only reached when the introsort depth limit is
        // exhausted, so exercise it directly rather than relying on a case
        // that happens to trigger the fallback.
        let mut values: Vec<_> = (0..40).rev().map(|id| (id, id)).collect();
        let len = values.len();
        partial_sort(&mut values, 0, len, &|value| value.0);
        assert_eq!(ids(&values), (0..40).collect::<Vec<_>>());

        // Equal keys still follow the old heap permutation.
        let mut equal: Vec<_> = (0..5).map(|id| (0u32, id)).collect();
        let len = equal.len();
        partial_sort(&mut equal, 0, len, &|value| value.0);
        assert_eq!(ids(&equal), vec![3, 4, 1, 0, 2]);
    }
}

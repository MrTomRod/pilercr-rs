/// Consensus string and the optional conserved span it was taken from.
pub type ConservedConsensus = (Vec<u8>, Option<(usize, usize)>);

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Msa {
    pub rows: Vec<Vec<u8>>,
    pub consensus: Vec<u8>,
}

#[derive(Clone, Default)]
struct Position {
    counts: [f32; 4],
    scores: [f32; 4],
    transitions: [f32; 4],
    open: f32,
    close: f32,
}

struct Profile {
    positions: Vec<Position>,
    rows: Vec<(usize, Vec<u8>)>,
}

fn letter(c: u8) -> usize {
    match c.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        _ => 3,
    }
}

fn gap(c: u8) -> bool {
    c == b'-' || c == b'.'
}

fn score_positions(positions: &mut [Position]) {
    for i in 0..positions.len() {
        let close = if i + 1 < positions.len() {
            positions[i + 1].transitions[2]
        } else {
            positions[i].transitions[3] + positions[i].transitions[1]
        };
        let p = &mut positions[i];
        p.open = ((1.0 - p.transitions[1] as f64) * -400.0 / 2.0) as f32;
        p.close = ((1.0 - close as f64) * -400.0 / 2.0) as f32;
        for a in 0..4 {
            p.scores[a] = 0.0;
            for b in 0..4 {
                p.scores[a] += p.counts[b] * if a == b { 50.0 } else { -100.0 };
            }
        }
    }
}

fn leaf(id: usize, row: Vec<u8>) -> Profile {
    let mut positions = Vec::with_capacity(row.len());
    for (i, &c) in row.iter().enumerate() {
        let mut p = Position::default();
        if !gap(c) {
            if c.eq_ignore_ascii_case(&b'N') {
                p.counts = [0.05; 4];
            } else {
                p.counts[letter(c)] = 1.0;
            }
        }
        let prev_gap = i > 0 && gap(row[i - 1]);
        p.transitions[match (prev_gap, gap(c)) {
            (false, false) => 0,
            (false, true) => 1,
            (true, false) => 2,
            (true, true) => 3,
        }] = 1.0;
        positions.push(p);
    }
    score_positions(&mut positions);
    Profile {
        positions,
        rows: vec![(id, row)],
    }
}

fn tuple_counts(row: &[u8]) -> [u8; 256] {
    let mut counts = [0u8; 256];
    for start in 0..row.len().saturating_sub(5) {
        let mut tuple = 0;
        for &c in &row[start..start + 4] {
            tuple = tuple * 4 + letter(c);
        }
        counts[tuple] = counts[tuple].wrapping_add(1);
    }
    counts
}

fn distances(repeats: &[Vec<u8>]) -> Vec<Vec<f32>> {
    let counts: Vec<_> = repeats.iter().map(|r| tuple_counts(r)).collect();
    let totals: Vec<u32> = counts
        .iter()
        .map(|c| c.iter().map(|&x| x as u32).sum())
        .collect();
    let mut distances = vec![vec![0.0; repeats.len()]; repeats.len()];
    for i in 1..repeats.len() {
        for j in 0..i {
            let shared: u32 = counts[i]
                .iter()
                .zip(&counts[j])
                .map(|(&a, &b)| a.min(b) as u32)
                .sum();
            let a = totals[i].max(1) as f64;
            let b = totals[j].max(1) as f64;
            let distance = ((a - shared as f64) / a).min((b - shared as f64) / b) as f32;
            distances[i][j] = distance;
            distances[j][i] = distance;
        }
    }
    distances
}

fn guide_tree(repeats: &[Vec<u8>]) -> Vec<(usize, usize)> {
    let n = repeats.len();
    let mut distances = distances(repeats);
    let mut nodes: Vec<_> = (0..n).map(Some).collect();
    let mut nearest = vec![0; n];
    let mut minima = vec![f32::MAX; n];
    for i in 1..n {
        for j in 0..i {
            let d = distances[i][j];
            if d < minima[i] {
                minima[i] = d;
                nearest[i] = j;
            }
            if d < minima[j] {
                minima[j] = d;
                nearest[j] = i;
            }
        }
    }
    let mut joins = Vec::with_capacity(n.saturating_sub(1));
    for step in 0..n.saturating_sub(1) {
        let mut left = 0;
        let mut best = f32::MAX;
        for j in 0..n {
            if nodes[j].is_some() && minima[j] < best {
                left = j;
                best = minima[j];
            }
        }
        let right = nearest[left];
        let mut new_min = f32::MAX;
        let mut new_nearest = 0;
        for j in 0..n {
            if j == left || j == right || nodes[j].is_none() {
                continue;
            }
            let d = distances[left][j].min(distances[right][j]);
            if nearest[j] == right {
                nearest[j] = left;
            }
            distances[left][j] = d;
            distances[j][left] = d;
            if d < new_min {
                new_min = d;
                new_nearest = j;
            }
        }
        joins.push((nodes[left].unwrap(), nodes[right].unwrap()));
        nodes[left] = Some(n + step);
        nodes[right] = None;
        nearest[left] = new_nearest;
        minima[left] = new_min;
    }
    joins
}

fn column_score(a: &Position, b: &Position) -> f32 {
    let mut order = [0, 1, 2, 3];
    order.sort_by(|&i, &j| a.counts[j].partial_cmp(&a.counts[i]).unwrap());
    let mut score = 0.0;
    for i in order {
        if a.counts[i] == 0.0 {
            break;
        }
        score += a.counts[i] * b.scores[i];
    }
    score
}

fn alignment_path(a: &mut [Position], b: &mut [Position]) -> Result<Vec<usize>, String> {
    let n = a.len();
    let m = b.len();
    a[0].open *= -1.0;
    b[0].open *= -1.0;
    if n > 1 {
        a[n - 1].close *= -1.0;
    }
    if m > 1 {
        b[m - 1].close *= -1.0;
    }
    let width = m + 1;
    let cells = (n + 1)
        .checked_mul(width)
        .ok_or("MSA matrix size overflow")?;
    let mut trace = Vec::new();
    trace
        .try_reserve_exact(cells)
        .map_err(|_| "Cannot allocate MSA traceback")?;
    trace.resize(cells, [0u8; 3]);
    let neg = -1e37f32;
    let mut previous = vec![[neg; 3]; width];
    previous[0][0] = 0.0;
    for j in 1..=m {
        previous[j][2] = b[0].open + (j - 1) as f32 * -25.0;
        trace[j][2] = if j == 1 { 0 } else { 2 };
    }
    for i in 1..=n {
        let mut current = vec![[neg; 3]; width];
        if i < n {
            current[0][1] = a[0].open + (i - 1) as f32 * -25.0;
            trace[i * width][1] = if i == 1 { 0 } else { 1 };
        }
        for j in 1..=m {
            let k = i * width + j;
            let score = column_score(&a[i - 1], &b[j - 1]);
            if i == 1 {
                current[j][0] = if j == 1 {
                    score
                } else {
                    score + b[0].open + (j - 2) as f32 * -25.0 + b[j - 2].close
                };
                trace[k][0] = if j == 1 { 0 } else { 2 };
            } else if j == 1 {
                current[j][0] = score + a[0].open + (i - 2) as f32 * -25.0 + a[i - 2].close;
                trace[k][0] = 1;
            } else {
                let mm = previous[j - 1][0];
                let dm = previous[j - 1][1] + a[i - 2].close;
                let im = previous[j - 1][2] + b[j - 2].close;
                let state = if mm >= dm && mm >= im {
                    0
                } else if dm >= mm && dm >= im {
                    1
                } else {
                    2
                };
                current[j][0] = score + [mm, dm, im][state];
                trace[k][0] = state as u8;
            }
            if i == n && j == 1 {
                if n == 1 {
                    current[j][0] = score + a[0].open + a[0].close;
                }
                trace[k][0] = 1;
            }
            let dd = previous[j][1] - 25.0;
            let md = previous[j][0] + a[i - 1].open;
            let state = if dd > md { 1 } else { 0 };
            current[j][1] = if state == 1 { dd } else { md };
            trace[k][1] = state;
            let ii = current[j - 1][2] - 25.0;
            let mi = current[j - 1][0] + b[j - 1].open;
            let state = if mi >= ii { 0 } else { 2 };
            current[j][2] = if state == 0 { mi } else { ii };
            trace[k][2] = state;
        }
        previous = current;
    }
    let mut state = 0;
    if previous[m][1] > previous[m][state] {
        state = 1;
    }
    if previous[m][2] > previous[m][state] {
        state = 2;
    }
    let (mut i, mut j) = (n, m);
    let mut path = Vec::with_capacity(n + m);
    while i != 0 || j != 0 {
        if (state != 2 && i == 0) || (state != 1 && j == 0) {
            return Err("Reference MSA traceback reached an invalid boundary".into());
        }
        path.push(state);
        let next = trace[i * width + j][state] as usize;
        if state != 2 {
            i -= 1;
        }
        if state != 1 {
            j -= 1;
        }
        state = next;
    }
    path.reverse();
    Ok(path)
}

fn transition(p: &Position, consume: bool, consumed_previous: bool) -> [f32; 4] {
    let [ll, lg, gl, gg] = p.transitions;
    match (consume, consumed_previous) {
        (true, true) => [ll, lg, gl, gg],
        (true, false) => [0.0, 0.0, ll + gl, lg + gg],
        (false, true) => [0.0, ll + gl, 0.0, lg + gg],
        (false, false) => [0.0, 0.0, 0.0, 1.0],
    }
}

fn merge(mut a: Profile, mut b: Profile) -> Result<Profile, String> {
    let path = alignment_path(&mut a.positions, &mut b.positions)?;
    let start = Position {
        transitions: [1.0, 0.0, 0.0, 0.0],
        ..Position::default()
    };
    let mut positions = Vec::with_capacity(path.len());
    let (mut i, mut j, mut previous) = (0, 0, 0);
    for &state in &path {
        let ca = state != 2;
        let cb = state != 1;
        if ca {
            i += 1;
        }
        if cb {
            j += 1;
        }
        let pa = if i == 0 { &start } else { &a.positions[i - 1] };
        let pb = if j == 0 { &start } else { &b.positions[j - 1] };
        let ta = transition(pa, ca, previous != 2);
        let tb = transition(pb, cb, previous != 1);
        let mut p = Position::default();
        for k in 0..4 {
            p.counts[k] = if ca { 0.5 * pa.counts[k] } else { 0.0 }
                + if cb { 0.5 * pb.counts[k] } else { 0.0 };
            p.transitions[k] = 0.5 * ta[k] + 0.5 * tb[k];
        }
        positions.push(p);
        previous = state;
    }
    score_positions(&mut positions);
    for (side, rows) in [(2, &mut a.rows), (1, &mut b.rows)] {
        for (_, row) in rows {
            let mut index = 0;
            let mut out = Vec::with_capacity(path.len());
            for &state in &path {
                if state == side {
                    out.push(b'-');
                } else {
                    out.push(row[index]);
                    index += 1;
                }
            }
            *row = out;
        }
    }
    a.rows.extend(b.rows);
    Ok(Profile {
        positions,
        rows: a.rows,
    })
}

fn column_consensus(rows: &[Vec<u8>], col: usize) -> (u8, usize) {
    let mut counts = [0; 5];
    for row in rows {
        match row[col].to_ascii_uppercase() {
            b'A' => counts[0] += 1,
            b'C' => counts[1] += 1,
            b'G' => counts[2] += 1,
            b'T' => counts[3] += 1,
            b'-' => counts[4] += 1,
            _ => (),
        }
    }
    let mut best = 0;
    for i in 1..5 {
        if counts[i] > counts[best] {
            best = i;
        }
    }
    (b"ACGT-"[best], counts[best])
}

impl Msa {
    pub fn ungapped_consensus(&self) -> Vec<u8> {
        self.consensus
            .iter()
            .copied()
            .filter(|&c| !gap(c))
            .collect()
    }

    pub fn match_rows(&self) -> Vec<Vec<u8>> {
        self.rows
            .iter()
            .map(|row| {
                row.iter()
                    .zip(&self.consensus)
                    .map(|(&c, &cons)| if c == cons && !gap(c) { b'.' } else { c })
                    .collect()
            })
            .collect()
    }

    pub fn conserved_consensus(&self, min_conservation: f64) -> Result<ConservedConsensus, String> {
        if !min_conservation.is_finite() || !(0.0..=1.0).contains(&min_conservation) {
            return Err("Conservation threshold must be between zero and one".into());
        }
        let cols = self.rows.first().map_or(0, Vec::len);
        if self.rows.iter().any(|r| r.len() != cols) {
            return Err("MSA rows have unequal lengths".into());
        }
        if cols == 0 {
            return Ok((Vec::new(), None));
        }
        let calls: Vec<_> = (0..cols).map(|c| column_consensus(&self.rows, c)).collect();
        let mut conserved: Vec<_> = calls
            .iter()
            .map(|&(_, count)| count as f64 / self.rows.len() as f64 >= min_conservation)
            .collect();
        // Indexes neighbours (col-2..=col+2) on purpose; an iterator form obscures that.
        #[allow(clippy::needless_range_loop)]
        for col in 2..cols.saturating_sub(2) {
            if !conserved[col]
                && conserved[col - 1]
                && conserved[col - 2]
                && conserved[col + 1]
                && conserved[col + 2]
            {
                conserved[col] = true;
            }
        }
        let (mut best_start, mut best_end) = (0, 0);
        let mut start = None;
        #[allow(clippy::needless_range_loop)] // reads conserved[col] alongside the span logic
        for col in 0..=cols {
            if col < cols && conserved[col] {
                if start.is_none() {
                    start = Some(col);
                }
            } else if let Some(s) = start.take()
                && col - s > best_end - best_start + 1
            {
                best_start = s;
                best_end = col - 1;
            }
        }
        let mut sequence: Vec<_> = calls[best_start..=best_end]
            .iter()
            .map(|&(c, _)| c)
            .filter(|&c| c != b'-')
            .collect();
        let mut next = b'-';
        for row in &self.rows {
            for (col, &c) in row.iter().enumerate().take(cols).skip(best_end + 1) {
                if c == b'-' {
                    continue;
                }
                if col == best_end + 1 {
                    next = c;
                } else if next != c {
                    return Ok((sequence, Some((best_start, best_end))));
                }
                break;
            }
        }
        if next != b'-' {
            sequence.push(next);
        }
        Ok((sequence, Some((best_start, best_end))))
    }
}

pub fn align(repeats: &[Vec<u8>]) -> Result<Msa, String> {
    let call_index = crate::trace::enabled().then(|| crate::trace::next_id("msa"));
    if repeats.is_empty() {
        return Ok(Msa::default());
    }
    if repeats.len() > 4096 {
        return Err("Repeat count exceeds the reference MSA limit of 4096".into());
    }
    for (id, row) in repeats.iter().enumerate() {
        if row.is_empty() {
            return Err(format!("Repeat {id} is empty"));
        }
        if row.len() >= 100000 {
            return Err(format!(
                "Repeat {id} exceeds the reference MSA length limit"
            ));
        }
        if row.iter().any(|c| {
            !matches!(
                c.to_ascii_uppercase(),
                b'A' | b'C' | b'G' | b'T' | b'U' | b'N' | b'-' | b'.'
            )
        }) {
            return Err(format!("Repeat {id} contains an unsupported nucleotide"));
        }
    }
    let joins = guide_tree(repeats);
    let mut profiles: Vec<Option<Profile>> = repeats
        .iter()
        .enumerate()
        .map(|(i, row)| Some(leaf(i, row.clone())))
        .collect();
    for (left, right) in joins {
        let a = profiles[left].take().unwrap();
        let b = profiles[right].take().unwrap();
        profiles.push(Some(merge(a, b)?));
    }
    let mut profile = profiles.pop().unwrap().unwrap();
    profile.rows.sort_by_key(|(id, _)| *id);
    let rows: Vec<_> = profile.rows.into_iter().map(|(_, row)| row).collect();
    let consensus = (0..rows[0].len())
        .map(|col| column_consensus(&rows, col).0)
        .collect();
    if let Some(call_index) = call_index {
        for (row_index, row) in rows.iter().enumerate() {
            crate::trace::emit(
                "13-msa",
                &[
                    call_index.to_string(),
                    row_index.to_string(),
                    String::from_utf8_lossy(row).into_owned(),
                ],
            );
        }
    }
    Ok(Msa { rows, consensus })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(rows: &[&[u8]]) -> Msa {
        align(&rows.iter().map(|r| r.to_vec()).collect::<Vec<_>>()).unwrap()
    }

    #[test]
    fn trace_rows_follow_input_order_and_call_ids_reset() {
        let directory = std::env::temp_dir().join(format!(
            "pilercr-msa-trace-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        for _ in 0..2 {
            crate::trace::with_directory(&directory, || {
                align(&[])?;
                let alignment = align(&[b"ACGT".to_vec(), b"ACGT".to_vec()])?;
                assert_eq!(alignment.rows, vec![b"ACGT".to_vec(); 2]);
                align(&[b"TGCA".to_vec()])?;
                Ok(())
            })
            .unwrap();
            let trace = std::fs::read_to_string(directory.join("13-msa.tsv")).unwrap();
            assert_eq!(
                trace.lines().skip(1).collect::<Vec<_>>(),
                ["1\t0\tACGT", "1\t1\tACGT", "2\t0\tTGCA"]
            );
            std::fs::remove_dir_all(&directory).unwrap();
        }
        assert!(!crate::trace::enabled());
    }

    #[test]
    fn exact_and_variable_repeats() {
        let exact = run(&[b"ACGTACGT", b"ACGTACGT", b"ACGTACGT"]);
        assert_eq!(exact.consensus, b"ACGTACGT");
        assert_eq!(exact.rows, vec![b"ACGTACGT".to_vec(); 3]);
        let variable = run(&[b"ACGTACGT", b"ACGTTCGT", b"ACGTACGT"]);
        assert_eq!(variable.consensus, b"ACGTACGT");
        assert_eq!(variable.match_rows()[1], b"....T...");
    }

    #[test]
    fn gaps_preserve_rows() {
        let inputs: &[&[u8]] = &[
            b"ACGTCAGTACGATCG",
            b"ACGTCAGTACGATCG",
            b"ACGTCAGGGGTACGATCG",
            b"ACGTCAGGGGTACGATCG",
        ];
        let msa = run(inputs);
        assert_eq!(
            msa.rows,
            vec![
                b"ACGTCA---GTACGATCG".to_vec(),
                b"ACGTCA---GTACGATCG".to_vec(),
                b"ACGTCAGGGGTACGATCG".to_vec(),
                b"ACGTCAGGGGTACGATCG".to_vec(),
            ]
        );
        assert_eq!(msa.ungapped_consensus(), b"ACGTCAGGGGTACGATCG");
        for (row, input) in msa.rows.iter().zip(inputs) {
            assert_eq!(row.len(), msa.consensus.len());
            assert_eq!(
                row.iter()
                    .copied()
                    .filter(|&c| c != b'-')
                    .collect::<Vec<_>>(),
                *input
            );
        }
    }

    #[test]
    fn consensus_ties_and_empty_inputs() {
        let tied = run(&[b"AAAAAA", b"CCCCCC"]);
        assert_eq!(
            tied.rows,
            vec![b"-----AAAAAA".to_vec(), b"CCCCCC-----".to_vec()]
        );
        assert_eq!(tied.consensus, b"CCCCCAAAAAA");
        assert_eq!(
            column_consensus(&[b"A".to_vec(), b"C".to_vec()], 0),
            (b'A', 1)
        );
        assert_eq!(run(&[b"-", b"A"]).consensus, b"A");
        assert_eq!(
            run(&[b"A", b"AA"]).rows,
            vec![b"-A".to_vec(), b"AA".to_vec()]
        );
        assert_eq!(align(&[]).unwrap(), Msa::default());
        assert!(align(&[Vec::new()]).is_err());
        assert!(align(&[b"AX".to_vec()]).is_err());
        assert_eq!(run(&[b"acgtN"]).rows[0], b"acgtN");
    }

    #[test]
    fn reference_tuple_windows_and_linkage() {
        assert_eq!(tuple_counts(b"AAAAA"), [0; 256]);
        assert_eq!(tuple_counts(b"AAAAAA")[0], 1);
        assert_eq!(tuple_counts(&vec![b'A'; 261])[0], 0);
        let inputs = vec![
            b"AAAAAA".to_vec(),
            b"CCCCCC".to_vec(),
            b"AAAAAA".to_vec(),
            b"CCCCCC".to_vec(),
        ];
        assert_eq!(guide_tree(&inputs), vec![(0, 2), (1, 3), (4, 5)]);
        assert_eq!(
            guide_tree(&vec![b"AAAAAA".to_vec(); 4]),
            vec![(0, 1), (4, 2), (5, 3)]
        );
    }

    #[test]
    fn conserved_block_and_smoothing() {
        let msa = run(&[b"ACGTACGT", b"ACGTTCGT", b"ACGTACGT"]);
        assert_eq!(
            msa.conserved_consensus(0.9).unwrap(),
            (b"ACGTACGT".to_vec(), Some((0, 7)))
        );
        assert!(msa.conserved_consensus(f64::NAN).is_err());
        assert_eq!(
            Msa::default().conserved_consensus(0.9).unwrap(),
            (vec![], None)
        );
    }
}

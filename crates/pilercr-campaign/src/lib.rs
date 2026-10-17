//! Deterministic synthetic FASTA generation and report classification.
//! The Rust PRNG stream intentionally differs from Python's random.Random stream.

pub const COMPARABLE: &str = "comparable";
pub const CRASH_EXEMPT: &str = "reference-crash-exempt";
pub const CRASH_MISMATCH: &str = "reference-crash-mismatch";

#[derive(Clone, Debug, PartialEq)]
pub struct CaseMeta {
    pub positive: bool,
    pub repeat_spans: Vec<(usize, usize)>,
    pub kind: String,
}

#[derive(Clone)]
pub struct Rng(u64);
impl Rng {
    pub fn new(seed: u64) -> Self {
        let state = seed ^ 0x9e3779b97f4a7c15;
        Self(if state == 0 {
            0xa0761d6478bd642f
        } else {
            state
        })
    }
    pub fn next_u64(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    pub fn range(&mut self, low: usize, high: usize) -> usize {
        let width = (high - low + 1) as u64;
        let threshold = width.wrapping_neg() % width;
        loop {
            let value = self.next_u64();
            if value >= threshold {
                return low + (value % width) as usize;
            }
        }
    }
    pub fn chance(&mut self, percent: u64) -> bool {
        self.next_u64() % 100 < percent
    }
    pub fn dna(&mut self, n: usize) -> String {
        (0..n).map(|_| b"ACGT"[self.range(0, 3)] as char).collect()
    }
}
pub fn derive_seed(master: u64, index: usize) -> u64 {
    let mut x = master ^ (index as u64).wrapping_mul(0x9e3779b97f4a7c15);
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94d049bb133111eb);
    x ^ (x >> 31)
}
fn array(r: &mut Rng, seq: &mut String, spans: &mut Vec<(usize, usize)>) {
    let len = r.range(28, 37);
    let consensus = r.dna(len);
    let n = r.range(3, 20);
    let spacer = r.range(30, 40);
    {
        let n = r.range(60, 300);
        seq.push_str(&r.dna(n));
    }
    for i in 0..n {
        let start = seq.len();
        let mut repeat = consensus.clone();
        if r.chance(30) {
            let p = r.range(0, len - 1);
            let mut b = r.dna(1);
            while b.as_bytes()[0] == repeat.as_bytes()[p] {
                b = r.dna(1);
            }
            repeat.replace_range(p..p + 1, &b);
        }
        seq.push_str(&repeat);
        spans.push((start, seq.len()));
        if i + 1 < n {
            let delta = r.range(0, 4) as isize - 2;
            let spacer_len = spacer.saturating_add_signed(delta).max(8);
            seq.push_str(&r.dna(spacer_len));
        }
    }
    {
        let n = r.range(60, 300);
        seq.push_str(&r.dna(n));
    }
}
pub fn generate_case(master: u64, index: usize, negative: bool) -> (u64, Vec<u8>, CaseMeta) {
    let seed = derive_seed(master, index);
    let mut r = Rng::new(seed);
    let mut seq = String::new();
    let mut spans = vec![];
    let mut kind = "positive".to_string();
    if negative {
        let k = r.range(0, 4);
        kind = ["background", "homopolymer", "tandem", "short", "ambiguous"][k].into();
        match k {
            0 => {
                let n = r.range(200, 2000);
                seq = r.dna(n);
            }
            1 => {
                seq = r.dna(80);
                seq.push_str(
                    &std::iter::repeat_n(b"ACGT"[r.range(0, 3)] as char, r.range(60, 300))
                        .collect::<String>(),
                );
                seq.push_str(&r.dna(80));
            }
            2 => {
                let n = r.range(2, 6);
                let unit = r.dna(n);
                seq = r.dna(60);
                for _ in 0..r.range(40, 120) {
                    seq.push_str(&unit)
                }
                seq.push_str(&r.dna(60));
            }
            3 => {
                let n = r.range(1, 40);
                seq = r.dna(n);
            }
            _ => {
                seq = (0..r.range(100, 400))
                    .map(|_| b"ACGTN"[r.range(0, 4)] as char)
                    .collect()
            }
        }
    } else {
        array(&mut r, &mut seq, &mut spans);
        if r.chance(30) {
            {
                let n = r.range(150, 400);
                seq.push_str(&r.dna(n));
            }
            array(&mut r, &mut seq, &mut spans);
        }
        if r.chance(20) {
            seq.make_ascii_lowercase();
        }
    }
    let label = format!("synthetic_{index:06}_{seed:016x}");
    (
        seed,
        format!(">{label}\n{seq}\n").into_bytes(),
        CaseMeta {
            positive: !negative,
            repeat_spans: spans,
            kind,
        },
    )
}
pub fn validated_array_count(report: &[u8], input: &str) -> Option<u64> {
    let lines: Vec<_> = report
        .split(|b| *b == b'\n')
        .filter(|l| !l.is_empty())
        .collect();
    let headers: Vec<_> = lines
        .iter()
        .filter_map(|line| {
            let s = std::str::from_utf8(line).ok()?.trim_end_matches('\r');
            let (path, tail) = s.rsplit_once(": ")?;
            let digits = tail.strip_suffix(" putative CRISPR arrays found.")?;
            if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
                return None;
            }
            let count = digits.parse().ok()?;
            Some((path, count))
        })
        .collect();
    if headers.len() != 1
        || std::path::Path::new(headers[0].0).file_name()?.to_str()?
            != std::path::Path::new(input).file_name()?.to_str()?
    {
        None
    } else {
        Some(headers[0].1)
    }
}
pub fn reference_crash(code: Option<i32>) -> bool {
    code.is_some_and(|c| c < 0)
}
pub fn classify(code: Option<i32>, rust_ok: bool, rust_report: &[u8], input: &str) -> &'static str {
    if !reference_crash(code) {
        COMPARABLE
    } else if rust_ok && validated_array_count(rust_report, input) == Some(0) {
        CRASH_EXEMPT
    } else {
        CRASH_MISMATCH
    }
}
pub fn classify_signal(
    signal: Option<i32>,
    rust_ok: bool,
    rust_report: &[u8],
    input: &str,
) -> &'static str {
    if signal.is_none() {
        COMPARABLE
    } else if rust_ok && validated_array_count(rust_report, input) == Some(0) {
        CRASH_EXEMPT
    } else {
        CRASH_MISMATCH
    }
}
pub fn mutation_step_substitutions(total: usize, steps: usize) -> Vec<usize> {
    if total == 0 || steps == 0 {
        return vec![];
    }
    (1..=steps).map(|s| (s * total).div_ceil(steps)).collect()
}
pub fn mutate_repeats(
    sequence: &str,
    spans: &[(usize, usize)],
    substitutions: usize,
    seed: u64,
) -> String {
    let mut positions: Vec<_> = spans.iter().flat_map(|&(a, b)| a..b).collect();
    let mut r = Rng::new(seed);
    for i in (1..positions.len()).rev() {
        let j = r.range(0, i);
        positions.swap(i, j)
    }
    let mut bytes = sequence.as_bytes().to_vec();
    for p in positions.into_iter().take(substitutions) {
        let upper = bytes[p].to_ascii_uppercase();
        let mut b = b"ACGT"[r.range(0, 3)];
        while b == upper {
            b = b"ACGT"[r.range(0, 3)]
        }
        bytes[p] = if bytes[p].is_ascii_lowercase() {
            b.to_ascii_lowercase()
        } else {
            b
        };
    }
    String::from_utf8(bytes).unwrap()
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn generation_is_repeatable_and_spans_bounded() {
        let a = generate_case(99, 1, false);
        assert_eq!(a, generate_case(99, 1, false));
        let n = a.1.iter().filter(|b| **b == b'\n').count();
        assert_eq!(n, 2);
        for &(s, e) in &a.2.repeat_spans {
            assert!(s < e && e < a.1.len())
        }
    }
    #[test]
    fn zero_xorshift_state_is_remapped() {
        let mut rng = Rng::new(0x9e3779b97f4a7c15);
        assert_ne!(rng.next_u64(), 0);
        let (_, fasta, _) = generate_case(17679207994573498142, 1, false);
        assert!(fasta.len() > 100);
    }
    #[test]
    fn count_header_requires_ascii_digits() {
        for count in ["+0", "-0", " 0", "0 "] {
            let report = format!("/tmp/a.fa: {count} putative CRISPR arrays found.\n");
            assert_eq!(validated_array_count(report.as_bytes(), "a.fa"), None);
            assert_eq!(
                classify_signal(Some(6), true, report.as_bytes(), "a.fa"),
                CRASH_MISMATCH
            );
        }
    }
    #[test]
    fn validation_and_crashes() {
        let r = b"/tmp/a.fa: 0 putative CRISPR arrays found.\n";
        assert_eq!(validated_array_count(r, "a.fa"), Some(0));
        assert_eq!(
            validated_array_count(&[r.as_slice(), r.as_slice()].concat(), "a.fa"),
            None
        );
        assert_eq!(classify(Some(-6), true, r, "a.fa"), CRASH_EXEMPT);
        assert_eq!(classify(Some(1), true, r, "a.fa"), COMPARABLE);
        assert_eq!(classify(Some(-6), false, r, "a.fa"), CRASH_MISMATCH)
    }
    #[test]
    fn mutation_caps() {
        let counts = mutation_step_substitutions(40, 4);
        assert_eq!(counts, vec![10, 20, 30, 40]);
        let s = "ACGTACGT";
        let m = mutate_repeats(s, &[(2, 6)], 2, 4);
        assert_eq!(s.len(), m.len());
        assert_eq!(&s[..2], &m[..2]);
        assert_eq!(&s[6..], &m[6..]);
    }
}

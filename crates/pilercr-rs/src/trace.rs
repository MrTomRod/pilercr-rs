use std::{
    cell::RefCell,
    fs::{self, File},
    io::Write,
    path::Path,
};

pub const STAGES: [&str; 16] = [
    "01-parameters",
    "02-contigs",
    "03-filter-hits",
    "04-trapezoids",
    "05-dp-hits",
    "06-images",
    "07-piles",
    "08-expanded-piles",
    "09-hit-piles",
    "10-edges",
    "11-components",
    "12-draft-arrays",
    "13-msa",
    "14-bounds",
    "15-decisions",
    "16-final-arrays",
];

struct Sink {
    files: Vec<File>,
    error: Option<String>,
    counters: std::collections::HashMap<String, usize>,
}

thread_local! {
    static SINK: RefCell<Option<Sink>> = const { RefCell::new(None) };
}

struct Scope(Option<Sink>);

impl Drop for Scope {
    fn drop(&mut self) {
        SINK.with(|sink| *sink.borrow_mut() = self.0.take());
    }
}

pub fn enabled() -> bool {
    SINK.with(|sink| sink.borrow().is_some())
}

pub fn next_id(name: &str) -> usize {
    SINK.with(|sink| {
        let mut sink = sink.borrow_mut();
        let Some(sink) = sink.as_mut() else { return 0 };
        let counter = sink.counters.entry(name.to_owned()).or_default();
        let id = *counter;
        *counter += 1;
        id
    })
}

pub fn emit(stage: &str, fields: &[String]) {
    SINK.with(|sink| {
        let mut sink = sink.borrow_mut();
        let Some(sink) = sink.as_mut() else { return };
        if sink.error.is_some() {
            return;
        }
        let Some(index) = STAGES.iter().position(|name| *name == stage) else {
            sink.error = Some(format!("unknown trace stage: {stage}"));
            return;
        };
        let mut row = String::new();
        for (index, field) in fields.iter().enumerate() {
            if index != 0 {
                row.push('\t');
            }
            for character in field.chars() {
                match character {
                    '\\' => row.push_str("\\\\"),
                    '\t' => row.push_str("\\t"),
                    '\r' => row.push_str("\\r"),
                    '\n' => row.push_str("\\n"),
                    other => row.push(other),
                }
            }
        }
        row.push('\n');
        if let Err(error) = sink.files[index].write_all(row.as_bytes()) {
            sink.error = Some(format!("trace {stage}: {error}"));
        }
    });
}

pub fn with_directory<T>(
    path: &Path,
    action: impl FnOnce() -> Result<T, String>,
) -> Result<T, String> {
    match fs::create_dir(path) {
        Ok(()) => {}
        Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {
            let mut entries = fs::read_dir(path).map_err(|e| format!("{}: {e}", path.display()))?;
            if entries
                .next()
                .transpose()
                .map_err(|e| format!("{}: {e}", path.display()))?
                .is_some()
            {
                return Err(format!("trace directory must be empty: {}", path.display()));
            }
        }
        Err(error) => return Err(format!("{}: {error}", path.display())),
    }
    let mut files = Vec::with_capacity(STAGES.len());
    for stage in STAGES {
        let filename = path.join(format!("{stage}.tsv"));
        let mut file = File::options()
            .write(true)
            .create_new(true)
            .open(&filename)
            .map_err(|e| format!("{}: {e}", filename.display()))?;
        writeln!(
            file,
            "# pilercr-trace-v1\t{stage}\tzero-based;native-endpoints"
        )
        .map_err(|e| format!("{}: {e}", filename.display()))?;
        files.push(file);
    }
    let _scope = Scope(SINK.with(|sink| {
        sink.replace(Some(Sink {
            files,
            error: None,
            counters: std::collections::HashMap::new(),
        }))
    }));
    let result = action();
    let error = SINK.with(|sink| sink.borrow_mut().as_mut().unwrap().error.take());
    match error {
        Some(error) => Err(error),
        None => result,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn directory() -> std::path::PathBuf {
        std::env::temp_dir().join(format!(
            "pilercr-trace-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ))
    }

    #[test]
    fn escapes_fields_and_restores_nested_scopes_on_error_and_unwind() {
        let outer = directory();
        let inner = directory();
        let panic_dir = directory();
        assert!(!enabled());
        emit("01-parameters", &[]);
        with_directory(&outer, || {
            assert_eq!(next_id("msa"), 0);
            assert_eq!(next_id("msa"), 1);
            let result: Result<(), String> = with_directory(&inner, || {
                assert_eq!(next_id("msa"), 0);
                emit("02-contigs", &["a\\b\tc\rd\ne".into()]);
                Err("action failed".into())
            });
            assert_eq!(result.unwrap_err(), "action failed");
            assert!(enabled());
            assert!(
                std::panic::catch_unwind(|| {
                    let _: Result<(), String> = with_directory(&panic_dir, || panic!("test"));
                })
                .is_err()
            );
            assert!(enabled());
            assert_eq!(next_id("msa"), 2);
            emit("01-parameters", &["outer".into()]);
            Ok(())
        })
        .unwrap();
        assert!(!enabled());
        assert_eq!(
            fs::read_to_string(inner.join("02-contigs.tsv")).unwrap(),
            "# pilercr-trace-v1\t02-contigs\tzero-based;native-endpoints\na\\\\b\\tc\\rd\\ne\n"
        );
        assert!(
            fs::read_to_string(outer.join("01-parameters.tsv"))
                .unwrap()
                .ends_with("outer\n")
        );
        assert!(with_directory(&outer, || Ok(())).is_err());
        for path in [outer, inner, panic_dir] {
            fs::remove_dir_all(path).unwrap();
        }
    }

    #[test]
    fn creation_and_write_errors_propagate() {
        let path = directory();
        fs::write(&path, "keep").unwrap();
        assert!(with_directory(&path, || Ok(())).is_err());
        assert_eq!(fs::read_to_string(&path).unwrap(), "keep");
        fs::remove_file(&path).unwrap();
        assert!(with_directory(&path.join("missing").join("trace"), || Ok(())).is_err());
        #[cfg(unix)]
        {
            with_directory(&path, || {
                SINK.with(|sink| {
                    sink.borrow_mut().as_mut().unwrap().files[0] =
                        File::open(path.join("01-parameters.tsv")).unwrap()
                });
                emit("01-parameters", &["not writable".into()]);
                emit("01-parameters", &["ignored after first error".into()]);
                Ok(())
            })
            .unwrap_err();
            assert!(!enabled());
            fs::remove_dir_all(path).unwrap();
        }
    }
}

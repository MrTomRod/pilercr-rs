import argparse
import math
from pathlib import Path
import sys
import tempfile
import unittest


STAGES = (
    "01-parameters", "02-contigs", "03-filter-hits", "04-trapezoids",
    "05-dp-hits", "06-images", "07-piles", "08-expanded-piles",
    "09-hit-piles", "10-edges", "11-components", "12-draft-arrays",
    "13-msa", "14-bounds", "15-decisions", "16-final-arrays",
)
FLOAT_COLUMNS = {"01-parameters": 6, "05-dp-hits": 8}
ESCAPES = {"\\": "\\", "t": "\t", "n": "\n", "r": "\r"}


class InvalidTrace(ValueError):
    pass


def header(stage):
    return f"# pilercr-trace-v1\t{stage}\tzero-based;native-endpoints\n"


def unescape(value):
    result = []
    index = 0
    while index < len(value):
        char = value[index]
        if char == "\\":
            index += 1
            if index == len(value) or value[index] not in ESCAPES:
                raise InvalidTrace(f"malformed escape in {value!r}")
            char = ESCAPES[value[index]]
        result.append(char)
        index += 1
    return "".join(result)


def record(stream, path, stage, number):
    line = stream.readline()
    if not line:
        return None
    location = f"{path}: record {number} (line {number + 1})"
    if not line.endswith("\n"):
        raise InvalidTrace(f"{location}: truncated record (missing final LF)")
    if "\r" in line:
        raise InvalidTrace(f"{location}: unescaped CR")
    fields = line[:-1].split("\t")
    try:
        for field in fields:
            unescape(field)
        column = FLOAT_COLUMNS.get(stage)
        if column is not None:
            if len(fields) <= column:
                raise InvalidTrace(f"missing float column {column + 1} (index {column})")
            if not math.isfinite(float(fields[column])):
                raise InvalidTrace(f"non-finite float in column {column + 1}")
    except ValueError as error:
        raise InvalidTrace(f"{location}: {error}") from error
    return fields


def compare(left, right, abs_tol=0.0, rel_tol=0.0):
    try:
        if any(not math.isfinite(tol) or tol < 0 for tol in (abs_tol, rel_tol)):
            raise InvalidTrace("tolerances must be finite and nonnegative")
        paths = [(stage, Path(left) / f"{stage}.tsv", Path(right) / f"{stage}.tsv")
                 for stage in STAGES]
        for stage, left_path, right_path in paths:
            for path in (left_path, right_path):
                if not path.is_file():
                    raise InvalidTrace(f"missing stage file: {path}")
        for stage, left_path, right_path in paths:
            for path in (left_path, right_path):
                with path.open(encoding="utf-8", newline="") as stream:
                    actual = stream.readline()
                if actual != header(stage):
                    raise InvalidTrace(f"{path}: invalid header: {actual!r}; expected {header(stage)!r}")
        for stage, left_path, right_path in paths:
            with left_path.open(encoding="utf-8", newline="") as lhs, right_path.open(encoding="utf-8", newline="") as rhs:
                lhs.readline()
                rhs.readline()
                number = 0
                while True:
                    number += 1
                    a = record(lhs, left_path, stage, number)
                    b = record(rhs, right_path, stage, number)
                    if a is None and b is None:
                        break
                    location = f"{stage}: record {number}"
                    if a is None or b is None:
                        side = "LEFT" if a is None else "RIGHT"
                        present = b if a is None else a
                        return 1, f"{location}: missing row in {side}; present row: {[unescape(v) for v in present]!r}"
                    for column in range(max(len(a), len(b))):
                        x = a[column] if column < len(a) else None
                        y = b[column] if column < len(b) else None
                        where = f"{location}, column {column + 1} (index {column})"
                        if x is None or y is None:
                            left_value = "<missing>" if x is None else repr(unescape(x))
                            right_value = "<missing>" if y is None else repr(unescape(y))
                            return 1, f"{where}: field difference: LEFT={left_value}, RIGHT={right_value}"
                        if column == FLOAT_COLUMNS.get(stage):
                            if not math.isclose(float(x), float(y), abs_tol=abs_tol, rel_tol=rel_tol):
                                return 1, (f"{where}: float difference: LEFT={x!r}, RIGHT={y!r}, "
                                           f"absolute delta={abs(float(x) - float(y))!r}, "
                                           f"abs_tol={abs_tol}, rel_tol={rel_tol}")
                        elif x != y:
                            return 1, f"{where}: field difference: LEFT={unescape(x)!r}, RIGHT={unescape(y)!r}"
        return 0, "All 16 trace stages match."
    except (OSError, UnicodeError, InvalidTrace) as error:
        return 2, f"Invalid trace: {error}"


class ComparatorTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.left = Path(self.temp.name) / "left"
        self.right = Path(self.temp.name) / "right"
        self.left.mkdir()
        self.right.mkdir()
        for stage in STAGES:
            for directory in (self.left, self.right):
                self.put(directory, stage, "")

    def put(self, directory, stage, rows):
        (directory / f"{stage}.tsv").write_text(header(stage) + rows, encoding="utf-8")

    def result(self, expected, **kwargs):
        code, message = compare(self.left, self.right, **kwargs)
        self.assertEqual(code, expected, message)
        return message

    def test_exact_and_escapes(self):
        for directory in (self.left, self.right):
            self.put(directory, "02-contigs", "id\ta\\tb\\nc\\rd\\\\e\n")
        self.result(0)
        self.assertEqual(unescape(r"a\tb\nc\rd\\e"), "a\tb\nc\rd\\e")

    def test_float_tolerance(self):
        for stage, column in FLOAT_COLUMNS.items():
            with self.subTest(stage=stage):
                prefix = "x\t" * column
                self.put(self.left, stage, prefix + "10\n")
                self.put(self.right, stage, prefix + "10.01\n")
                self.assertIn("float difference", self.result(1))
                self.result(0, abs_tol=0.02)
                self.result(0, rel_tol=0.002)
                self.put(self.right, stage, prefix + "10.0\n")
                self.result(0)

    def test_field_order_and_display(self):
        self.put(self.left, "02-contigs", "a\\tb\nsecond\n")
        self.put(self.right, "02-contigs", "second\na\\tb\n")
        message = self.result(1)
        self.assertIn("02-contigs: record 1, column 1", message)
        self.assertIn(repr("a\tb"), message)

    def test_header(self):
        path = self.right / "16-final-arrays.tsv"
        path.write_text(header("15-decisions"), encoding="utf-8")
        self.assertIn("invalid header", self.result(2))

    def test_missing_file_before_earlier_difference(self):
        self.put(self.left, "02-contigs", "different\n")
        (self.right / "16-final-arrays.tsv").unlink()
        self.assertIn("missing stage file", self.result(2))

    def test_missing_row(self):
        self.put(self.left, "02-contigs", "one\ntwo\n")
        self.put(self.right, "02-contigs", "one\n")
        self.assertIn("record 2: missing row in RIGHT", self.result(1))

    def test_missing_column(self):
        self.put(self.left, "02-contigs", "one\ttwo\n")
        self.put(self.right, "02-contigs", "one\n")
        self.assertIn("column 2 (index 1)", self.result(1))

    def test_truncation(self):
        self.put(self.right, "02-contigs", "unfinished")
        self.assertIn("truncated record", self.result(2))

    def test_malformed_escapes(self):
        for value in ("bad\\q\n", "bad\\\n"):
            with self.subTest(value=value):
                self.put(self.right, "02-contigs", value)
                self.assertIn("malformed escape", self.result(2))

    def test_invalid_float(self):
        for value in ("nan", "inf", "not-a-number"):
            with self.subTest(value=value):
                self.put(self.right, "05-dp-hits", "x\t" * 8 + value + "\n")
                self.result(2)

    def test_invalid_tolerance(self):
        for value in (-1, math.inf, math.nan):
            with self.subTest(value=value):
                self.result(2, abs_tol=value)
                self.result(2, rel_tol=value)


def tolerance(value):
    try:
        number = float(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("expected a finite nonnegative number") from error
    if not math.isfinite(number) or number < 0:
        raise argparse.ArgumentTypeError("expected a finite nonnegative number")
    return number


def main(argv=None):
    parser = argparse.ArgumentParser(description="Compare ordered, escaped TSV traces across all 16 stages. Exit: 0 match, 1 difference, 2 invalid input.")
    parser.add_argument("leftdir", nargs="?")
    parser.add_argument("rightdir", nargs="?")
    parser.add_argument("--abs-tol", type=tolerance, default=0.0)
    parser.add_argument("--rel-tol", type=tolerance, default=0.0)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args(argv)
    if args.self_test:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(ComparatorTests)
        return 0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1
    if args.leftdir is None or args.rightdir is None:
        parser.error("LEFTDIR and RIGHTDIR are required unless --self-test is used")
    code, message = compare(args.leftdir, args.rightdir, args.abs_tol, args.rel_tol)
    print(message, file=sys.stderr if code == 2 else sys.stdout)
    return code


if __name__ == "__main__":
    sys.exit(main())

"""korean_gate.py -- Hangul scan over source files.

A Python character-range loop, never `grep -P`: the PCRE class can fail by
locale, and a `|| echo PASS` fallback then prints a false PASS. Here a file
that cannot be read is an ERROR, not a pass, and a run that scanned zero
files is a failure -- "no match" and "failed to run" must not look alike.
"""
import sys, os

RANGES = [(0x1100, 0x11FF, "Hangul Jamo"), (0x3130, 0x318F, "Compatibility Jamo"),
          (0xA960, 0xA97F, "Jamo Extended-A"), (0xAC00, 0xD7A3, "Hangul Syllables"),
          (0xD7B0, 0xD7FF, "Jamo Extended-B"), (0xFFA0, 0xFFDC, "Halfwidth Hangul")]

def block(ch):
    cp = ord(ch)
    for lo, hi, name in RANGES:
        if lo <= cp <= hi:
            return name
    return None

def main(roots, exts):
    files, chars, hits, errors = 0, 0, [], []
    for root in roots:
        for dirpath, dirnames, filenames in os.walk(root):
            dirnames[:] = [d for d in dirnames if d not in (".git", "node_modules")]
            for fn in sorted(filenames):
                if not any(fn.endswith(e) for e in exts):
                    continue
                path = os.path.join(dirpath, fn)
                try:
                    text = open(path, encoding="utf-8").read()
                except Exception as exc:
                    errors.append((path, repr(exc)))
                    continue
                files += 1
                chars += len(text)
                for lineno, line in enumerate(text.split("\n"), 1):
                    for col, ch in enumerate(line, 1):
                        b = block(ch)
                        if b:
                            hits.append((path, lineno, col, ch, b))
    print(f"scanned {files} file(s), {chars} character(s), "
          f"extensions {' '.join(exts)}")
    for path, exc in errors:
        print(f"  ERROR unreadable: {path}: {exc}")
    for path, lineno, col, ch, b in hits:
        print(f"  HANGUL {path}:{lineno}:{col} U+{ord(ch):04X} ({b})")
    if files == 0:
        print("GATE FAILED: zero files scanned - a gate that inspected nothing "
              "has not passed")
        return 1
    if errors or hits:
        print(f"GATE FAILED: {len(hits)} Hangul character(s), "
              f"{len(errors)} unreadable file(s)")
        return 1
    print("GATE PASSED: no Hangul in the scanned sources")
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[2:] or ["."], sys.argv[1].split(",")))

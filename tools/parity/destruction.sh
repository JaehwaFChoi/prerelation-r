#!/bin/bash
# destruction.sh -- deliberately break each harness input and confirm the
# check fails. For every check the question asked is: what would this input
# look like if the defect were present? If the answer is "the same", the
# check is decoration (Finding recIm4AOaEX20IyUq).
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
RLIB="${R_PRERELATION_LIB:-$HOME/rlib}"
GOLD="${PRERELATION_GOLDEN:?set PRERELATION_GOLDEN to the reference package tests/golden directory}"
PASS=0; FAIL=0
WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT
cd "$(dirname "$0")"

# The harness outputs these two files. Earlier this script assumed they were
# already on disk from a previous run; when they were not, D1/D3/D6 failed at
# the *read* and the non-zero exit was scored as the gate firing. A check that
# errors is not a check that passed (Finding recVJAACoUcnnfL9C), so the inputs
# are produced here and their existence is asserted.
R_PRERELATION_LIB="$RLIB" Rscript parity_r.R "$GOLD" "$WORK"/parity_r_out.csv >/dev/null 2>&1
R_PRERELATION_LIB="$RLIB" Rscript condense_r.R graphs.json > "$WORK"/condense_r.txt 2>/dev/null
for f in "$WORK"/parity_r_out.csv "$WORK"/condense_r.txt; do
  [ -s "$f" ] || { echo "SETUP FAILED: $f absent or empty - the destruction tests below would report false passes"; exit 2; }
done
echo "setup: "$WORK"/parity_r_out.csv $(wc -l < "$WORK"/parity_r_out.csv) lines, "$WORK"/condense_r.txt $(wc -l < "$WORK"/condense_r.txt) lines"
echo
verdict () { if [ "$1" = "$2" ]; then echo "  [OK]   $3"; PASS=$((PASS+1)); else echo "  [BAD]  $3 (expected exit $2, got $1)"; FAIL=$((FAIL+1)); fi }

echo "D1 perturb one R value (product/PI + 1e-9): the tolerance check must fire"
sed 's|^product,PI,0.61690185543994114$|product,PI,0.6169018564399411|' "$WORK"/parity_r_out.csv > /tmp/d1.csv
diff -q "$WORK"/parity_r_out.csv /tmp/d1.csv >/dev/null && echo "  !! mutation did not change the file - test is void"
python parity_driver.py $GOLD /tmp/d1.csv >/tmp/d1.log 2>&1; verdict $? 1 "driver rejects a perturbed R value"
grep -E "^product +PI" /tmp/d1.log | head -1

echo
echo "D2 foreign permutation pool, on the SCAN layer: the exact p-value check must fire"
# This test used to run on the six golden fixtures, and it could not fire
# there for a reason worth recording. Their permutation p-values are
# saturated (count 0 of 199, or 199 of 199), so substituting a different
# valid permutation SET leaves the count -- and therefore the p-value --
# unchanged. The break lands on the input and the statistic has no dynamic
# range to register it (Finding recbTS1KoQldEFfrs). The scan layer is where
# the evidence lives: its twelve ordered pairs carry interior p-values.
SCAN="$WORK/scan"; mkdir -p "$SCAN"
python scan_driver.py --setup "$SCAN" >/dev/null 2>&1
R_PRERELATION_LIB="$RLIB" Rscript scan_r.R "$SCAN" >/dev/null 2>&1
cp "$SCAN"/perm_scan.csv "$WORK"/perm_scan_orig.csv
python - "$WORK"/perm_scan_orig.csv "$SCAN"/perm_scan.csv <<'PY2'
import sys, random
src, dst = sys.argv[1], sys.argv[2]
rng = random.Random(20260828)
rows = [l.strip() for l in open(src) if l.strip()]
out = []
for r in rows:
    v = [int(t) for t in r.split(",")]
    rng.shuffle(v)
    out.append(",".join(str(t) for t in v))
open(dst, "w").write("\n".join(out) + "\n")
PY2
# Prove the break landed, before reading any verdict.
cmp -s "$WORK"/perm_scan_orig.csv "$SCAN"/perm_scan.csv && echo "  !! pool unchanged - test is void"
python - "$WORK"/perm_scan_orig.csv "$SCAN"/perm_scan.csv <<'PY2'
import sys
a = sorted(open(sys.argv[1]).read().split())
b = sorted(open(sys.argv[2]).read().split())
print("  row multiset differs:", a != b, "(False would mean the test is void)")
assert a != b, "mutation left the permutation SET unchanged - test is void"
PY2
# R recomputes from the mutated pool; Python compares against its own,
# computed from the pool it read before the mutation.
R_PRERELATION_LIB="$RLIB" Rscript scan_r.R "$SCAN" >/dev/null 2>&1
cp "$WORK"/perm_scan_orig.csv "$SCAN"/perm_scan.csv
python scan_driver.py --compare "$SCAN" >"$WORK"/d2.log 2>&1; verdict $? 1 "scan driver rejects a foreign permutation pool"
grep -c "FAIL" "$WORK"/d2.log | sed 's/^/  FAIL lines: /'
grep -E "p_value.*FAIL" "$WORK"/d2.log | head -2

echo "D3 drop a quantity from the R output: coverage must fire, not pass silently"
grep -v ",A1," "$WORK"/parity_r_out.csv > /tmp/d3.csv
python parity_driver.py $GOLD /tmp/d3.csv >/tmp/d3.log 2>&1; verdict $? 1 "driver rejects a missing quantity"
grep "missing in R output" /tmp/d3.log | head -2

echo
echo "D4 perturb the fixture data: the package golden test must fire"
rm -rf /tmp/gold_bad2 && cp -r $GOLD /tmp/gold_bad2
python - <<'PY'
p = "/tmp/gold_bad2/fixture_product.csv"
lines = open(p).read().split("\n")
lines[1] = "0.9,0.05"          # one person moved to a corner
open(p, "w").write("\n".join(lines))
print("  fixture row 1 replaced with 0.9,0.05")
PY
R_PRERELATION_LIB="$RLIB" Rscript "$REPO"/tests/run_tests.R /tmp/gold_bad2 >/tmp/d4.log 2>&1; verdict $? 1 "package golden test rejects perturbed fixture data"
grep -E "^  FAIL  product" /tmp/d4.log | head -3
tail -2 /tmp/d4.log

echo
echo "D5 no golden directory: the golden check must report SKIP, not PASS"
R_PRERELATION_LIB="$RLIB" Rscript "$REPO"/tests/run_tests.R >/tmp/d5.log 2>&1; verdict $? 0 "run without fixtures still exits 0"
grep -E "SKIP|skipped" /tmp/d5.log | head -3

echo
echo "D6 change one graph edge for the JS side only: the condense diff must fire"
sed 's|\["a1","a3"\]|["a3","a1"]|' graphs.json > /tmp/graphs_bad.json
cmp -s graphs.json /tmp/graphs_bad.json && echo "  !! graph unchanged - test is void"
node condense_js.mjs /tmp/graphs_bad.json > /tmp/d6_js.txt 2>&1
diff -q "$WORK"/condense_r.txt /tmp/d6_js.txt >/dev/null; verdict $? 1 "condense diff detects a changed graph"
diff "$WORK"/condense_r.txt /tmp/d6_js.txt | head -4

echo
echo "destruction summary: $PASS behaved as required, $FAIL did not"
[ $FAIL -eq 0 ] || exit 1

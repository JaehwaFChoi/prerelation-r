// condense_js.mjs -- the reference condense from prerelation-js, printed in
// the same canonical form. JS class ids are zero-based; +1 is applied here
// so the two renderings are directly comparable (the same off-by-one the
// shared permutation matrices carry).
import { readFileSync } from "node:fs";
const JS_SRC = process.env.PRERELATION_JS_SRC || "../../../prerelation-js/src/scan.mjs";
const { condense } = await import(JS_SRC);
const graphs = JSON.parse(readFileSync(process.argv[2], "utf8"));
for (const [name, g] of Object.entries(graphs)) {
  const r = condense(g.nodes, g.edges);
  console.log(`graph=${name}`);
  console.log(`  classes=${r.classes.map((c) => c.join("+")).join("|")}`);
  console.log(`  class_of=${g.nodes.map((u) => `${u}:${r.classOf.get(u) + 1}`).join(",")}`);
  const qe = r.quotientEdges.length === 0 ? "-" : r.quotientEdges.map(([a, b]) => `${a + 1}>${b + 1}`).join(",");
  const he = r.hasseEdges.length === 0 ? "-" : r.hasseEdges.map(([a, b]) => `${a + 1}>${b + 1}`).join(",");
  console.log(`  quotient=${qe}`);
  console.log(`  hasse=${he}`);
}

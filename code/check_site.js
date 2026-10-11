// Load the documentation site in a simulated browser and report what it builds.
// Any script error, failed fetch or empty view is printed; nothing is silently swallowed.
const fs = require("fs");
const path = require("path");
const { JSDOM, VirtualConsole } = require("jsdom");

const DOCS = path.resolve(__dirname, "..", "docs");
const errors = [];
const vc = new VirtualConsole();
vc.on("jsdomError", e => errors.push("jsdomError: " + (e.stack || e.message)));
vc.on("error", (...a) => errors.push("console.error: " + a.join(" ")));
vc.on("warn", (...a) => errors.push("console.warn: " + a.join(" ")));

const dom = new JSDOM(fs.readFileSync(path.join(DOCS, "index.html"), "utf8"), {
  runScripts: "outside-only",
  pretendToBeVisual: true,
  url: "https://s-matysik.github.io/glassbox-ftir/",
  virtualConsole: vc,
});
const w = dom.window;

// fetch served straight off disk, exactly what GitHub Pages returns
w.fetch = (url) => {
  const rel = String(url).replace(/^https?:\/\/[^/]+\/glassbox-ftir\//, "");
  const p = path.join(DOCS, rel);
  if (!fs.existsSync(p)) return Promise.resolve({ ok: false, status: 404 });
  const body = fs.readFileSync(p, "utf8");
  return Promise.resolve({ ok: true, status: 200, json: () => Promise.resolve(JSON.parse(body)),
                           text: () => Promise.resolve(body) });
};
w.matchMedia = w.matchMedia || (() => ({ matches: false, addListener() {}, addEventListener() {} }));

for (const f of ["assets/charts.js", "assets/app.js"]) {
  try {
    w.eval(fs.readFileSync(path.join(DOCS, f), "utf8"));
  } catch (e) {
    errors.push(`EXCEPTION while evaluating ${f}: ${e.stack || e.message}`);
  }
}

const VIEWS = ["overview", "spectra", "decision", "parity", "rules", "overlap",
               "learning", "stability", "acceptance", "data"];

function render(view) {
  return new Promise(resolve => {
    const before = errors.length;
    w.location.hash = "#" + view;
    w.dispatchEvent(new w.Event("hashchange"));
    setTimeout(() => {
      const d = w.document;
      const views = d.querySelector("#views");
      resolve({
        view,
        svg: views ? views.querySelectorAll("svg").length : 0,
        paths: views ? views.querySelectorAll("svg path, svg rect, svg circle, svg line").length : 0,
        controls: views ? views.querySelectorAll("select, input, button").length : 0,
        chars: views ? views.textContent.trim().length : 0,
        newErrors: errors.length - before,
      });
    }, 900);
  });
}

function census(view) {
  const svg = w.document.querySelector("#views svg");
  if (!svg) return "no svg";
  const c = {};
  svg.querySelectorAll("*").forEach(n => { c[n.tagName] = (c[n.tagName] || 0) + 1; });
  return Object.entries(c).sort((a, b) => b[1] - a[1])
    .map(([k, v]) => `${k}:${v}`).join(" ");
}

async function interact(view) {
  const sel = w.document.querySelector("#views select, #views input[type=range]");
  if (!sel) return "no control";
  const before = w.document.querySelector("#views svg").innerHTML.length;
  if (sel.tagName === "SELECT" && sel.options.length > 1) {
    sel.selectedIndex = sel.options.length - 1;
  } else {
    sel.value = sel.max || sel.value;
  }
  sel.dispatchEvent(new w.Event("change", { bubbles: true }));
  sel.dispatchEvent(new w.Event("input", { bubbles: true }));
  await new Promise(r => setTimeout(r, 700));
  const svg = w.document.querySelector("#views svg");
  const after = svg ? svg.innerHTML.length : 0;
  return after === 0 ? "chart disappeared" : (after === before ? "no change" : "redrew");
}

setTimeout(async () => {
  const d = w.document;
  console.log("nav entries: " + d.querySelectorAll("nav a, nav button").length);
  const bad = [];
  for (const v of VIEWS) {
    const r = await render(v);
    if (v !== "overview") {
      console.log("     " + v.padEnd(11) + " elements: " + census(v));
      console.log("     " + v.padEnd(11) + " control:  " + await interact(v));
    }
    const ok = r.view === "overview" ? r.chars > 300 : (r.svg > 0 && r.paths > 3);
    if (!ok) bad.push(r.view);
    console.log(`${ok ? "ok  " : "FAIL"} ${r.view.padEnd(11)} svg=${r.svg} marks=${String(r.paths).padEnd(5)} ` +
                `controls=${String(r.controls).padEnd(3)} chars=${String(r.chars).padEnd(6)} newErrors=${r.newErrors}`);
  }
  console.log("views failing: " + (bad.length ? bad.join(", ") : "none"));
  console.log("total errors: " + errors.length);
  errors.slice(0, 8).forEach(e => console.log("  - " + e.split("\n").slice(0, 3).join("\n    ")));
  process.exit(bad.length || errors.length ? 1 : 0);
}, 3000);

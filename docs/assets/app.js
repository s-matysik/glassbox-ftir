/* glassbox-ftir documentation site: view switching, cross-filtering and
   download of exactly what is on screen. Every number comes from docs/data/,
   exported by code/20_export_web.py from the tables in results/. */
(function () {
  "use strict";

  var POLY = ["HDPE", "LDPE", "PET", "PP", "PS", "PVC"];
  var ATTRS = ["A_W1", "C_W1", "O_W1", "A_W2", "C_W2", "O_W2", "A_W3", "C_W3", "O_W3"];
  var KIND = { A: "amplitude", C: "shape correlation", O: "baseline offset" };
  var WIN = { W1: "C-H stretch", W2: "carbonyl", W3: "CH deformation" };
  var PCOL = { HDPE: "#1f4e79", LDPE: "#4a90c2", PET: "#2e7d32", PP: "#e0a458", PS: "#8e5bab", PVC: "#c0392b" };
  var D = {}, state = { view: "overview", poly: POLY.slice(), k: 4 };

  function attrLabel(a) { var p = a.split("_"); return KIND[p[0]] + ", " + WIN[p[1]]; }
  function col(f) { return GB.css(f); }
  function toObjs(t) {
    return t.rows.map(function (r) {
      var o = {}; t.columns.forEach(function (c, i) { o[c] = r[i]; }); return o;
    });
  }
  function q(sel, root) { return (root || document).querySelector(sel); }

  function panel(opts) {
    var d = document.createElement("section");
    d.className = "panel";
    d.innerHTML =
      '<header><div><h2>' + opts.title + '</h2>' +
      (opts.sub ? '<div class="cap">' + opts.sub + '</div>' : '') + '</div>' +
      '<div class="actions">' +
      '<button data-act="png">PNG</button><button data-act="csv">CSV</button></div></header>' +
      '<div class="controls"></div><div class="plot"></div>' +
      (opts.cap ? '<p class="cap">' + opts.cap + '</p>' : '');
    q("#views").appendChild(d);
    d.querySelector('[data-act="png"]').onclick = function () { GB.savePNG(d.querySelector(".plot"), opts.file + ".png"); };
    d.querySelector('[data-act="csv"]').onclick = function () {
      var t = opts.csv();
      GB.download(opts.file + ".csv", GB.toCSV(t.columns, t.rows));
    };
    return { root: d, plot: d.querySelector(".plot"), controls: d.querySelector(".controls") };
  }
  function control(host, label, node) {
    var w = document.createElement("label");
    w.className = "control";
    w.innerHTML = "<span>" + label + "</span>";
    w.appendChild(node);
    host.appendChild(w);
    return node;
  }
  function select(options, value, onchange) {
    var s = document.createElement("select");
    options.forEach(function (o) {
      var op = document.createElement("option");
      op.value = o.v !== undefined ? o.v : o; op.textContent = o.k !== undefined ? o.k : o;
      s.appendChild(op);
    });
    s.value = value; s.onchange = function () { onchange(s.value); };
    return s;
  }
  function slider(min, max, step, value, oninput) {
    var wrap = document.createElement("span");
    var s = document.createElement("input");
    s.type = "range"; s.min = min; s.max = max; s.step = step; s.value = value;
    var out = document.createElement("b");
    out.style.fontFamily = "var(--mono)"; out.style.marginLeft = "7px";
    out.textContent = value;
    s.oninput = function () { out.textContent = s.value; oninput(+s.value); };
    wrap.appendChild(s); wrap.appendChild(out);
    return wrap;
  }

  /* ---------------------------------------------------------------- views */

  var VIEWS = {};

  VIEWS.overview = function () {
    var par = toObjs(D.parity), st = toObjs(D.stability.summary);
    var cov = par.filter(function (r) { return r.method === "glassbox_covered"; })[0];
    var svm = par.filter(function (r) { return r.method === "svm_rbf"; })[0];
    var gset = toObjs(D.rough.by_set)[0];
    var sh = st.filter(function (r) { return /post-hoc/.test(r.representation); })[0];
    var sm = toObjs(D.selective.stability_matched);
    var gbsel = sm.filter(function (r) { return /selected attributes/.test(r.selector); })[0];
    var shap2 = sm.filter(function (r) { return /attribution/.test(r.selector) && r.set_size === 2; })[0];
    var cards = [
      { v: GB.nice(cov.bacc_mean, 4), k: "balanced accuracy on accepted specimens", c: "ok" },
      { v: (cov.coverage_mean * 100).toFixed(1) + "%", k: "specimens the rule accepts" },
      { v: GB.nice(svm.bacc_mean, 4), k: "strongest black box, all specimens" },
      { v: GB.nice(gset.gamma, 4), k: "cells unambiguous in training, best pair", c: "ok" },
      { v: Math.round(gbsel.modal_share * gbsel.repeats) + " / " + gbsel.repeats,
        k: "repeats with an unchanged attribute pair", c: "ok" },
      { v: Math.round(shap2.modal_share * shap2.repeats) + " / " + shap2.repeats,
        k: "repeats with an unchanged attribution top pair" }
    ];
    var g = document.createElement("div");
    g.className = "stat-grid";
    g.innerHTML = cards.map(function (c) {
      return '<div class="stat ' + (c.c || "") + '"><div class="v">' + c.v + '</div><div class="k">' + c.k + '</div></div>';
    }).join("");
    q("#views").appendChild(g);

    var d = document.createElement("section");
    d.className = "panel";
    d.innerHTML =
      "<h2>What this site is</h2>" +
      "<p class='cap'>Supporting material for the article <b>A glass-box classifier for FTIR identification of plastics</b>. " +
      "Every chart reads a table produced by the analysis scripts, so nothing here can disagree with the manuscript. " +
      "Each panel offers the figure as a print-resolution PNG and the data behind it as CSV.</p>" +
      "<div class='note'>The headline result is not that the interpretable model is more accurate. It is not. " +
      "It reaches parity <b>on the specimens it accepts</b>, and the fraction it refers for a second test is stated " +
      "rather than hidden. That referral rate is the price of auditability, and this site lets you inspect what it buys.</div>" +
      "<p class='cap'>Use the polymer chips at the top to cross-filter every view. Views keep their settings in the address bar, " +
      "so a specific chart can be linked directly from the paper.</p>";
    q("#views").appendChild(d);
  };

  VIEWS.spectra = function () {
    var S = D.spectra;
    var mode = "mean";
    var p = panel({
      title: "Spectrum explorer",
      sub: "Mean spectrum per polymer with the three analysis windows shaded.",
      file: "spectra",
      cap: "The three windows are fixed before any data is seen: <b>W1</b> 2800 to 3000 cm-1 carbon-hydrogen stretching, " +
        "<b>W2</b> 1650 to 1800 cm-1 carbonyl, <b>W3</b> 1350 to 1500 cm-1 methylene and methyl deformation. " +
        "Nine attributes are the amplitude, the shape correlation and the baseline offset of each window. " +
        "Corresponds to Figure 1 of the article.",
      csv: function () {
        var cols = ["wavenumber_cm-1"].concat(state.poly.map(function (x) { return x + "_mean"; }));
        return {
          columns: cols,
          rows: S.grid.map(function (g, i) {
            return [g].concat(state.poly.map(function (x) { return S.polymers[x].mean[i]; }));
          })
        };
      }
    });
    control(p.controls, "show", select(
      [{ v: "mean", k: "mean with 5 to 95 percent band" }, { v: "samples", k: "individual specimens" }],
      mode, function (v) { mode = v; draw(); }));
    function draw() {
      var ymin = Infinity, ymax = -Infinity;
      state.poly.forEach(function (name) {
        var o = S.polymers[name];
        (mode === "mean" ? [o.p05, o.p95] : o.samples).forEach(function (arr) {
          arr.forEach(function (v) { if (v < ymin) ymin = v; if (v > ymax) ymax = v; });
        });
      });
      if (!isFinite(ymin)) { p.plot.innerHTML = "<p class='cap'>No polymer selected.</p>"; return; }
      var pad = (ymax - ymin) * 0.06;
      var f = new GB.Frame(p.plot, { height: 360 });
      f.scales([S.grid[0], S.grid[S.grid.length - 1]], [ymin - pad, ymax + pad]);
      Object.keys(S.windows).forEach(function (w) {
        f.vspan(S.windows[w][0], S.windows[w][1], col("--accent"), w);
      });
      f.axes("wavenumber (cm-1)", "transmittance", { yfmt: function (v) { return v.toFixed(0); } });
      state.poly.forEach(function (name) {
        var o = S.polymers[name], c = PCOL[name];
        if (mode === "mean") {
          f.band(S.grid.map(function (g, i) { return [g, o.p05[i], o.p95[i]]; }), c);
          f.line(S.grid.map(function (g, i) { return [g, o.mean[i]]; }), c, { w: 1.7 });
        } else {
          o.samples.forEach(function (s) {
            f.line(S.grid.map(function (g, i) { return [g, s[i]]; }), c, { w: 0.8 });
          });
        }
      });
      GB.legend(p.plot, state.poly.map(function (x) { return { k: x, c: PCOL[x] }; }), function () { });
    }
    draw();
  };

  VIEWS.decision = function () {
    var sets = toObjs(D.rough.by_set), byp = toObjs(D.rough.by_polymer);
    var chosen = sets[0].attribute_set;
    var p = panel({
      title: "Where the decision is certain",
      sub: "Rough-set approximation: how much of the sample a named attribute set decides without ambiguity.",
      file: "decision",
      cap: "A specimen is <b>decided</b> when its cell of the attribute grid contains only one polymer, and <b>referred</b> otherwise. " +
        "Cell width is the attribute's own detection threshold, so two specimens share a cell exactly when the instrument cannot tell them apart. " +
        "Corresponds to Table 3 of the article.",
      csv: function () { return D.rough.by_polymer; }
    });
    control(p.controls, "attribute set", select(
      sets.slice(0, 12).map(function (s) {
        return { v: s.attribute_set, k: s.attribute_set.split("+").map(attrLabel).join("  +  ") + "   (gamma " + GB.nice(s.gamma, 3) + ")" };
      }), chosen, function (v) { chosen = v; draw(); }));
    function draw() {
      p.plot.innerHTML = "";
      var set = sets.filter(function (s) { return s.attribute_set === chosen; })[0];
      var rows = byp.filter(function (r) { return r.attribute_set === chosen && state.poly.indexOf(r.polymer) >= 0; });
      var head = document.createElement("div");
      head.className = "stat-grid";
      head.innerHTML =
        '<div class="stat ok"><div class="v">' + GB.nice(set.gamma, 4) + '</div><div class="k">quality of approximation, all six polymers</div></div>' +
        '<div class="stat"><div class="v">' + set.n_cells + '</div><div class="k">cells in the attribute grid</div></div>' +
        '<div class="stat"><div class="v">' + set.n_attributes + '</div><div class="k">attributes used</div></div>';
      p.plot.appendChild(head);
      var host = document.createElement("div");
      p.plot.appendChild(host);
      GB.barsH(host, rows.map(function (r) {
        return {
          k: r.polymer, v: r.alpha, c: PCOL[r.polymer],
          label: GB.nice(r.alpha, 3) + "   (" + r.boundary_size + " referred)",
          tip: "<b>" + r.polymer + "</b><br>accuracy of approximation " + GB.nice(r.alpha, 4) +
            "<br>boundary region " + r.boundary_size + " specimens"
        };
      }), { max: 1, xlab: "accuracy of approximation (1 = no boundary region)", left: 86 });
    }
    draw();
  };

  VIEWS.parity = function () {
    var rows = toObjs(D.parity);
    var nice = {
      glassbox_forced: "glass box, forced choice", glassbox_covered: "glass box, on accepted",
      svm_rbf: "SVM", random_forest: "random forest", plsda: "PLS-DA", cnn1d: "1-D CNN"
    };
    var p = panel({
      title: "Accuracy parity on accepted specimens",
      sub: "Balanced accuracy under five-fold cross-validation, with the interval across folds.",
      file: "parity",
      cap: "The glass box appears twice. Forced to decide every specimen it is <b>0.052</b> behind the strongest baseline, " +
        "a difference that is statistically significant. Allowed to refer, it reaches parity on the <b>75.3 percent</b> it accepts. " +
        "Corresponds to Figure 3 and Table 2 of the article.",
      csv: function () { return D.parity; }
    });
    rows.sort(function (a, b) { return a.bacc_mean - b.bacc_mean; });
    GB.barsH(p.plot, rows.map(function (r) {
      return {
        k: nice[r.method] || r.method, v: r.bacc_mean,
        c: /glass/.test(r.method) ? col("--warn") : col("--muted"),
        label: GB.nice(r.bacc_mean, 4) + (r.coverage_mean < 0.999 ? "   (" + (r.coverage_mean * 100).toFixed(1) + "% decided)" : ""),
        tip: "<b>" + (nice[r.method] || r.method) + "</b><br>balanced accuracy " + GB.nice(r.bacc_mean, 4) +
          "<br>95% interval " + GB.nice(r.ci_low, 4) + " to " + GB.nice(Math.min(1, r.ci_high), 4) +
          "<br>coverage " + GB.nice(r.coverage_mean, 3)
      };
    }), { min: 0.9, max: 1.0, xlab: "balanced accuracy", left: 150, height: 250 });
  };

  VIEWS.rules = function () {
    var rows = toObjs(D.rules);
    var pairs = [];
    rows.forEach(function (r) { if (pairs.indexOf(r.body) < 0) pairs.push(r.body); });
    var p = panel({
      title: "Which attribute responds to which polymer change",
      sub: "Confidence that an attribute moves beyond its detection threshold when the polymer changes.",
      file: "rules",
      cap: "Each cell is a rule mined over pairs of specimens that differ in polymer and agree on everything else. " +
        "Hover for the Wilson interval, the pair count and the false-discovery-rate q value. " +
        "Corresponds to Figure 4 and Table 4 of the article.",
      csv: function () { return D.rules; }
    });
    function draw() {
      p.plot.innerHTML = "";
      var shown = pairs.filter(function (b) {
        var two = b.split(" -> ");
        if (two.length < 2) return true;   // the pooled "any polymer change" rule
        return state.poly.indexOf(two[0]) >= 0 && state.poly.indexOf(two[1]) >= 0;
      });
      if (shown.length < 2) { p.plot.innerHTML = "<p class='cap'>Select at least two polymers.</p>"; return; }
      var idx = {};
      rows.forEach(function (r) { idx[r.attribute + "|" + r.body] = r; });
      var order = ATTRS.slice().sort(function (a, b) {
        var m = function (x) {
          var s = 0, n = 0;
          shown.forEach(function (bb) { var r = idx[x + "|" + bb]; if (r) { s += r.conf_change; n++; } });
          return n ? s / n : 0;
        };
        return m(b) - m(a);
      });
      GB.heatmap(p.plot, {
        rows: order, cols: shown, min: 0.2, max: 1,
        rowLabel: attrLabel,
        colLabel: function (c) { return c.replace(" -> ", " / ").replace("any polymer change", "pooled"); },
        value: function (a, b) { var r = idx[a + "|" + b]; return r ? r.conf_change : null; },
        tip: function (a, b, v) {
          var r = idx[a + "|" + b];
          if (!r) return "no rule";
          return "<b>" + attrLabel(a) + "</b><br>" + b + "<br>confidence " + GB.nice(v, 3) +
            "<br>Wilson " + GB.nice(r.wilson_low, 3) + " to " + GB.nice(r.wilson_high, 3) +
            "<br>pairs " + r.n_pairs + "<br>q = " + GB.nice(r.q_fdr, 4) +
            (r.significant_q05 ? "<br>significant at q &lt; 0.05" : "<br>not significant");
        }
      });
    }
    draw();
    p.root._refresh = draw;
  };

  VIEWS.overlap = function () {
    var rows = toObjs(D.ovl);
    var p = panel({
      title: "How well each attribute separates polymers",
      sub: "Overlap coefficient between the distributions of two polymers, by attribute.",
      file: "overlap",
      cap: "Zero means the two distributions do not overlap at all; one means they are indistinguishable. " +
        "Corresponds to Figure 5 of the article.",
      csv: function () { return D.ovl; }
    });
    var mode = "median";
    control(p.controls, "show", select(
      [{ v: "median", k: "median over polymer pairs" }, { v: "matrix", k: "every attribute by pair" }],
      mode, function (v) { mode = v; draw(); }));
    function draw() {
      p.plot.innerHTML = "";
      var shown = rows.filter(function (r) {
        var two = r.comparison.split(" vs ");
        return state.poly.indexOf(two[0]) >= 0 && state.poly.indexOf(two[1]) >= 0;
      });
      if (!shown.length) { p.plot.innerHTML = "<p class='cap'>Select at least two polymers.</p>"; return; }
      if (mode === "median") {
        var med = ATTRS.map(function (a) {
          var v = shown.filter(function (r) { return r.attribute === a; }).map(function (r) { return r.ovl; }).sort(function (x, y) { return x - y; });
          return { k: attrLabel(a), v: v.length ? v[Math.floor(v.length / 2)] : 0, raw: a };
        }).sort(function (x, y) { return x.v - y.v; });
        GB.barsH(p.plot, med.map(function (m) {
          return {
            k: m.k, v: m.v, c: m.v < 0.3 ? col("--accent") : col("--muted"),
            tip: "<b>" + m.k + "</b><br>median overlap " + GB.nice(m.v, 4)
          };
        }), { max: 1, xlab: "median overlap coefficient (lower is better)", left: 186 });
      } else {
        var cmp = [];
        shown.forEach(function (r) { if (cmp.indexOf(r.comparison) < 0) cmp.push(r.comparison); });
        var idx = {};
        shown.forEach(function (r) { idx[r.attribute + "|" + r.comparison] = r.ovl; });
        GB.heatmap(p.plot, {
          rows: ATTRS, cols: cmp, min: 0, max: 1, rowLabel: attrLabel,
          colLabel: function (c) { return c.replace(" vs ", " / "); },
          value: function (a, c) { var v = idx[a + "|" + c]; return v === undefined ? null : v; },
          tip: function (a, c, v) { return "<b>" + attrLabel(a) + "</b><br>" + c + "<br>overlap " + GB.nice(v, 4); }
        });
      }
    }
    draw();
    p.root._refresh = draw;
  };

  VIEWS.learning = function () {
    var rows = toObjs(D.learning);
    var methods = [];
    rows.forEach(function (r) { if (methods.indexOf(r.method) < 0) methods.push(r.method); });
    var hidden = {};
    var MC = {
      glassbox_forced: col("--warn"), svm_rbf: "#1f4e79", plsda: "#4a90c2",
      random_forest: "#6b6b6b", cnn1d: "#2e7d32"
    };
    var p = panel({
      title: "Data efficiency: the hypothesis is not supported",
      sub: "Balanced accuracy against training set size, three repetitions of five-fold cross-validation at each size.",
      file: "learning",
      cap: "Interpretable models are often claimed to need less data. Here the glass box is behind at every size and " +
        "there is no crossover: the gap widens slightly from 1 to 2 percent and then narrows steadily. " +
        "Corresponds to Figure 6 and Table 6 of the article.",
      csv: function () { return D.learning; }
    });
    function draw() {
      p.plot.innerHTML = "";
      var vis = methods.filter(function (m) { return !hidden[m]; });
      var ys = [], xs = [];
      rows.forEach(function (r) {
        if (vis.indexOf(r.method) < 0) return;
        ys.push(r.bacc_ci_low, r.bacc_ci_high); xs.push(r.size_pct);
      });
      if (!ys.length) { p.plot.innerHTML = "<p class='cap'>No method selected.</p>"; return; }
      var f = new GB.Frame(p.plot, { height: 340, logx: true });
      f.scales([Math.min.apply(null, xs), Math.max.apply(null, xs)],
        [Math.min.apply(null, ys) - 0.02, Math.min(1.005, Math.max.apply(null, ys) + 0.015)]);
      f.axes("training specimens (% of the full set)", "balanced accuracy",
        { yfmt: function (v) { return v.toFixed(2); }, xfmt: function (v) { return v + "%"; } });
      vis.forEach(function (m) {
        var g = rows.filter(function (r) { return r.method === m; }).sort(function (a, b) { return a.size_pct - b.size_pct; });
        var c = MC[m] || col("--muted");
        f.band(g.map(function (r) { return [r.size_pct, r.bacc_ci_low, r.bacc_ci_high]; }), c);
        f.line(g.map(function (r) { return [r.size_pct, r.balanced_accuracy]; }), c, {
          w: m === "glassbox_forced" ? 2.3 : 1.5, markers: true,
          tip: function (pt) {
            var r = g.filter(function (x) { return x.size_pct === pt[0]; })[0];
            return "<b>" + m.replace(/_/g, " ") + "</b><br>" + r.size_pct + "% of specimens, " +
              r.n_per_class + " per class<br>balanced accuracy " + GB.nice(r.balanced_accuracy, 4) +
              "<br>95% interval " + GB.nice(r.bacc_ci_low, 4) + " to " + GB.nice(r.bacc_ci_high, 4);
          }
        });
      });
      GB.legend(p.plot, methods.map(function (m) { return { k: m.replace(/_/g, " "), c: MC[m] || col("--muted") }; }),
        function (k, on) {
          var m = methods.filter(function (x) { return x.replace(/_/g, " ") === k; })[0];
          hidden[m] = !on; draw();
        });
    }
    draw();
  };

  VIEWS.stability = function () {
    var ranks = toObjs(D.stability.ranks), sum = toObjs(D.stability.summary);
    var p = panel({
      title: "Does the explanation stay the same when the model is refitted",
      sub: "Twenty fits on independently drawn training samples.",
      file: "stability",
      cap: "Left: the rank each attribute receives from post-hoc attribution, fit by fit. Right: the attribute set the rule selects. " +
        "The attribution ranking is reasonably stable overall, but its top three changes in one repeat in ten, " +
        "while the rule selects the same pair in 19 of 20 fits. Corresponds to Figure 7 and Table 7 of the article.",
      csv: function () { return D.stability.ranks; }
    });
    var sh = sum.filter(function (r) { return /post-hoc/.test(r.representation); })[0];
    var gbRow = sum.filter(function (r) { return /selected attributes/.test(r.representation); })[0];
    var gbSets = gbRow ? gbRow.n_distinct_sets : 2;
    var cards = document.createElement("div");
    cards.className = "stat-grid";
    cards.innerHTML =
      '<div class="stat warn"><div class="v">' + sh.n_distinct_sets + '</div><div class="k">distinct attribution top-three sets</div></div>' +
      '<div class="stat warn"><div class="v">' + (sh.modal_set_share * 100).toFixed(0) + '%</div><div class="k">repeats sharing the modal set</div></div>' +
      '<div class="stat"><div class="v">' + GB.nice(sh.kendall_tau_mean, 3) + '</div><div class="k">mean rank correlation between fits</div></div>' +
      '<div class="stat ok"><div class="v">' + gbSets + '</div><div class="k">distinct glass-box attribute sets</div></div>';
    p.plot.appendChild(cards);
    var host = document.createElement("div");
    p.plot.appendChild(host);
    var f = new GB.Frame(host, { height: 320, margin: { t: 14, r: 170, b: 42, l: 54 } });
    f.scales([1, ranks.length], [ATTRS.length + 0.4, 0.6]);
    f.axes("repeated fit", "rank (1 = most important)", { ny: ATTRS.length, yfmt: function (v) { return v.toFixed(0); } });
    var pal = ["#1f4e79", "#4a90c2", "#2e7d32", "#e0a458", "#8e5bab", "#c0392b", "#6b6b6b", "#17a2b8", "#b5651d"];
    ATTRS.forEach(function (a, i) {
      f.line(ranks.map(function (r) { return [r.repeat, r[a]]; }), pal[i], {
        w: 1.5, markers: true,
        tip: function (pt) { return "<b>" + attrLabel(a) + "</b><br>fit " + pt[0] + ", rank " + pt[1]; }
      });
    });
    var lg = GB.el("g", { transform: "translate(" + (f.iw + 14) + ",6)" }, f.g);
    ATTRS.forEach(function (a, i) {
      GB.el("rect", { x: 0, y: i * 17, width: 11, height: 3, fill: pal[i], rx: 1.5 }, lg);
      GB.el("text", { x: 17, y: i * 17 + 4, fill: f.muted, "font-size": 10.5 }, lg).textContent = attrLabel(a);
    });
    var note = document.createElement("div");
    note.className = "note";
    var modal = {};
    ranks.forEach(function (r) { modal[r.glassbox_set] = (modal[r.glassbox_set] || 0) + 1; });
    var best = Object.keys(modal).sort(function (a, b) { return modal[b] - modal[a]; })[0];
    note.innerHTML = "The rule selected <b>" + best.split(", ").map(attrLabel).join("</b> and <b>") +
      "</b> in " + modal[best] + " of " + ranks.length + " fits.";
    p.plot.appendChild(note);
  };

  VIEWS.referral = function () {
    var matched = toObjs(D.selective.matched);
    var pairs = toObjs(D.selective.pairs);
    var per = toObjs(D.selective.per_class);
    var nice = {
      glassbox_forced: "glass box", glassbox: "glass box",
      svm_rbf: "SVM", random_forest: "random forest", plsda: "PLS-DA", cnn1d: "1-D CNN"
    };
    var own = matched.filter(function (r) { return /own reject/.test(r.regime); });
    var p = panel({
      title: "What the referral buys, and what it costs",
      sub: "Every classifier held to the same referral budget of 75.3 percent coverage.",
      file: "referral",
      cap: "Given a reject option of their own, the margin-based baselines are perfect on the specimens they keep, " +
        "while the rule reaches <b>0.9973</b>. Selective accuracy is therefore not the advantage: the advantage is " +
        "that the refusal itself is stated in named spectral quantities. The measured cost at equal coverage is " +
        "<b>0.0027</b> balanced accuracy. Corresponds to Table 9 of the article.",
      csv: function () { return D.selective.matched; }
    });
    own.sort(function (a, b) { return a.balanced_accuracy - b.balanced_accuracy; });
    GB.barsH(p.plot, own.map(function (r) {
      return {
        k: nice[r.method] || r.method, v: r.balanced_accuracy,
        c: /glass/.test(r.method) ? col("--warn") : col("--muted"),
        label: GB.nice(r.balanced_accuracy, 4),
        tip: "<b>" + (nice[r.method] || r.method) + "</b><br>balanced accuracy " +
          GB.nice(r.balanced_accuracy, 4) + "<br>MCC " + GB.nice(r.mcc, 4) +
          "<br>coverage " + GB.nice(r.coverage, 3)
      };
    }), { min: 0.99, max: 1.0, xlab: "balanced accuracy on accepted specimens", left: 150, height: 190 });

    var h2 = document.createElement("h3");
    h2.textContent = "Who gets referred";
    h2.style.margin = "22px 0 6px";
    p.plot.appendChild(h2);
    var hostA = document.createElement("div");
    p.plot.appendChild(hostA);
    per.sort(function (a, b) { return a.coverage - b.coverage; });
    GB.barsH(hostA, per.map(function (r) {
      return {
        k: r.polymer, v: r.coverage, c: PCOL[r.polymer] || col("--accent"),
        label: (r.coverage * 100).toFixed(1) + "%",
        tip: "<b>" + r.polymer + "</b><br>" + r.n_accepted + " of 500 accepted" +
          "<br>recall on accepted " + GB.nice(r.recall_accepted, 4)
      };
    }), { min: 0, max: 1, xlab: "share of specimens the rule accepts", left: 90, height: 210 });

    var h3 = document.createElement("h3");
    h3.textContent = "Which attribute pair, judged out of fold";
    h3.style.margin = "22px 0 6px";
    p.plot.appendChild(h3);
    var hostB = document.createElement("div");
    p.plot.appendChild(hostB);
    var f = new GB.Frame(hostB, { height: 290, margin: { t: 14, r: 20, b: 44, l: 58 } });
    var xs = pairs.map(function (r) { return r.coverage_out_of_fold; });
    var ys = pairs.map(function (r) { return r.bacc_covered; });
    f.scales([Math.min.apply(null, xs) - 0.03, Math.max.apply(null, xs) + 0.03],
             [Math.min.apply(null, ys) - 0.01, Math.max.apply(null, ys) + 0.005]);
    f.axes("coverage, out of fold", "accuracy on accepted specimens");
    var tipFor = function (pt) {
      var r = pt[2];
      return "<b>" + r.attributes.replace("+", " and ") + "</b><br>coverage " +
        GB.nice(r.coverage_out_of_fold, 4) + "<br>accepted " + GB.nice(r.bacc_covered, 4) +
        "<br>forced " + GB.nice(r.bacc_forced, 4) +
        "<br>training quality of approximation " + GB.nice(r.gamma_train_mean, 4);
    };
    f.dots(pairs.map(function (r) {
      return [r.coverage_out_of_fold, r.bacc_covered, r,
              r.attributes === "A_W2+C_W3" ? col("--warn") : col("--accent")];
    }), col("--accent"), { r: 5.5, alpha: 0.95, tip: tipFor });
    pairs.forEach(function (r) {
      GB.el("text", {
        x: f.x(r.coverage_out_of_fold) + 8, y: f.y(r.bacc_covered) + 3.5,
        fill: f.muted, "font-size": 10.5
      }, f.g).textContent = r.attributes.replace("+", "/");
    });
    var note = document.createElement("div");
    note.className = "note";
    note.innerHTML = "The highlighted point is the pair the classifier selects by inner cross-validation. " +
      "It is first on coverage and on accuracy at once. The pair with the highest training-set quality of " +
      "approximation, <b>O_W1 and O_W2</b>, sits below and to the left: that criterion overfits the training partition.";
    p.plot.appendChild(note);
  };

  VIEWS.acceptance = function () {
    var rows = toObjs(D.oc);
    var plans = [];
    rows.forEach(function (r) { var k = r.n + "|" + r.c; if (plans.indexOf(k) < 0) plans.push(k); });
    var plan = plans[0];
    var p = panel({
      title: "What the referral rate is worth in acceptance sampling",
      sub: "Operating-characteristic curves under a perfect inspector and under the measured one.",
      file: "acceptance",
      cap: "An imperfect test shifts the whole curve. At a one percent defective fraction the supplier of conforming material " +
        "carries almost the entire cost, while the buyer's protection barely changes. A referral outcome routes the ambiguous " +
        "specimen to a second test instead of rejecting the lot. Corresponds to Figure 10 and Table 8 of the article.",
      csv: function () { return D.oc; }
    });
    control(p.controls, "sampling plan", select(
      plans.map(function (k) { var s = k.split("|"); return { v: k, k: "n = " + s[0] + ", accept on c = " + s[1] }; }),
      plan, function (v) { plan = v; draw(); }));
    function draw() {
      p.plot.innerHTML = "";
      var s = plan.split("|");
      var sel = rows.filter(function (r) { return String(r.n) === s[0] && String(r.c) === s[1]; });
      var f = new GB.Frame(p.plot, { height: 340 });
      f.scales([0, Math.max.apply(null, sel.map(function (r) { return r.lot_fraction_defective; })) * 100], [0, 1.02]);
      f.axes("lot fraction defective (%)", "probability of acceptance", { yfmt: function (v) { return v.toFixed(1); } });
      [["perfect inspector", col("--muted"), "4 3"], ["observed inspector", col("--warn"), null]].forEach(function (cfg) {
        var g = sel.filter(function (r) { return r.inspector === cfg[0]; })
          .sort(function (a, b) { return a.lot_fraction_defective - b.lot_fraction_defective; });
        f.line(g.map(function (r) { return [r.lot_fraction_defective * 100, r.prob_accept]; }), cfg[1],
          { w: 2, dash: cfg[2] });
      });
      [1, 10].forEach(function (v) {
        GB.el("line", { x1: f.x(v), x2: f.x(v), y1: 0, y2: f.ih, stroke: f.muted, "stroke-width": 1, "stroke-dasharray": "2 3" }, f.g);
      });
      var tbl = document.createElement("table");
      var risk = [1, 10].map(function (pct) {
        var row = ["at " + pct + "% defective"];
        ["perfect inspector", "observed inspector"].forEach(function (ins) {
          var g = sel.filter(function (r) { return r.inspector === ins; });
          var near = g.reduce(function (a, b) {
            return Math.abs(b.lot_fraction_defective * 100 - pct) < Math.abs(a.lot_fraction_defective * 100 - pct) ? b : a;
          });
          row.push(GB.nice(near.prob_accept, 4));
        });
        return row;
      });
      tbl.innerHTML = "<tr><th>Acceptance probability</th><th class='num'>perfect inspector</th><th class='num'>observed inspector</th></tr>" +
        risk.map(function (r) { return "<tr><td>" + r[0] + "</td><td class='num'>" + r[1] + "</td><td class='num'>" + r[2] + "</td></tr>"; }).join("");
      p.plot.appendChild(tbl);
      GB.legend(p.plot, [{ k: "perfect inspector", c: col("--muted") }, { k: "observed inspector", c: col("--warn") }], function () { });
    }
    draw();
  };

  VIEWS.data = function () {
    var ent = toObjs(D.entry), sens = toObjs(D.sensitivity);
    var p = panel({
      title: "Entry condition",
      sub: "Between-specimen spread divided by the instrument noise standard deviation.",
      file: "entry",
      cap: "Below about one the method cannot work at all, because real differences are smaller than what the instrument " +
        "produces by itself. Corresponds to Figure 2 of the article.",
      csv: function () { return D.entry; }
    });
    GB.barsH(p.plot, ent.filter(function (r) { return state.poly.indexOf(r.stratum) >= 0; }).map(function (r) {
      return {
        k: r.stratum, v: r.spread_over_noise, c: PCOL[r.stratum], label: GB.nice(r.spread_over_noise, 1),
        tip: "<b>" + r.stratum + "</b><br>spread " + GB.nice(r.sd_between, 4) +
          "<br>noise " + GB.nice(r.noise_sd, 4) +
          "<br>spread / noise " + GB.nice(r.spread_over_noise, 1) +
          "<br>spread / threshold " + GB.nice(r.spread_over_threshold, 1)
      };
    }), { xlab: "between-specimen spread / instrument noise", left: 86 });

    var p2 = panel({
      title: "Sensitivity to the protocol's open choices",
      sub: "Headline rule confidence under every setting the frozen protocol left variable.",
      file: "sensitivity",
      cap: "Varying the pair sampling cap, the threshold multiplier k and the spectral resolution moves the result by " +
        "less than two percentage points. Corresponds to Table 5 of the article.",
      csv: function () { return D.sensitivity; }
    });
    var conf = sens.filter(function (r) { return r.metric === "conf_change" && r.attribute === "A_W1"; });
    var t = document.createElement("table");
    t.innerHTML = "<tr><th>Varied</th><th>Setting</th><th class='num'>Confidence</th></tr>" +
      conf.map(function (r) {
        return "<tr><td>" + r.variant_type + "</td><td>" + r.variant_value + "</td><td class='num'>" + GB.nice(r.value, 4) + "</td></tr>";
      }).join("");
    p2.plot.appendChild(t);

    var d = document.createElement("section");
    d.className = "panel";
    d.innerHTML =
      "<h2>Data, licences and what they can support</h2>" +
      "<table><tr><th>Role</th><th>Collection</th><th>Licence</th></tr>" +
      "<tr><td>Primary</td><td>FTIR-Plastics, Zenodo record 10736650, 3000 specimens of six industrial polymers</td><td>CC BY 4.0</td></tr>" +
      "<tr><td>Transfer case</td><td>Recyclate FTIR archive, Zenodo 21822279, 356 spectra over 54 particles</td><td>CC BY 4.0</td></tr></table>" +
      "<div class='note warn'><b>Read before reusing the transfer archive.</b> Its polymer, morphology and colour fields record the " +
      "reference-library entry that was matched to each spectrum, not an independently established property of the particle, and its " +
      "pass or fail verdict comes from the archive's own automated thresholds rather than from inspection. Only 17 of its 54 particles " +
      "carry a polymer label, which is why the methodological claims in this work rest on the primary collection instead.</div>" +
      "<h2 style='margin-top:18px'>Limitations</h2>" +
      "<p class='cap'>The evaluation completed one repetition of five-fold cross-validation rather than the twenty the protocol " +
      "specifies, so the intervals shown here are computed across five folds and are narrower than repeated cross-validation would give. " +
      "The archive does not document whether the 500 specimens per polymer are distinct physical items. The detection threshold is " +
      "estimated from filter residuals rather than from replicate acquisitions, so it is a lower bound and the referral rate with it. " +
      "Mined rules are symmetric by construction and pairs are not independent, so a pair count is not a sample size. " +
      "Acceptance-sampling results are a design calculation, not a lot-level validation. The full list is in the article and in " +
      "<code>PROTOCOL_v1.1.md</code>.</p>";
    q("#views").appendChild(d);
  };

  /* -------------------------------------------------------------- routing */

  var NAV = [
    ["overview", "Overview"], ["spectra", "Spectrum explorer"], ["decision", "Where the decision is certain"],
    ["parity", "Accuracy parity"], ["rules", "Rule matrix"], ["overlap", "Separation"],
    ["learning", "Data efficiency"], ["stability", "Explanation stability"],
    ["referral", "What referral costs"], ["acceptance", "Acceptance sampling"], ["data", "Data and limits"]
  ];

  function render() {
    q("#views").innerHTML = "";
    (VIEWS[state.view] || VIEWS.overview)();
    document.querySelectorAll("nav a").forEach(function (a) {
      a.classList.toggle("active", a.dataset.view === state.view);
    });
    q("#chips").style.display = ["overview", "parity", "learning", "stability", "referral"].indexOf(state.view) >= 0 ? "none" : "flex";
    location.replace("#" + state.view);
  }

  function boot() {
    var nav = q("nav");
    NAV.forEach(function (v, i) {
      var a = document.createElement("a");
      a.href = "#" + v[0]; a.dataset.view = v[0];
      a.innerHTML = '<span class="n">' + String(i).padStart(2, "0") + '</span>' + v[1];
      a.onclick = function (e) { e.preventDefault(); state.view = v[0]; render(); };
      nav.appendChild(a);
    });
    var chips = q("#chips");
    POLY.forEach(function (x) {
      var c = document.createElement("span");
      c.className = "chip on"; c.textContent = x;
      c.onclick = function () {
        c.classList.toggle("on");
        state.poly = POLY.filter(function (y) { return q('.chip[data-p="' + y + '"]').classList.contains("on"); });
        render();
      };
      c.dataset.p = x;
      chips.appendChild(c);
    });
    q("#theme").onclick = function () {
      var dark = document.documentElement.getAttribute("data-theme") === "dark";
      document.documentElement.setAttribute("data-theme", dark ? "light" : "dark");
      try { localStorage.setItem("gb-theme", dark ? "light" : "dark"); } catch (e) { }
      render();
    };
    try {
      var saved = localStorage.getItem("gb-theme");
      if (saved) document.documentElement.setAttribute("data-theme", saved);
    } catch (e) { }
    var h = location.hash.replace("#", "");
    if (VIEWS[h]) state.view = h;
    addEventListener("hashchange", function () {
      var v = location.hash.replace("#", "");
      if (VIEWS[v] && v !== state.view) { state.view = v; render(); }
    });
    addEventListener("resize", function () { clearTimeout(boot._t); boot._t = setTimeout(render, 180); });
    render();
  }

  var FILES = ["spectra", "parity", "rough", "rules", "ovl", "learning", "oc", "stability", "entry", "sensitivity", "selective"];
  Promise.all(FILES.map(function (f) {
    return fetch("data/" + f + ".json").then(function (r) {
      if (!r.ok) throw new Error(f + ": HTTP " + r.status);
      return r.json();
    });
  })).then(function (res) {
    FILES.forEach(function (f, i) { D[f] = res[i]; });
    boot();
  }).catch(function (err) {
    q("#views").innerHTML = "<div class='panel'><h2>Could not load the data files</h2>" +
      "<p class='cap'>" + err.message + ". When opening the site from a local copy, serve the folder over HTTP " +
      "(<code>python -m http.server</code> inside <code>docs/</code>) rather than opening the file directly, " +
      "because browsers block local fetch requests.</p></div>";
  });
})();

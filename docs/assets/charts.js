/* Minimal SVG chart layer for the glassbox-ftir documentation site.
   No external dependency: the site must render from a downloaded copy of the
   Zenodo record with no network access. Supports line, band, scatter, bar,
   horizontal bar and heatmap, with hover read-out, legend toggling and
   PNG / CSV download of exactly what is on screen. */
(function (global) {
  "use strict";

  var NS = "http://www.w3.org/2000/svg";
  var tip = null;

  function el(name, attrs, parent) {
    var n = document.createElementNS(NS, name);
    for (var k in attrs) if (attrs[k] !== undefined && attrs[k] !== null) n.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(n);
    return n;
  }
  function css(name) {
    return getComputedStyle(document.body).getPropertyValue(name).trim();
  }
  function ensureTip() {
    if (!tip) { tip = document.createElement("div"); tip.className = "tip"; document.body.appendChild(tip); }
    return tip;
  }
  function showTip(e, html) {
    var t = ensureTip();
    t.innerHTML = html; t.classList.add("on");
    var x = e.clientX + 14, y = e.clientY + 14;
    if (x + t.offsetWidth > innerWidth - 8) x = e.clientX - t.offsetWidth - 12;
    if (y + t.offsetHeight > innerHeight - 8) y = e.clientY - t.offsetHeight - 12;
    t.style.left = x + "px"; t.style.top = y + "px";
  }
  function hideTip() { if (tip) tip.classList.remove("on"); }

  function nice(v, d) {
    if (v === null || v === undefined || isNaN(v)) return "n/a";
    return Number(v).toFixed(d === undefined ? 3 : d);
  }
  function ticks(lo, hi, n) {
    if (lo === hi) { lo -= 0.5; hi += 0.5; }
    var span = hi - lo, step = Math.pow(10, Math.floor(Math.log10(span / n)));
    var err = span / (n * step);
    if (err >= 7.5) step *= 10; else if (err >= 3.5) step *= 5; else if (err >= 1.5) step *= 2;
    var out = [], t = Math.ceil(lo / step) * step;
    for (; t <= hi + step * 1e-9; t += step) out.push(+t.toFixed(10));
    return out;
  }
  function logTicks(lo, hi) {
    var out = [];
    for (var e = Math.floor(Math.log10(lo)); e <= Math.ceil(Math.log10(hi)); e++)
      for (var m = 1; m <= 9; m++) {
        var v = m * Math.pow(10, e);
        if (v >= lo * 0.999 && v <= hi * 1.001 && (m === 1 || m === 2 || m === 5)) out.push(v);
      }
    return out;
  }

  function Frame(host, opts) {
    opts = opts || {};
    host.innerHTML = "";
    var W = host.clientWidth || 760, H = opts.height || 330;
    var m = Object.assign({ t: 16, r: 18, b: 42, l: 62 }, opts.margin || {});
    var svg = el("svg", { class: "chart", viewBox: "0 0 " + W + " " + H, width: "100%", height: H }, host);
    this.svg = svg; this.W = W; this.H = H; this.m = m;
    this.iw = W - m.l - m.r; this.ih = H - m.t - m.b;
    this.g = el("g", { transform: "translate(" + m.l + "," + m.t + ")" }, svg);
    this.host = host;
    // Named rule rather than `line`: an own property of that name would shadow
    // Frame.prototype.line, the method that draws a series.
    this.ink = css("--ink"); this.muted = css("--muted"); this.rule = css("--line");
    this.logx = !!opts.logx;
  }
  Frame.prototype.scales = function (xd, yd) {
    var self = this;
    this.xd = xd; this.yd = yd;
    this.x = this.logx
      ? function (v) { return (Math.log10(v) - Math.log10(xd[0])) / (Math.log10(xd[1]) - Math.log10(xd[0])) * self.iw; }
      : function (v) { return (v - xd[0]) / (xd[1] - xd[0]) * self.iw; };
    this.y = function (v) { return self.ih - (v - yd[0]) / (yd[1] - yd[0]) * self.ih; };
    return this;
  };
  Frame.prototype.axes = function (xlab, ylab, o) {
    o = o || {};
    var self = this, xs = this.logx ? logTicks(this.xd[0], this.xd[1]) : ticks(this.xd[0], this.xd[1], o.nx || 6);
    var ys = ticks(this.yd[0], this.yd[1], o.ny || 5);
    ys.forEach(function (v) {
      el("line", { x1: 0, x2: self.iw, y1: self.y(v), y2: self.y(v), stroke: self.rule, "stroke-width": 1 }, self.g);
      el("text", { x: -9, y: self.y(v) + 3.5, "text-anchor": "end", fill: self.muted, "font-size": 11 }, self.g)
        .textContent = o.yfmt ? o.yfmt(v) : v;
    });
    xs.forEach(function (v) {
      el("text", { x: self.x(v), y: self.ih + 18, "text-anchor": "middle", fill: self.muted, "font-size": 11 }, self.g)
        .textContent = o.xfmt ? o.xfmt(v) : v;
    });
    el("line", { x1: 0, x2: self.iw, y1: self.ih, y2: self.ih, stroke: self.muted, "stroke-width": 1 }, this.g);
    if (xlab) el("text", { x: this.iw / 2, y: this.ih + 36, "text-anchor": "middle", fill: this.muted, "font-size": 12 }, this.g).textContent = xlab;
    if (ylab) el("text", { transform: "translate(" + (-this.m.l + 14) + "," + this.ih / 2 + ") rotate(-90)", "text-anchor": "middle", fill: this.muted, "font-size": 12 }, this.g).textContent = ylab;
    return this;
  };
  Frame.prototype.band = function (pts, colour) {
    var self = this, d = pts.map(function (p, i) { return (i ? "L" : "M") + self.x(p[0]) + "," + self.y(p[1]); }).join("");
    for (var i = pts.length - 1; i >= 0; i--) d += "L" + self.x(pts[i][0]) + "," + self.y(pts[i][2]);
    el("path", { d: d + "Z", fill: colour, "fill-opacity": .14, stroke: "none" }, this.g);
    return this;
  };
  Frame.prototype.line = function (pts, colour, o) {
    o = o || {};
    var self = this;
    var d = pts.map(function (p, i) { return (i ? "L" : "M") + self.x(p[0]) + "," + self.y(p[1]); }).join("");
    var p = el("path", { d: d, fill: "none", stroke: colour, "stroke-width": o.w || 1.9, "stroke-linejoin": "round", "stroke-dasharray": o.dash || null }, this.g);
    if (o.markers) pts.forEach(function (pt) {
      var c = el("circle", { cx: self.x(pt[0]), cy: self.y(pt[1]), r: 3.4, fill: colour }, self.g);
      if (o.tip) { c.style.cursor = "pointer"; c.addEventListener("mousemove", function (e) { showTip(e, o.tip(pt)); }); c.addEventListener("mouseleave", hideTip); }
    });
    return p;
  };
  Frame.prototype.hline = function (v, colour) {
    el("line", { x1: 0, x2: this.iw, y1: this.y(v), y2: this.y(v), stroke: colour || this.muted, "stroke-width": 1.2 }, this.g);
    return this;
  };
  Frame.prototype.vspan = function (a, b, colour, label) {
    el("rect", { x: this.x(a), y: 0, width: Math.abs(this.x(b) - this.x(a)), height: this.ih, fill: colour, "fill-opacity": .09 }, this.g);
    if (label) el("text", { x: (this.x(a) + this.x(b)) / 2, y: 12, "text-anchor": "middle", fill: this.muted, "font-size": 10 }, this.g).textContent = label;
    return this;
  };
  Frame.prototype.dots = function (pts, colour, o) {
    o = o || {};
    var self = this;
    pts.forEach(function (p) {
      var c = el("circle", { cx: self.x(p[0]), cy: self.y(p[1]), r: o.r || 3.2, fill: p[3] || colour, "fill-opacity": o.alpha || .8 }, self.g);
      if (o.tip) { c.style.cursor = "pointer"; c.addEventListener("mousemove", function (e) { showTip(e, o.tip(p)); }); c.addEventListener("mouseleave", hideTip); }
    });
    return this;
  };

  function barsH(host, rows, o) {
    o = o || {};
    var f = new Frame(host, { height: o.height || (rows.length * 30 + 56), margin: { t: 10, r: 24, b: 38, l: o.left || 170 } });
    var max = o.max !== undefined ? o.max : Math.max.apply(null, rows.map(function (r) { return r.v; }));
    var min = o.min !== undefined ? o.min : 0;
    f.scales([min, max], [0, rows.length]);
    var bh = f.ih / rows.length * 0.66;
    rows.forEach(function (r, i) {
      var y = f.ih - (i + 0.5) / rows.length * f.ih - bh / 2;
      var rect = el("rect", { x: 0, y: y, width: Math.max(1, f.x(r.v)), height: bh, fill: r.c || css("--accent"), rx: 3 }, f.g);
      rect.style.cursor = "pointer";
      rect.addEventListener("mousemove", function (e) { showTip(e, r.tip || (r.k + "<br>" + nice(r.v))); });
      rect.addEventListener("mouseleave", hideTip);
      el("text", { x: -9, y: y + bh / 2 + 4, "text-anchor": "end", fill: f.ink, "font-size": 11.5 }, f.g).textContent = r.k;
      el("text", { x: f.x(r.v) + 6, y: y + bh / 2 + 4, fill: f.muted, "font-size": 11, "font-family": "var(--mono)" }, f.g)
        .textContent = r.label !== undefined ? r.label : nice(r.v);
    });
    ticks(min, max, 5).forEach(function (v) {
      el("text", { x: f.x(v), y: f.ih + 18, "text-anchor": "middle", fill: f.muted, "font-size": 11 }, f.g).textContent = v;
    });
    if (o.xlab) el("text", { x: f.iw / 2, y: f.ih + 34, "text-anchor": "middle", fill: f.muted, "font-size": 12 }, f.g).textContent = o.xlab;
    return f;
  }

  function heatmap(host, cfg) {
    var rows = cfg.rows, cols = cfg.cols, val = cfg.value;
    var left = cfg.left || 168, bottom = cfg.bottom || 96;
    var f = new Frame(host, { height: rows.length * 26 + bottom + 24, margin: { t: 8, r: 70, b: bottom, l: left } });
    var cw = f.iw / cols.length, ch = f.ih / rows.length;
    var lo = cfg.min, hi = cfg.max;
    rows.forEach(function (r, i) {
      cols.forEach(function (c, j) {
        var v = val(r, c);
        var rect = el("rect", { x: j * cw, y: i * ch, width: cw - 1, height: ch - 1, fill: v === null ? css("--line") : colourOf(v, lo, hi) }, f.g);
        rect.style.cursor = "pointer";
        rect.addEventListener("mousemove", function (e) { showTip(e, cfg.tip(r, c, v)); });
        rect.addEventListener("mouseleave", hideTip);
        if (cfg.onclick) { rect.addEventListener("click", function () { cfg.onclick(r, c, v); }); }
      });
      el("text", { x: -8, y: i * ch + ch / 2 + 4, "text-anchor": "end", fill: f.ink, "font-size": 11 }, f.g).textContent = cfg.rowLabel ? cfg.rowLabel(r) : r;
    });
    cols.forEach(function (c, j) {
      el("text", { transform: "translate(" + (j * cw + cw / 2) + "," + (f.ih + 8) + ") rotate(55)", fill: f.muted, "font-size": 10 }, f.g)
        .textContent = cfg.colLabel ? cfg.colLabel(c) : c;
    });
    var lg = el("g", { transform: "translate(" + (f.iw + 14) + ",0)" }, f.g);
    for (var s = 0; s < 40; s++) {
      el("rect", { x: 0, y: s / 40 * f.ih, width: 11, height: f.ih / 40 + 1, fill: colourOf(hi - (s / 39) * (hi - lo), lo, hi) }, lg);
    }
    el("text", { x: 15, y: 9, fill: f.muted, "font-size": 10 }, lg).textContent = nice(hi, 2);
    el("text", { x: 15, y: f.ih, fill: f.muted, "font-size": 10 }, lg).textContent = nice(lo, 2);
    return f;
  }
  function colourOf(v, lo, hi) {
    var t = Math.max(0, Math.min(1, (v - lo) / (hi - lo)));
    var stops = [[68, 1, 84], [59, 82, 139], [33, 145, 140], [94, 201, 98], [253, 231, 37]];
    var p = t * (stops.length - 1), i = Math.min(stops.length - 2, Math.floor(p)), f = p - i;
    var a = stops[i], b = stops[i + 1];
    return "rgb(" + Math.round(a[0] + (b[0] - a[0]) * f) + "," + Math.round(a[1] + (b[1] - a[1]) * f) + "," + Math.round(a[2] + (b[2] - a[2]) * f) + ")";
  }

  function legend(host, items, onToggle) {
    var box = document.createElement("div");
    box.className = "legend";
    items.forEach(function (it) {
      var s = document.createElement("span");
      s.innerHTML = '<i style="background:' + it.c + '"></i>' + it.k;
      s.addEventListener("click", function () { s.classList.toggle("off"); onToggle(it.k, !s.classList.contains("off")); });
      box.appendChild(s);
    });
    host.appendChild(box);
    return box;
  }

  function toCSV(columns, rows) {
    var esc = function (v) {
      v = v === null || v === undefined ? "" : String(v);
      return /[",\n;]/.test(v) ? '"' + v.replace(/"/g, '""') + '"' : v;
    };
    return [columns.map(esc).join(",")].concat(rows.map(function (r) { return r.map(esc).join(","); })).join("\n");
  }
  function download(name, text, mime) {
    var b = new Blob([text], { type: mime || "text/csv;charset=utf-8" });
    var a = document.createElement("a");
    a.href = URL.createObjectURL(b); a.download = name; a.click();
    setTimeout(function () { URL.revokeObjectURL(a.href); }, 1500);
  }
  function savePNG(host, name, scale) {
    var svg = host.querySelector("svg");
    if (!svg) return;
    var clone = svg.cloneNode(true);
    clone.setAttribute("xmlns", NS);
    clone.querySelectorAll("text").forEach(function (t) {
      t.setAttribute("font-family", "Helvetica, Arial, sans-serif");
      if (!t.getAttribute("fill")) t.setAttribute("fill", css("--ink"));
    });
    var vb = (svg.getAttribute("viewBox") || "").split(/\s+/).map(Number);
    var w = vb[2] || svg.clientWidth, h = vb[3] || svg.clientHeight, k = scale || 3;
    var bg = el("rect", { x: 0, y: 0, width: w, height: h, fill: css("--bg") });
    clone.insertBefore(bg, clone.firstChild);
    var src = new XMLSerializer().serializeToString(clone);
    var img = new Image();
    img.onload = function () {
      var cv = document.createElement("canvas");
      cv.width = w * k; cv.height = h * k;
      var cx = cv.getContext("2d"); cx.scale(k, k); cx.drawImage(img, 0, 0);
      cv.toBlob(function (b) {
        var a = document.createElement("a");
        a.href = URL.createObjectURL(b); a.download = name; a.click();
        setTimeout(function () { URL.revokeObjectURL(a.href); }, 1500);
      });
    };
    img.src = "data:image/svg+xml;base64," + btoa(unescape(encodeURIComponent(src)));
  }

  global.GB = {
    Frame: Frame, barsH: barsH, heatmap: heatmap, legend: legend, el: el, css: css,
    nice: nice, ticks: ticks, toCSV: toCSV, download: download, savePNG: savePNG,
    showTip: showTip, hideTip: hideTip, colourOf: colourOf
  };
})(window);

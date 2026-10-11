"""Build native Word equation objects (OMML) from a restricted LaTeX subset.

Word stores equations as Office MathML, so an equation written this way is a
real equation object: it renders in the equation editor and can be edited by a
co-author, unlike an image or a run of italic text.

Supported: \\frac{}{}, _{}, ^{}, \\sum_{}^{}, \\sqrt{}, \\overline{}, \\hat{},
\\text{}, \\mathrm{}, Greek names, and the relation and set symbols the article
needs. Anything outside that set is emitted as literal text rather than being
silently dropped.
"""
import re

M = "http://schemas.openxmlformats.org/officeDocument/2006/math"
W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"

SYM = {
    "alpha": "\u03b1", "beta": "\u03b2", "gamma": "\u03b3", "delta": "\u03b4",
    "epsilon": "\u03b5", "theta": "\u03b8", "kappa": "\u03ba", "lambda": "\u03bb",
    "mu": "\u03bc", "nu": "\u03bd", "pi": "\u03c0", "rho": "\u03c1",
    "sigma": "\u03c3", "tau": "\u03c4", "phi": "\u03c6", "chi": "\u03c7",
    "psi": "\u03c8", "omega": "\u03c9", "Delta": "\u0394", "Sigma": "\u03a3",
    "Pi": "\u03a0", "Lambda": "\u039b", "Gamma": "\u0393", "Theta": "\u0398",
    "varepsilon": "\u03b5", "varphi": "\u03c6", "vartheta": "\u03d1", "ell": "\u2113",
    "sim": "\u223c", "propto": "\u221d", "equiv": "\u2261", "forall": "\u2200",
    "exists": "\u2203", "circ": "\u2218", "ast": "\u2217",
    "in": "\u2208", "notin": "\u2209", "subset": "\u2282", "subseteq": "\u2286",
    "cup": "\u222a", "cap": "\u2229", "emptyset": "\u2205",
    "bigcup": "\u22c3", "bigcap": "\u22c2", "setminus": "\u2216", "quad": "   ",
    "leq": "\u2264", "geq": "\u2265", "neq": "\u2260", "approx": "\u2248",
    "pm": "\u00b1", "times": "\u00d7", "cdot": "\u00b7",
    "to": "\u2192", "rightarrow": "\u2192", "infty": "\u221e",
    "ldots": "\u2026", "mid": "|", "uparrow": "\u2191", "downarrow": "\u2193",
}


def _esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _run(text, italic=True):
    if text == "":
        return ""
    props = "" if italic else '<m:rPr><m:sty m:val="p"/></m:rPr>'
    return '<m:r>' + props + '<m:t xml:space="preserve">' + _esc(text) + '</m:t></m:r>'


class _P:
    def __init__(self, s):
        self.s = s
        self.i = 0

    def eof(self):
        return self.i >= len(self.s)

    def peek(self):
        return self.s[self.i] if not self.eof() else ""

    def group(self):
        while self.peek() == " ":
            self.i += 1
        if self.peek() == "{":
            depth, start = 0, self.i
            while not self.eof():
                if self.s[self.i] == "{":
                    depth += 1
                elif self.s[self.i] == "}":
                    depth -= 1
                    if depth == 0:
                        inner = self.s[start + 1:self.i]
                        self.i += 1
                        return parse(inner)
                self.i += 1
            return parse(self.s[start + 1:])
        if self.peek() == "\\":
            j = self.i
            self.i += 1
            while not self.eof() and self.peek().isalpha():
                self.i += 1
            return parse(self.s[j:self.i])
        ch = self.peek()
        self.i += 1
        return _run(ch)


def parse(s):
    out = []
    p = _P(s)
    buf = ""

    def flush():
        nonlocal buf
        if buf:
            out.append(_run(buf))
            buf = ""

    while not p.eof():
        c = p.peek()
        if c == "\\":
            flush()
            j = p.i + 1
            k = j
            while k < len(s) and s[k].isalpha():
                k += 1
            name = s[j:k]
            p.i = k if k > j else j + 1
            if name == "frac":
                num, den = p.group(), p.group()
                out.append('<m:f><m:fPr><m:type m:val="bar"/></m:fPr>'
                           '<m:num>' + num + '</m:num><m:den>' + den + '</m:den></m:f>')
            elif name == "binom":
                top, bot = p.group(), p.group()
                out.append('<m:d><m:dPr><m:begChr m:val="("/><m:endChr m:val=")"/>'
                           '</m:dPr><m:e><m:f><m:fPr><m:type m:val="noBar"/></m:fPr>'
                           '<m:num>' + top + '</m:num><m:den>' + bot +
                           '</m:den></m:f></m:e></m:d>')
            elif name == "underline":
                e = p.group()
                out.append('<m:bar><m:barPr><m:pos m:val="bot"/></m:barPr>'
                           '<m:e>' + e + '</m:e></m:bar>')
            elif name == "sqrt":
                e = p.group()
                out.append('<m:rad><m:radPr><m:degHide m:val="1"/></m:radPr>'
                           '<m:deg/><m:e>' + e + '</m:e></m:rad>')
            elif name in ("text", "mathrm", "operatorname"):
                e = p.group()
                out.append(e.replace("<m:r>", '<m:r><m:rPr><m:sty m:val="p"/></m:rPr>'))
            elif name in ("overline", "bar"):
                e = p.group()
                out.append('<m:bar><m:barPr><m:pos m:val="top"/></m:barPr>'
                           '<m:e>' + e + '</m:e></m:bar>')
            elif name == "hat":
                e = p.group()
                out.append('<m:acc><m:accPr><m:chr m:val="\u0302"/></m:accPr>'
                           '<m:e>' + e + '</m:e></m:acc>')
            elif name == "int":
                sub = sup = ""
                while p.peek() in ("_", "^"):
                    which = p.peek()
                    p.i += 1
                    g = p.group()
                    if which == "_":
                        sub = g
                    else:
                        sup = g
                rest = parse(s[p.i:])
                p.i = len(s)
                hide = ""
                if not sub:
                    hide += '<m:subHide m:val="1"/>'
                if not sup:
                    hide += '<m:supHide m:val="1"/>'
                out.append('<m:nary><m:naryPr><m:chr m:val="\u222b"/>'
                           '<m:limLoc m:val="subSup"/><m:grow m:val="1"/>' + hide +
                           '</m:naryPr><m:sub>' + sub + '</m:sub><m:sup>' + sup +
                           '</m:sup><m:e>' + rest + '</m:e></m:nary>')
            elif name == "sum":
                sub = sup = ""
                while p.peek() in ("_", "^"):
                    which = p.peek()
                    p.i += 1
                    g = p.group()
                    if which == "_":
                        sub = g
                    else:
                        sup = g
                rest = parse(s[p.i:])
                p.i = len(s)
                hide = ""
                if not sub:
                    hide += '<m:subHide m:val="1"/>'
                if not sup:
                    hide += '<m:supHide m:val="1"/>'
                out.append('<m:nary><m:naryPr><m:chr m:val="\u2211"/>'
                           '<m:limLoc m:val="undOvr"/><m:grow m:val="1"/>' + hide +
                           '</m:naryPr><m:sub>' + sub + '</m:sub><m:sup>' + sup +
                           '</m:sup><m:e>' + rest + '</m:e></m:nary>')
            elif name in ("left", "right", "big", "Big"):
                continue
            elif name in SYM:
                out.append(_run(SYM[name], italic=False))
            elif name == "":
                # A backslash followed by punctuation: the LaTeX spacing commands
                # \, \; \: \! are thin spaces, not literal characters.
                ch = s[p.i - 1] if p.i - 1 < len(s) else ""
                out.append(_run(" " if ch in ",;:!" else ch, italic=False))
            else:
                out.append(_run(name, italic=False))
        elif c in ("_", "^"):
            if buf:
                out.append(_run(buf[:-1]))
                base = _run(buf[-1])
                buf = ""
            else:
                base = out.pop() if out else _run("")
            p.i += 1
            first = p.group()
            second = None
            if p.peek() in ("_", "^"):
                nxt = p.peek()
                p.i += 1
                second = p.group()
            if second is None:
                tag = "sSub" if c == "_" else "sSup"
                inner = "sub" if c == "_" else "sup"
                out.append('<m:' + tag + '><m:e>' + base + '</m:e><m:' + inner + '>' +
                           first + '</m:' + inner + '></m:' + tag + '>')
            else:
                sub, sup = (first, second) if c == "_" else (second, first)
                out.append('<m:sSubSup><m:e>' + base + '</m:e><m:sub>' + sub +
                           '</m:sub><m:sup>' + sup + '</m:sup></m:sSubSup>')
        elif c == "{":
            flush()
            out.append(p.group())
        elif c == "}":
            p.i += 1
        elif c.isdigit() or c in "+-=()[]|,.:;/<>!? ":
            flush()
            j = p.i
            while not p.eof() and (p.peek().isdigit() or p.peek() in "+-=()[]|,.:;/<>!? "):
                p.i += 1
            out.append(_run(s[j:p.i], italic=False))
        else:
            buf += c
            p.i += 1
    flush()
    return "".join(out)


def omml_xml(latex):
    return '<m:oMath xmlns:m="' + M + '" xmlns:w="' + W + '">' + parse(latex.strip()) + '</m:oMath>'


def add_equation(doc, latex, number=None, indent_in=0.25):
    """Append a display equation; the number sits at the right margin."""
    from docx.oxml import parse_xml
    from docx.enum.text import WD_TAB_ALIGNMENT
    from docx.shared import Inches, Pt

    try:
        p = doc.add_paragraph(style="Equation")   # the template's own equation style
    except KeyError:
        p = doc.add_paragraph()
    pf = p.paragraph_format
    pf.left_indent = Inches(indent_in)
    pf.space_before = Pt(9)
    pf.space_after = Pt(9)
    pf.tab_stops.add_tab_stop(Inches(3.33), WD_TAB_ALIGNMENT.RIGHT)
    p._p.append(parse_xml(omml_xml(latex)))
    if number is not None:
        r = p.add_run("\t(" + str(number) + ")")
        r.italic = False
    return p

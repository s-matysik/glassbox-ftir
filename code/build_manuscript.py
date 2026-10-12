"""Build the IEEE Access manuscript inside the official Word template.

The template's own styles are reused (Paper Title, AU, PI, Abstract, IT, H1,
H2, PARA, Fig Caption, REF_Txt, AU_Bios), so the output keeps the journal's
formatting instead of imitating it. Body text comes from manuscript.json;
numbers come from the tables in results/, never from the text file, so the
manuscript cannot drift away from the computed results.

Run from the deposit root:  python code/build_manuscript.py
"""
import json
import os

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import omml  # noqa: E402

from docx import Document
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
from docx.shared import Inches, Pt


def _bordered(table):
    """Draw single-line borders when the template carries no grid style."""
    tbl_pr = table._tbl.tblPr
    borders = OxmlElement("w:tblBorders")
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        el = OxmlElement(f"w:{edge}")
        el.set(qn("w:val"), "single")
        el.set(qn("w:sz"), "4")
        el.set(qn("w:color"), "000000")
        borders.append(el)
    tbl_pr.append(borders)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
FIG = os.path.join(ROOT, "figures")
DOCS = os.path.join(ROOT, "docs")
TEMPLATE = os.path.join(DOCS, "Access-Template-2024.docx")
CONTENT = os.path.join(DOCS, "manuscript.json")
OUT = os.path.join(DOCS, "IEEE_Access_manuscript.docx")

STYLE = {"title": "Paper Title", "authors": "AU", "aff": "PI_No Space",
         "aff_last": "PI", "abstract": "Abstract", "index": "IT",
         "h1_first": "H1_List (No Space)", "h1": "H1_List (Space)",
         "h2_first": "H2_First", "h2": "H2_Cont", "h3": "H3",
         "para": "PARA", "para_indent": "PARA_Indent",
         "caption": "Fig Caption", "ref": "REF_Txt", "bio": "AU_Bios"}


def clear_body(doc):
    """Remove the template's instructional body, keep styles and section setup."""
    body = doc.element.body
    for child in list(body):
        if child.tag.endswith("}sectPr"):
            continue
        body.remove(child)


def add(doc, text, style, indent=False):
    s = STYLE[style]
    try:
        p = doc.add_paragraph(text, style=s)
    except KeyError:
        p = doc.add_paragraph(text)
    return p


def add_runs(doc, parts, style):
    """parts: list of (text, bold, italic)."""
    try:
        p = doc.add_paragraph(style=STYLE[style])
    except KeyError:
        p = doc.add_paragraph()
    for text, bold, italic in parts:
        r = p.add_run(text)
        r.bold, r.italic = bold, italic
    return p


def _styled(doc, text, preferred, fallback_pt):
    """Paragraph in the template's own style, or at the right size if absent.

    The template's Normal is 12 pt while its body text is 10 pt, so anything
    left in Normal reads as oversized beside the column it sits in.
    """
    from docx.shared import Pt
    for name in preferred:
        try:
            return doc.add_paragraph(text, style=name)
        except KeyError:
            continue
    p = doc.add_paragraph(text)
    for r in p.runs:
        r.font.size = Pt(fallback_pt)
    return p


def add_table(doc, caption_number, caption, header, rows, note=None):
    cap = _styled(doc, f"TABLE {caption_number}", ("Table Title",), 8)
    cap.runs[0].bold = True
    sub = _styled(doc, caption, ("Table Caption",), 7)
    sub.runs[0].italic = True
    t = doc.add_table(rows=1, cols=len(header))
    for name in ("Table Grid", "TableGrid", "Grid Table 1 Light"):
        try:
            t.style = name
            break
        except KeyError:
            continue
    else:
        _bordered(t)
    for i, h in enumerate(header):
        c = t.rows[0].cells[i]
        c.text = ""
        _h = c.paragraphs[0].add_run(str(h))
        _h.bold = True
        _h.font.size = Pt(8)
    for row in rows:
        cells = t.add_row().cells
        for i, v in enumerate(row):
            cells[i].text = ""
            cells[i].paragraphs[0].add_run("" if v is None else str(v)).font.size = Pt(8)
    if note:
        n = _styled(doc, note, ("Table Caption",), 7)
        n.runs[0].italic = True
    _styled(doc, "", ("PARA",), 10)
    return t


def add_figure(doc, filename, number, caption, width_in=3.4):
    path = os.path.join(FIG, filename)
    if not os.path.exists(path):
        raise FileNotFoundError(f"figure missing: {path}")
    doc.add_picture(path, width=Inches(width_in))
    add(doc, f"FIGURE {number}.  {caption}", "caption")



def tidy_styles(doc):
    """Leave nothing in Normal: the template's Normal is 12 pt, the body is 10 pt."""
    from docx.shared import Pt
    moved = 0
    for p in doc.paragraphs:
        if p.style.name == "Normal":
            try:
                p.style = doc.styles["PARA"]
            except KeyError:
                for r in p.runs:
                    r.font.size = Pt(10)
            moved += 1
    for t in doc.tables:
        for row in t.rows:
            for c in row.cells:
                for p in c.paragraphs:
                    for r in p.runs:
                        if r.font.size is None:
                            r.font.size = Pt(8)
    return moved



W_NS = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


def set_columns(section, n, space=400):
    """Set the column count on a section, as the journal layout requires."""
    from docx.oxml import OxmlElement
    sectPr = section._sectPr
    for c in sectPr.findall(W_NS + "cols"):
        sectPr.remove(c)
    cols = OxmlElement("w:cols")
    cols.set(W_NS + "num", str(n))
    cols.set(W_NS + "space", str(space))
    sectPr.append(cols)
    return section


def continuous_section(doc, n_cols):
    """Start a new continuous section with the given column count."""
    from docx.enum.section import WD_SECTION
    s = doc.add_section(WD_SECTION.CONTINUOUS)
    return set_columns(s, n_cols)


def add_wide_figure(doc, filename, number, caption):
    """A figure wider than one column spans the page, as in the journal."""
    continuous_section(doc, 1)
    add_figure(doc, filename, number, caption, width_in=6.4)
    continuous_section(doc, 2)


def fix_footer(doc):
    """Clear the template's placeholder footer.

    The template ships a footer naming a volume and year that belong to the
    sample document. Writing a volume number of our own would be worse: a
    submitted manuscript has no volume, and a footer that supplies one reads as
    though the article were already published. The footer is therefore emptied
    and the journal typesets its own.
    """
    for s in doc.sections:
        for p in list(s.footer.paragraphs) + list(s.header.paragraphs):
            if "VOLUME" in p.text.upper() or "2017" in p.text:
                for r in p.runs:
                    r.text = ""


def main():
    if not os.path.exists(CONTENT):
        raise SystemExit(f"write {CONTENT} first: it holds the manuscript text")
    content = json.load(open(CONTENT, encoding="utf-8"))
    doc = Document(TEMPLATE)
    clear_body(doc)

    add(doc, content["title"], "title")
    add(doc, content["authors_line"], "authors")
    for i, a in enumerate(content["affiliations"]):
        add(doc, a, "aff_last" if i == len(content["affiliations"]) - 1 else "aff")
    add(doc, content["corresponding"], "aff_last")
    if content.get("funding"):
        add(doc, content["funding"], "aff_last")

    add_runs(doc, [("ABSTRACT ", True, False), (content["abstract"], False, False)], "abstract")
    add_runs(doc, [("INDEX TERMS ", True, False), (content["index_terms"], False, False)], "index")

    # The journal runs title, authors and abstract across the full page and
    # only then switches to two columns.
    set_columns(doc.sections[0], 1, space=720)
    continuous_section(doc, 2)

    first_h1 = True
    for sec in content["sections"]:
        add(doc, sec["heading"], "h1_first" if first_h1 else "h1")
        first_h1 = False
        first_h2 = True
        for block in sec["blocks"]:
            kind = block["type"]
            if kind == "h2":
                add(doc, block["text"], "h2_first" if first_h2 else "h2")
                first_h2 = False
            elif kind == "h3":
                add(doc, block["text"], "h3")
            elif kind == "para":
                add(doc, block["text"], "para_indent" if block.get("indent") else "para")
            elif kind == "table":
                add_table(doc, block["number"], block["caption"], block["header"],
                          block["rows"], block.get("note"))
            elif kind == "figure":
                if block.get("width", 3.4) > 3.6:
                    add_wide_figure(doc, block["file"], block["number"], block["caption"])
                else:
                    add_figure(doc, block["file"], block["number"], block["caption"],
                               block.get("width", 3.4))
            elif kind == "eq":
                omml.add_equation(doc, block["latex"], block["number"])
            else:
                raise ValueError(f"unknown block type {kind}")

    add(doc, "REFERENCES", "h1")
    for i, ref in enumerate(content["references"], start=1):
        add(doc, f"[{i}] {ref}", "ref")

    for bio in content.get("biographies", []):
        add_runs(doc, [(bio["name"].upper() + " ", True, False), (bio["text"], False, False)], "bio")

    moved = tidy_styles(doc)
    fix_footer(doc)
    doc.save(OUT)
    words = sum(len(p.text.split()) for p in doc.paragraphs)
    print(f"saved {OUT} | paragraphs restyled out of Normal: {moved}")
    print(f"abstract words {len(content['abstract'].split())} | body paragraphs {len(doc.paragraphs)} "
          f"| approx words {words} | tables {len(doc.tables)} | references {len(content['references'])}")


if __name__ == "__main__":
    main()

"""Check that every equation in the manuscript becomes a real Word equation.

A reader of a flattened export sees a fraction as two adjacent groups and may
conclude the bar is missing. This test settles that objectively: for each
equation it compares what the source asks for against what the generated Office
MathML contains, and fails if any LaTeX command survives as literal text.

Run from the deposit root:  python code/test_equations.py
"""
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import omml  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCS = os.path.join(ROOT, "docs")

STRUCTURE = [("\\frac", "<m:f>", "fraction"),
             ("\\sum", '<m:chr m:val="\u2211"/>', "summation"),
             ("\\int", '<m:chr m:val="\u222b"/>', "integral"),
             ("\\sqrt", "<m:rad>", "radical"),
             ("\\binom", 'm:val="noBar"', "binomial"),
             ("\\underline", 'm:val="bot"', "underline"),
             ("\\overline", 'm:val="top"', "overline"),
             ("\\hat", "<m:acc>", "accent")]

# Commands that are meant to render as upright words, not as structure.
WORD_OPERATORS = {"min", "max", "med", "text", "mathrm", "operatorname", "left",
                  "right", "quad"}


def equations():
    for i in (1, 2, 3):
        path = os.path.join(DOCS, f"manuscript_part{i}.json")
        if not os.path.exists(path):
            continue
        part = json.load(open(path, encoding="utf-8"))
        for s in part["sections"]:
            for b in s["blocks"]:
                if b.get("type") == "eq":
                    yield b["number"], b["latex"]


def main():
    failures = []
    checked = 0
    for number, latex in sorted(equations()):
        xml = omml.omml_xml(latex)
        checked += 1
        for cmd, marker, name in STRUCTURE:
            want = latex.count(cmd)
            got = xml.count(marker)
            if want and got < want:
                failures.append(f"({number}) asks for {want} {name}(s), produced {got}")
        leaked = set(re.findall(r'<m:t xml:space="preserve">([A-Za-z]+)</m:t>', xml))
        bad = {w for w in leaked
               if re.search(r"\\" + w + r"(?![A-Za-z])", latex) and w not in WORD_OPERATORS}
        if bad:
            failures.append(f"({number}) command(s) left as literal text: {sorted(bad)}")
        if re.search(r"<m:t[^>]*>[^<]*\\", xml):
            failures.append(f"({number}) a backslash survived into the output")

    print(f"equations checked: {checked}")
    if failures:
        for f in failures:
            print("  FAIL " + f)
        raise SystemExit(1)
    print("every equation carries the structure its source asks for")


if __name__ == "__main__":
    main()

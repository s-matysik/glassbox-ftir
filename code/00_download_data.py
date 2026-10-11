"""Fetch the four open datasets this deposit is built on.

Everything comes from Zenodo under CC BY 4.0 and is addressed by record id, so
the download is reproducible and version-pinned. Re-running skips files that
are already present.

Run from the deposit root:  python code/00_download_data.py
"""
import json
import os
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data_raw")

RECORDS = [
    ("21822279", "ftir",
     "FTIR spectra of plastic recovered from organic fertiliser, compost-amended "
     "soil and wastewater-irrigated soil. Concepts 1 and 2."),
    ("20767361", "industrial_co2",
     "Measured versus modelled CO2 for industrial facilities: EUTL and EPA GHGRP "
     "against Climate TRACE. Concept 3."),
    ("21740286", "choice_experiment",
     "Discrete choice experiment on bioplastic food packaging, Germany. "
     "Screening of survey instruments."),
    ("3887306", "buildings",
     "Building Data Genome Project 2: hourly meters and metadata. Concept 4. "
     "Single archive of about 595 MB."),
]


def fetch(record_id, folder, note):
    out = os.path.join(RAW, folder)
    os.makedirs(out, exist_ok=True)
    url = f"https://zenodo.org/api/records/{record_id}"
    rec = json.load(urllib.request.urlopen(url, timeout=120))
    with open(os.path.join(out, "zenodo_record.json"), "w", encoding="utf-8") as fh:
        json.dump({"id": record_id, "title": rec["metadata"]["title"],
                   "doi": rec["metadata"].get("doi"),
                   "license": (rec["metadata"].get("license") or {}).get("id"),
                   "publication_date": rec["metadata"].get("publication_date"),
                   "note": note,
                   "files": [{"key": f["key"], "size": f["size"]} for f in rec["files"]]},
                  fh, ensure_ascii=False, indent=1)
    for f in rec["files"]:
        name = os.path.basename(f["key"])
        dst = os.path.join(out, name)
        if os.path.exists(dst) and os.path.getsize(dst) == f["size"]:
            continue
        urllib.request.urlretrieve(f["links"]["self"], dst)
    got = [n for n in os.listdir(out) if n != "zenodo_record.json"]
    size = sum(os.path.getsize(os.path.join(out, n)) for n in got)
    print(f"{folder}: {len(got)} files, {size / 1e6:.1f} MB")


def main():
    for record_id, folder, note in RECORDS:
        fetch(record_id, folder, note)


if __name__ == "__main__":
    main()

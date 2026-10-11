"""The entry condition: between-specimen spread against instrument noise.

The method is only applicable when the variation to be explained is large
relative to what the instrument contributes. This script measures both on the
primary dataset and writes results/entry_condition_primary.csv, which Fig. 2,
the supplement and the site all read.

Note on what the ratio is: it compares the between-specimen spread with the
NOISE estimate, not with the detection threshold. The threshold is k times a
noise standard deviation, with k fixed at 4 by the protocol, so dividing these
values by 4 gives the headroom against the threshold.

Run from the deposit root:  python code/09_entry_condition.py
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plastics_data as pdata  # noqa: E402
from ara_aco import instrument_noise, preprocess  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
K = 4.0


def main():
    os.makedirs(RES, exist_ok=True)
    grid, Y, meta = pdata.load("c8")
    y = meta["Polymer"].to_numpy().astype(str)
    Yp = preprocess(Y)
    sigma = instrument_noise(Yp)

    rows = []
    for polymer in sorted(set(y)):
        m = y == polymer
        spread = float(Yp[m].std(axis=0, ddof=1).mean())
        noise = float(np.median(sigma[m]))
        rows.append({"dataset": "FTIR-Plastics C8", "stratum": polymer,
                     "n": int(m.sum()), "sd_between": round(spread, 4),
                     "noise_sd": round(noise, 4),
                     "spread_over_noise": round(spread / noise, 1),
                     "spread_over_threshold": round(spread / (K * noise), 1)})
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(RES, "entry_condition_primary.csv"), index=False,
               encoding="utf-8-sig")
    print(out.to_string(index=False))
    print(f"spread over noise: {out.spread_over_noise.min():.1f} to "
          f"{out.spread_over_noise.max():.1f}")
    print(f"spread over threshold (k = {K:g}): {out.spread_over_threshold.min():.1f} to "
          f"{out.spread_over_threshold.max():.1f}")


if __name__ == "__main__":
    main()

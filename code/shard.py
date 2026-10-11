"""Parallelism by plain subprocesses.

The execution environment refuses ``CreateNamedPipe``, so ``multiprocessing``
and the joblib ``loky`` backend cannot start (``PermissionError: WinError 5``).
Each script therefore re-runs itself as N independent worker processes that
communicate only through files: the parent sets ``ARA_SHARD_INDEX`` and
``ARA_SHARD_COUNT``, each worker computes the jobs whose position is congruent
to its index, pickles its results into ``logs/`` and exits, and the parent
concatenates the pickles in shard order. Job order, and therefore every
result, is independent of how many workers are used.
"""
from __future__ import annotations

import os
import pickle
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOGDIR = os.path.join(ROOT, "logs")


def is_worker() -> bool:
    return os.environ.get("ARA_SHARD_INDEX") is not None


def worker_id() -> tuple[int, int]:
    return int(os.environ["ARA_SHARD_INDEX"]), int(os.environ["ARA_SHARD_COUNT"])


def shard_file(tag: str, index: int) -> str:
    os.makedirs(LOGDIR, exist_ok=True)
    return os.path.join(LOGDIR, f"_shard_{tag}_{index}.pkl")


def my_jobs(jobs: list, index: int, count: int) -> list[tuple[int, object]]:
    return [(i, j) for i, j in enumerate(jobs) if i % count == index]


def dump(tag: str, index: int, payload) -> None:
    with open(shard_file(tag, index), "wb") as fh:
        pickle.dump(payload, fh)


def run(script_path: str, tag: str, n_workers: int, extra_env: dict | None = None,
        threads: int = 2) -> list:
    """Spawn ``n_workers`` copies of ``script_path`` and return the merged,
    order-restored list of their results."""
    os.makedirs(LOGDIR, exist_ok=True)
    for i in range(n_workers):
        f = shard_file(tag, i)
        if os.path.exists(f):
            os.remove(f)
    procs = []
    for i in range(n_workers):
        env = dict(os.environ)
        env.update(extra_env or {})
        env["ARA_SHARD_INDEX"] = str(i)
        env["ARA_SHARD_COUNT"] = str(n_workers)
        for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                  "NUMEXPR_NUM_THREADS"):
            env[v] = str(threads)
        out = open(os.path.join(LOGDIR, f"_shard_{tag}_{i}.out.txt"), "w", encoding="utf-8")
        procs.append((subprocess.Popen([sys.executable, script_path], cwd=ROOT,
                                       stdout=out, stderr=subprocess.STDOUT, env=env), out))
    failed = []
    for i, (p, fh) in enumerate(procs):
        rc = p.wait()
        fh.close()
        if rc != 0:
            failed.append((i, rc))
    if failed:
        raise RuntimeError(f"worker shards failed: {failed}; see logs/_shard_{tag}_*.out.txt")
    merged = []
    for i in range(n_workers):
        with open(shard_file(tag, i), "rb") as fh:
            merged.extend(pickle.load(fh))
    merged.sort(key=lambda kv: kv[0])
    return [v for _, v in merged]

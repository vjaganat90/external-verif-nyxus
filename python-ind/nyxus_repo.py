"""Locate a Nyxus checkout for independent-verification scripts.

The verifier scripts parse fixtures and current Nyxus regression constants out
of the Nyxus C++ test headers. This helper deliberately prefers an existing
checkout and does not fetch, pull, checkout, or reset it. That keeps an oracle
run from quietly changing the source tree being validated.

Resolution order:

1. Explicit path argument.
2. ``NYXUS_REPO_DIR``.
3. Common local locations such as ``../nyxus`` and an existing ``.cache/nyxus``.
4. Optional one-time clone into ``.cache/nyxus`` only when
   ``NYXUS_ALLOW_CLONE=1`` is set.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


NYXUS_REPO_URL = os.environ.get("NYXUS_REPO_URL", "https://github.com/vjaganat90/nyxus.git")

_EXTERNAL_VERIF_ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_CACHE_DIR = _EXTERNAL_VERIF_ROOT / ".cache" / "nyxus"
_REQUIRED_HEADER = Path("tests") / "test_data.h"


def _run(cmd: list[str], cwd: Path | None = None) -> None:
    print(f"[nyxus_repo] $ {' '.join(cmd)}", file=sys.stderr)
    subprocess.run(cmd, cwd=str(cwd) if cwd else None, check=True)


def _allow_clone() -> bool:
    return os.environ.get("NYXUS_ALLOW_CLONE", "").strip().lower() in {"1", "true", "yes", "on"}


def _looks_like_nyxus_checkout(path: Path) -> bool:
    return path.is_dir() and (path / _REQUIRED_HEADER).is_file()


def _require_checkout(path: Path, source: str) -> Path:
    if _looks_like_nyxus_checkout(path):
        return path.resolve()
    raise RuntimeError(
        f"{source} points to {path}, but it does not look like a Nyxus checkout "
        f"(missing {_REQUIRED_HEADER})."
    )


def _local_candidates() -> list[Path]:
    cwd = Path.cwd().resolve()
    return [
        _EXTERNAL_VERIF_ROOT.parent / "nyxus",
        cwd,
        cwd.parent / "nyxus",
        _DEFAULT_CACHE_DIR,
    ]


def ensure_nyxus_repo(explicit_path: str | os.PathLike[str] | None = None) -> Path:
    """Return a path to a usable Nyxus checkout without mutating it."""
    if explicit_path:
        path = Path(explicit_path).resolve()
        return _require_checkout(path, "explicit Nyxus path")

    env_dir = os.environ.get("NYXUS_REPO_DIR")
    if env_dir:
        path = Path(env_dir).resolve()
        return _require_checkout(path, "NYXUS_REPO_DIR")

    for candidate in _local_candidates():
        if _looks_like_nyxus_checkout(candidate):
            return candidate.resolve()

    if _allow_clone():
        cache_dir = _DEFAULT_CACHE_DIR
        if cache_dir.exists():
            return _require_checkout(cache_dir, "existing .cache/nyxus")
        cache_dir.parent.mkdir(parents=True, exist_ok=True)
        cmd = ["git", "clone", "--depth", "1"]
        branch = os.environ.get("NYXUS_BRANCH")
        if branch:
            cmd.extend(["--branch", branch])
        cmd.extend([NYXUS_REPO_URL, str(cache_dir)])
        _run(cmd)
        return _require_checkout(cache_dir, "newly cloned .cache/nyxus")

    searched = "\n  - ".join(str(path) for path in _local_candidates())
    raise RuntimeError(
        "Could not find a Nyxus checkout for V&V fixture parsing.\n"
        "Pass the repo path as the first script argument, set NYXUS_REPO_DIR, "
        "place a checkout at ../nyxus, or set NYXUS_ALLOW_CLONE=1 to clone once "
        f"into {_DEFAULT_CACHE_DIR}.\n"
        f"Searched:\n  - {searched}"
    )


if __name__ == "__main__":
    print(ensure_nyxus_repo(sys.argv[1] if len(sys.argv) > 1 else None))

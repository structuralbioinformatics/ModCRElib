"""
Small helpers for running external programs without shell redirection.

Used by profile modules (FIMO/TOMTOM wrappers, scan pipelines) so behavior does not
depend on the login shell, ``PATH``, or ``> file`` redirects.
"""

from __future__ import annotations

import os
import shlex
import shutil
import subprocess
from typing import Iterable, List, Optional, Sequence


def append_file(src: str, dst: str) -> None:
    """Append the contents of ``src`` onto ``dst`` using Python I/O (no ``cat``)."""
    parent = os.path.dirname(os.path.abspath(dst))
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(dst, "a", encoding="utf-8", errors="replace") as out:
        with open(src, "r", encoding="utf-8", errors="replace") as inp:
            shutil.copyfileobj(inp, out)


def run_argv(
    argv: Sequence[str],
    *,
    cwd: Optional[str] = None,
    env: Optional[dict] = None,
    timeout: Optional[float] = None,
) -> subprocess.CompletedProcess:
    """Run a command with a proper argument list (no shell)."""
    return subprocess.run(
        list(argv),
        cwd=cwd,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        timeout=timeout,
        check=False,
    )


def run_python_with_parameters(
    python: str,
    script: str,
    parameters: str,
    *,
    cwd: Optional[str] = None,
) -> subprocess.CompletedProcess:
    """Run ``python script`` with a whitespace-separated parameter string."""
    argv = [python, script] + shlex.split(parameters)
    return run_argv(argv, cwd=cwd)


def capture_stdout(argv: Sequence[str], *, env: Optional[dict] = None) -> List[str]:
    """Run ``argv`` and return stdout split into lines; raise on failure or empty output."""
    proc = run_argv(argv, env=env)
    if proc.returncode != 0:
        err = (proc.stderr or proc.stdout or "").strip()
        raise OSError("command failed (exit %s): %s\n%s" % (proc.returncode, argv[0], err))
    if not proc.stdout or not proc.stdout.strip():
        raise OSError("command produced empty stdout: %s" % argv[0])
    return proc.stdout.splitlines()


def read_first_existing_file(output_dir: str, basenames: Iterable[str]) -> List[str]:
    """Read the first non-empty file found under ``output_dir``."""
    for name in basenames:
        path = os.path.join(output_dir, name)
        if os.path.isfile(path) and os.path.getsize(path) > 0:
            with open(path, "r", encoding="utf-8", errors="replace") as fh:
                return fh.read().splitlines()
    raise OSError(
        "No output in %s (looked for %s)"
        % (output_dir, ", ".join(basenames))
    )

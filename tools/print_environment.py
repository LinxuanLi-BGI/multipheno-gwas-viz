#!/usr/bin/env python3
"""Print runtime versions and protocol source identifiers for a run log.

Usage: python tools/print_environment.py > run_log.txt
Uses only the Python standard library; Git information is optional.
"""

from datetime import datetime
import hashlib
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
import platform
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
PACKAGES = ('pandas', 'numpy', 'scipy', 'matplotlib')


def git_output(*args):
    """Read Git metadata from the protocol directory, regardless of cwd."""
    try:
        result = subprocess.run(
            ['git', '-C', str(ROOT), *args],
            capture_output=True, text=True, encoding='utf-8', errors='replace',
            timeout=10,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    return result.stdout.strip() if result.returncode == 0 else None


def main():
    print('Run environment and protocol versions')
    print(f'Recorded at: {datetime.now().astimezone().isoformat(timespec="seconds")}')
    print(f'Operating system: {platform.platform()}')
    print(f'Architecture: {platform.machine()}')
    print(f'Python: {platform.python_version()} ({platform.python_implementation()})')
    print(f'Python build: {sys.version.replace(chr(10), " ")}')
    print(f'Python executable: {sys.executable}')
    print('\nPackage versions (active Python environment):')
    for package in PACKAGES:
        try:
            installed = version(package)
        except PackageNotFoundError:
            installed = 'not installed'
        print(f'{package}: {installed}')

    print(f'\nProtocol directory: {ROOT}')
    # Avoid reporting an enclosing repository for a downloaded source copy.
    git_root = git_output('rev-parse', '--show-toplevel')
    commit = None
    if git_root and Path(git_root).resolve() == ROOT:
        commit = git_output('rev-parse', 'HEAD')
    if commit:
        print(f'Git commit: {commit}')
        print(f'Git revision: {git_output("describe", "--tags", "--always", "--dirty") or "unavailable"}')
        status = git_output('status', '--short', '--untracked-files=normal')
        if status is None:
            print('Git working tree: unavailable')
        elif status:
            print('Git working tree: local changes present')
            print(status)
        else:
            print('Git working tree: clean')
    else:
        print('Git commit/release: unavailable (Git missing or no Git checkout).')
        print('Record the downloaded release tag or source commit separately.')

if __name__ == '__main__':
    main()

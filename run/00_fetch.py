#!/usr/bin/env python3
"""Download the deposited files that are too large to track in the repository.

By default this fetches the per-tool point clouds into data/pointclouds/, which
is all that run/04_iou.py needs. With --frames it also fetches the archive of
the 500 ensemble frames, which is needed only to re-run the detectors
themselves; see tools/README.md.

Every file is checked against the checksum the deposit publishes for it, so a
truncated or corrupted download is reported instead of being analysed silently.
A file already present with the right checksum is left alone, which makes the
script cheap to re-run and safe to interrupt.

The deposit is named in the ``deposit`` block of data/parameters.yaml and
nowhere else.

Writes data/pointclouds/pc_GLUT1_<tool>.npz
"""
import argparse
import hashlib
import json
import os
import shutil
import sys
import tempfile
import urllib.error
import urllib.request

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from src import io

CHUNK = 1 << 20


def record_id(doi: str) -> str:
    """The numeric record id from a Zenodo DOI, DOI URL or bare id."""
    text = str(doi).strip().rstrip("/")
    if "zenodo." in text:
        return text.rsplit("zenodo.", 1)[1]
    if text.isdigit():
        return text
    raise SystemExit(
        f"cannot read a Zenodo record id out of {doi!r}\n"
        "expected something like 10.5281/zenodo.1234567"
    )


def manifest(api: str, doi: str) -> dict[str, dict]:
    """Map filename -> file entry, from the deposit's own metadata."""
    url = f"{api.rstrip('/')}/{record_id(doi)}"
    print(f"  reading {url}")
    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            meta = json.load(response)
    except urllib.error.HTTPError as exc:
        raise SystemExit(f"the deposit returned HTTP {exc.code} for {url}\n"
                         "check deposit.doi in data/parameters.yaml")
    except urllib.error.URLError as exc:
        raise SystemExit(f"cannot reach {url}: {exc.reason}")
    entries = {}
    for entry in meta.get("files", []):
        name = entry.get("key") or entry.get("filename")
        if name:
            entries[name] = entry
    if not entries:
        raise SystemExit(f"the deposit at {url} lists no files")
    return entries


def digest(path: str, algorithm: str) -> str:
    h = hashlib.new(algorithm)
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(CHUNK), b""):
            h.update(block)
    return h.hexdigest()


def expected(entry: dict) -> tuple[str, str] | None:
    """(algorithm, hex) from a checksum field like 'md5:abc...', if present."""
    value = entry.get("checksum")
    if not value:
        return None
    algorithm, _, hexdigest = str(value).partition(":")
    if not hexdigest:
        return "md5", algorithm
    return algorithm, hexdigest


def verified(path: str, entry: dict) -> bool:
    """True when the local file matches the deposit. Size alone if no checksum."""
    if not os.path.exists(path):
        return False
    size = entry.get("size")
    if size is not None and os.path.getsize(path) != int(size):
        return False
    want = expected(entry)
    if want is None:
        return True
    algorithm, hexdigest = want
    try:
        return digest(path, algorithm) == hexdigest
    except ValueError:
        print(f"  cannot check {algorithm} digests with this Python; trusting size")
        return True


def source_url(entry: dict) -> str:
    links = entry.get("links") or {}
    url = links.get("self") or links.get("download") or links.get("content")
    if not url:
        raise SystemExit(f"the deposit gives no download link for {entry.get('key')}")
    return url


def download(entry: dict, path: str) -> None:
    url = source_url(entry)
    total = int(entry.get("size") or 0)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    handle, staged = tempfile.mkstemp(dir=os.path.dirname(path), suffix=".part")
    os.close(handle)
    try:
        with urllib.request.urlopen(url, timeout=120) as response, open(staged, "wb") as out:
            done = 0
            for block in iter(lambda: response.read(CHUNK), b""):
                out.write(block)
                done += len(block)
                if total:
                    print(f"\r    {os.path.basename(path)}  "
                          f"{done / 1e6:6.1f} / {total / 1e6:.1f} MB", end="")
            print()
        if not verified(staged, entry):
            raise SystemExit(f"{os.path.basename(path)} failed its checksum; "
                             "the download is incomplete or the deposit changed")
        shutil.move(staged, path)
    finally:
        if os.path.exists(staged):
            os.remove(staged)


def collect(names: list[str], destination: str, files: dict[str, dict]) -> int:
    fetched = 0
    for name in names:
        entry = files.get(name)
        if entry is None:
            print(f"  MISSING  {name} is not in the deposit")
            print(f"           it lists: {', '.join(sorted(files)) or '(nothing)'}")
            raise SystemExit("the deposit does not hold the expected filenames")
        path = os.path.join(destination, name)
        if verified(path, entry):
            print(f"  have     {name}")
            continue
        download(entry, path)
        fetched += 1
    return fetched


ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--frames", action="store_true",
                help="also fetch the ensemble frames (only needed to re-run the detectors)")
ap.add_argument("--doi", default=None,
                help="override deposit.doi from data/parameters.yaml")
A = ap.parse_args()

p = io.parameters()
D = p.get("deposit") or {}
doi = A.doi or D.get("doi")
if not doi:
    sys.exit("the deposit DOI is not recorded.\n"
             "set deposit.doi in data/parameters.yaml, or pass --doi, to the DOI in\n"
             "the Data Availability statement of the paper.")

io.banner("00  fetching the deposited files")
files = manifest(D.get("api", "https://zenodo.org/api/records"), doi)

clouds = os.path.join(io.DATA, "pointclouds")
wanted = [f"pc_GLUT1_{tool}.npz" for tool in p["tools"]]
n = collect(wanted, clouds, files)
print(f"  point clouds in {clouds}  ({n} downloaded, {len(wanted) - n} already present)")

if A.frames:
    archive = D.get("frames_archive")
    if not archive:
        sys.exit("set deposit.frames_archive in data/parameters.yaml to the name of the\n"
                 "frames archive in the deposit before using --frames.")
    collect([archive], io.DATA, files)
    print(f"  frames archive in {io.DATA}; unpack it where tools/README.md expects it")

print("\n  run/04_iou.py can now run.")

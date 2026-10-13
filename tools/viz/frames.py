"""Local cache of the ensemble frames, filled on demand from the deposit.

The frames are deposited as a single 60 MB archive, so fetching one frame by
downloading the archive would be wasteful for a viewer that typically looks at a
handful of them. Zenodo serves byte ranges, so this module reads the archive's
central directory once and then pulls individual members, about 0.6 MB each,
leaving the rest on the server. If the server turns out not to support ranges the
whole archive is downloaded once instead and used the same way.

Cached frames live in data/frames/ and are not tracked by git.

    from tools.viz import frames
    path = frames.ensure(42)        # downloads only if absent
"""
from __future__ import annotations

import json
import os
import struct
import sys
import urllib.error
import urllib.request
import zlib

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from src import io

RECORD = "21238931"                       # doi:10.5281/zenodo.21238931
ARCHIVE = "frames-GLUT1_MD.zip"
MEMBER = "frames/frame.%d.pdb"            # layout inside the archive
API = "https://zenodo.org/api/records"
CACHE = os.path.join(io.DATA, "frames")
INDEX = os.path.join(CACHE, ".archive_index.json")
LOCAL_ARCHIVE = os.path.join(CACHE, ARCHIVE)


def frame_path(frame: int) -> str:
    return os.path.join(CACHE, "frame.%d.pdb" % frame)


def _get(url: str, start: int | None = None, end: int | None = None) -> bytes:
    headers = {} if start is None else {"Range": "bytes=%d-%d" % (start, end)}
    request = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(request, timeout=120) as response:
        if start is not None and response.status != 206:
            raise OSError("the server ignored the range request")
        return response.read()


def _archive_url() -> tuple[str, int]:
    """(url, size) of the deposited archive, from the record metadata."""
    with urllib.request.urlopen("%s/%s" % (API, RECORD), timeout=60) as response:
        record = json.load(response)
    for entry in record.get("files", []):
        if (entry.get("key") or entry.get("filename")) == ARCHIVE:
            return (entry.get("links") or {})["self"], int(entry["size"])
    raise SystemExit("%s is not in the deposit record %s" % (ARCHIVE, RECORD))


def _parse_central_directory(blob: bytes) -> dict[str, list[int]]:
    """name -> [local header offset, compressed size, method] for each member."""
    entries, at = {}, 0
    while at + 46 <= len(blob) and blob[at:at + 4] == b"PK\x01\x02":
        method, = struct.unpack("<H", blob[at + 10:at + 12])
        csize, = struct.unpack("<I", blob[at + 20:at + 24])
        name_len, extra_len, comment_len = struct.unpack("<HHH", blob[at + 28:at + 34])
        offset, = struct.unpack("<I", blob[at + 42:at + 46])
        name = blob[at + 46:at + 46 + name_len].decode("utf-8", "replace")
        entries[name] = [offset, csize, method]
        at += 46 + name_len + extra_len + comment_len
    return entries


def _index() -> dict:
    """The archive's member table, read once and cached next to the frames."""
    if os.path.exists(INDEX):
        with open(INDEX) as fh:
            return json.load(fh)
    os.makedirs(CACHE, exist_ok=True)
    url, size = _archive_url()
    tail = _get(url, max(0, size - 131072), size - 1)
    at = tail.rfind(b"PK\x05\x06")
    if at < 0:
        raise SystemExit("cannot read the archive's directory; it may be zip64")
    cd_size, = struct.unpack("<I", tail[at + 12:at + 16])
    cd_offset, = struct.unpack("<I", tail[at + 16:at + 20])
    entries = _parse_central_directory(_get(url, cd_offset, cd_offset + cd_size - 1))
    index = {"url": url, "size": size, "entries": entries}
    with open(INDEX, "w") as fh:
        json.dump(index, fh)
    print("  read the archive directory: %d members" % len(entries))
    return index


def _member_bytes(index: dict, name: str) -> bytes:
    offset, csize, method = index["entries"][name]
    head = _get(index["url"], offset, offset + 29)
    name_len, extra_len = struct.unpack("<HH", head[26:30])
    start = offset + 30 + name_len + extra_len
    raw = _get(index["url"], start, start + csize - 1)
    if method == 0:
        return raw
    if method == 8:
        return zlib.decompress(raw, -zlib.MAX_WBITS)
    raise SystemExit("unsupported compression method %d for %s" % (method, name))


def _from_local_archive(name: str) -> bytes:
    import zipfile
    with zipfile.ZipFile(LOCAL_ARCHIVE) as zf:
        return zf.read(name)


def _download_archive(index: dict) -> None:
    print("  falling back to the whole archive (%.0f MB), once" % (index["size"] / 1e6))
    staged = LOCAL_ARCHIVE + ".part"
    with urllib.request.urlopen(index["url"], timeout=600) as response, open(staged, "wb") as out:
        done = 0
        for block in iter(lambda: response.read(1 << 20), b""):
            out.write(block)
            done += len(block)
            print("\r    %.0f / %.0f MB" % (done / 1e6, index["size"] / 1e6), end="")
    print()
    os.replace(staged, LOCAL_ARCHIVE)


def ensure(frame: int, quiet: bool = False) -> str:
    """Path to frame `frame`, fetching it from the deposit if not cached."""
    path = frame_path(frame)
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    os.makedirs(CACHE, exist_ok=True)
    index = _index()
    name = MEMBER % frame
    if name not in index["entries"]:
        raise SystemExit("frame %d is not in the archive (it holds %d frames)"
                         % (frame, len(index["entries"])))
    if not quiet:
        print("  fetching frame %d" % frame)
    if os.path.exists(LOCAL_ARCHIVE):
        data = _from_local_archive(name)
    else:
        try:
            data = _member_bytes(index, name)
        except (OSError, urllib.error.HTTPError):
            _download_archive(index)
            data = _from_local_archive(name)
    staged = path + ".part"
    with open(staged, "wb") as fh:
        fh.write(data)
    os.replace(staged, path)
    return path


def cached() -> list[int]:
    if not os.path.isdir(CACHE):
        return []
    out = []
    for name in os.listdir(CACHE):
        if name.startswith("frame.") and name.endswith(".pdb"):
            out.append(int(name[6:-4]))
    return sorted(out)

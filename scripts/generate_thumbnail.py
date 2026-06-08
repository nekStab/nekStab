#!/usr/bin/env python3
"""Extract a thumbnail PNG from an mp4 using ffmpeg.

Usage: scripts/generate_thumbnail.py <input.mp4> [--out thumb.png] [--frame 0]

Exit codes:
  0  success
  1  ffmpeg missing or invocation failed
  2  input file missing or invalid argument
"""
import argparse, shutil, subprocess, sys
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Extract thumbnail from mp4 via ffmpeg")
    parser.add_argument("input", help="path to input .mp4")
    parser.add_argument("--out", default=None, help="output PNG (default: input with .png)")
    parser.add_argument("--frame", type=float, default=0.0, help="timestamp in seconds (default 0)")
    parser.add_argument("--width", type=int, default=480, help="max width (default 480)")
    args = parser.parse_args()

    src = Path(args.input)
    if not src.is_file():
        print(f"error: input not found: {src}", file=sys.stderr)
        return 2

    out = Path(args.out) if args.out else src.with_suffix(".png")

    if shutil.which("ffmpeg") is None:
        print("error: ffmpeg not installed", file=sys.stderr)
        return 1

    cmd = ["ffmpeg", "-y", "-ss", str(args.frame), "-i", str(src),
           "-vframes", "1", "-vf", f"scale={args.width}:-2", str(out)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"ffmpeg failed: {res.stderr}", file=sys.stderr)
        return 1
    print(f"wrote {out}")
    return 0

if __name__ == "__main__":
    sys.exit(main())

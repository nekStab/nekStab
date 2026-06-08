#!/usr/bin/env bash
# encode_anim.sh — canonical nekStab mode-animation encoder.
#
# Usage:
#   scripts/encode_anim.sh <frames_dir> <output.mp4>
#
# Expects PNG frames named frame_NNNNN.png in <frames_dir>.
# Produces an mp4 under 5 MB for typical mode animations (24 fps, 720p).
# NOTE: ffmpeg must be installed separately (not bundled with nekStab).
set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <frames_dir> <output.mp4>" >&2
    exit 1
fi

frames_dir="$1"
output="$2"

if [ ! -d "$frames_dir" ]; then
    echo "Error: frames directory '$frames_dir' not found." >&2
    exit 1
fi

if ! command -v ffmpeg >/dev/null 2>&1; then
    echo "Error: ffmpeg is not installed. Install it before running this script." >&2
    exit 1
fi

ffmpeg -y -framerate 24 -i "$frames_dir/frame_%05d.png" \
       -c:v libx264 -preset slow -crf 30 \
       -vf "scale=720:-2,format=yuv420p" \
       -movflags +faststart \
       "$output"

MAX_SIZE=5242880
actual_size=$(stat -c '%s' "$output" 2>/dev/null || stat -f '%z' "$output")
echo "Output: $output ($actual_size bytes)"

if [ "$actual_size" -gt "$MAX_SIZE" ]; then
    echo "WARNING: $output exceeds 5 MB ($actual_size bytes). Consider increasing CRF or reducing resolution." >&2
fi

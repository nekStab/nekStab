"""Make the in-repo `nekstab` package importable when running pytest from the
repo root without installing it (PYTHONPATH-free convenience)."""
import sys
from pathlib import Path

# python/ is the parent of this tests/ dir; put it on sys.path so `import nekstab`
# resolves to the source tree.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

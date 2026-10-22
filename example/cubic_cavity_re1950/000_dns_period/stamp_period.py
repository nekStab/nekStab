#!/usr/bin/env python3
"""Copy a Nek5000 field file and write the orbit period T into its time stamp.

Usage: stamp_period.py <field file> <T> <output file>

nekStab reads the Floquet period from the time stamp of the base-flow file.
Only the 132-byte text header is changed; the field data are copied as they are.
"""
import re
import shutil
import sys

src, period, dst = sys.argv[1], float(sys.argv[2]), sys.argv[3]
shutil.copyfile(src, dst)
with open(dst, 'r+b') as fh:
    header = fh.read(132).decode('ascii')
    if not header.startswith('#std'):
        sys.exit('not a Nek5000 field file header')
    old = re.search(r'[0-9]\.[0-9]{13}E[+-][0-9]{2}', header)
    mant, exp = f'{period:.12E}'.split('E')
    new = f'{float(mant) / 10:.13f}E{int(exp) + 1:+03d}'  # Nek writes 0.dddE+xx
    assert len(new) == len(old.group()), (new, old.group())
    fh.seek(0)
    fh.write(header.replace(old.group(), new).encode('ascii'))
print(f'{dst}: time stamp {old.group()} -> {new}')

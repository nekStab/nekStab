#!/bin/bash
for f in *.f
do
	echo "Indenting ${f}"
	emacs -batch ${f} --eval '(indent-region (point-min) (point-max) nil)' -f save-buffer 2> /dev/null
done

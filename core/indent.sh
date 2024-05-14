#!/bin/bash

# Not working as we need to respect fixed format style due to legacy nek5000 code
# for ext in f90 f; do
#   for f in *."$ext"; do
#     [ -f "$f" ] || continue
#     echo "Indenting $f"
#     emacs -batch "$f" --eval '(indent-region (point-min) (point-max) nil)' -f save-buffer 2> /dev/null
#   done
# done

#conda install fprettify
#pip install fprettify

# Directory to process
dir="./"

# Find all .f90 files in the directory
find "$dir" -type f -name "*.f90" | while read file; do
    echo "Processing file: $file"

    #echo "Converting tabs to spaces..."
    expand -t 4 "$file" >"$file.expanded"
    mv "$file.expanded" "$file"

    #echo "Removing trailing spaces..."
    sed -i 's/[ \t]*$//' "$file"

    #echo "Removing leading whitespaces..."
    sed -i 's/^[ \t]*//' "$file"

    #echo "Running fprettify..."
    fprettify "$file" --case 1 1 1 1 --enable-decl --enable-replacements --c-relations

    # Fixed Form Fortran, which was used in older versions of Fortran (Fortran 77 and earlier)
    # the first 6 columns of each line have special meaning:
    # - Column 1: Used for a comment if the line starts with a `C` or `*`.
    # - Columns 2-5: Used for statement labels (for `GOTO` statements, etc.).
    # - Column 6: Used for line continuation if the line starts with any character.

    #echo "Adding 6 leading spaces to all lines..."
    sed -i 's/^/      /' "$file"

    #echo "Fixing the number of spaces before $ to 5 for lines starting with $..."
    awk '{sub(/[ \t]+\$/, "     $   "); print}' "$file" >temp && mv temp "$file"

    echo "Done processing file: $file"
    echo ""
done

echo "Deleting backup files..."
find "$dir" -type f -name "*.f90~" -delete

echo "All done."

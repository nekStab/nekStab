#!/bin/bash

# Determine the OS type for compatible sed in-place editing
if [[ "$(uname)" == "Darwin" ]]; then
    SED_INPLACE=('sed' '-i' '')
else
    SED_INPLACE=('sed' '-i')
fi

# Function to check if file is Fortran source
is_fortran_file() {
    local file="$1"
    # Check if file has .f90 extension or is a known include file
    if [[ "$file" == *.f90 ]] || \
       [[ "$file" == *"NEKSTAB"* ]]; then
        return 0
    fi
    return 1
}

# Function to process a single file
process_file() {
    local file="$1"
    echo "Processing file: $file"

    # Convert tabs to spaces
    expand -t 4 "$file" > "${file}.expanded"
    mv "${file}.expanded" "$file"

    # Remove trailing spaces and tabs
    "${SED_INPLACE[@]}" 's/[[:space:]]*$//' "$file"

    # Remove leading spaces and tabs
    "${SED_INPLACE[@]}" 's/^[[:space:]]*//' "$file"

    # Run fprettify
    fprettify "$file" --case 1 1 1 1 --enable-decl --enable-replacements --c-relations

    # Add 6 leading spaces to all lines
    "${SED_INPLACE[@]}" 's/^/      /' "$file"

    # Fix the number of spaces before $ to 5 for lines starting with $
    awk '{sub(/[[:space:]]+\$/, "     $   "); print}' "$file" > temp && mv temp "$file"

    echo "Done processing file: $file"
    echo ""
}

# Check if a filename is provided as an argument
if [[ $# -eq 1 ]]; then
    FILE="$1"
    if [[ -f "$FILE" ]] && is_fortran_file "$FILE"; then
        process_file "$FILE"
    else
        echo "Error: Provided file does not exist or is not a Fortran file."
        exit 1
    fi
else
    # Find and process all Fortran files in the current directory and subdirectories
    find . -type f | while read -r file; do
        if is_fortran_file "$file"; then
            process_file "$file"
        fi
    done
fi

echo "Deleting backup files..."
find . -type f -name "*~" -delete

echo "All done."
#!/usr/bin/env bash

_mks_has_case_file () {
    local stem=$1
    local ext
    for ext in re2 rea par usr sep map; do
        [[ -f "${stem}.${ext}" ]] && return 0
    done
    return 1
}

_mks_suggest_cases () {
    local -a files
    local f
    local stem

    {
    if [[ -f SESSION.INFO ]]; then
        awk 'NF && $1 !~ /^\// && $1 ~ /^[A-Za-z0-9._-]+$/ {print $1}' SESSION.INFO
    elif [[ -f SESSION.NAME ]]; then
        read -r stem < SESSION.NAME
        [[ "$stem" =~ ^[A-Za-z0-9._-]+$ ]] && printf '%s\n' "$stem"
    fi

    shopt -s nullglob
    files=(
        *.re2
        *.rea
        *.par
        *.usr
        *.sep
        *.map
    )
    shopt -u nullglob

    for f in "${files[@]}"; do
        stem=${f%.*}
        [[ -z "$stem" ]] && continue

        # Avoid obvious generated/backup variants.
        case "$stem" in
            *_old|*~|*/*)
                continue
                ;;
        esac

        printf '%s\n' "$stem"
    done
    } | sort -u
}

_mks_completion () {
    local cur token
    local -a candidates
    local i
    local has_case_files=false
    local first_pos=""
    local check

    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"

    for (( i=1; i<COMP_CWORD; i++ )); do
        token="${COMP_WORDS[i]}"
        if [[ "$token" != --* ]]; then
            first_pos="$token"
            break
        fi
    done

    if [[ -z "$first_pos" ]]; then
        if [[ "$cur" == -* ]]; then
            candidates=(--debug --fresh --help -h)
        else
            candidates=(clean full_clean)
            for check in *.re2 *.rea *.par *.usr *.sep *.map; do
                if [[ -e "$check" ]]; then
                    has_case_files=true
                    break
                fi
            done
            while IFS= read -r token; do
                [[ -n "$token" ]] || continue
                if [ "$has_case_files" = false ] || _mks_has_case_file "$token"; then
                    candidates+=("$token")
                fi
            done < <(_mks_suggest_cases)
        fi
    else
        if [[ "$cur" == -* ]]; then
            candidates=(--debug --fresh --help -h)
        else
            return 0
        fi
    fi

    COMPREPLY=( $(compgen -W "${candidates[*]}" -- "$cur") )

    if (( ${#COMPREPLY[@]} == 0 )) && [[ "$cur" == f* ]]; then
        COMPREPLY=(full_clean)
    fi
}

if [[ -n "$PS1" ]]; then
    complete -F _mks_completion mks makeneks
fi

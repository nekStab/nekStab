# Set the nekStab environment in the current terminal.
#
# Use (bash or zsh), from any directory:
#     source ~/nekStab/exports.sh
# To have it in every new terminal, add that same line once to ~/.bashrc
# (or ~/.zshrc). No other file needs editing.
#
# Running it as a program (bash exports.sh) does nothing useful: the variables
# would be set in a child shell that ends at once. The script checks for that.

# Directory of this file, also when it is sourced from another directory.
if [ -n "${BASH_SOURCE:-}" ]; then
    _nekstab_file=${BASH_SOURCE[0]}
elif [ -n "${ZSH_VERSION:-}" ]; then
    _nekstab_file=${(%):-%x}
else
    _nekstab_file=$0
fi
_nekstab_dir=$(cd "$(dirname "$_nekstab_file")" && pwd)

# Executed, not sourced? bash: the file is $0. zsh: the evaluation context lacks "file".
# Other shells (sh, dash): $0 is the file name when executed, the shell name when sourced.
_nekstab_executed=
if [ -n "${BASH_SOURCE:-}" ]; then
    [ "${BASH_SOURCE[0]}" = "$0" ] && _nekstab_executed=1
elif [ -n "${ZSH_VERSION:-}" ]; then
    case $ZSH_EVAL_CONTEXT in *:file*) ;; *) _nekstab_executed=1 ;; esac
else
    case ${0##*/} in exports.sh) _nekstab_executed=1 ;; esac
fi
if [ -n "$_nekstab_executed" ]; then
    echo "exports.sh: run it with 'source $_nekstab_dir/exports.sh', not 'bash exports.sh'."
    exit 1
fi

export NEKSTAB_SOURCE_ROOT=$_nekstab_dir
# Keep a Nek5000 that the user set before (an existing installation).
export NEK_SOURCE_ROOT=${NEK_SOURCE_ROOT:-$NEKSTAB_SOURCE_ROOT/Nek5000}
export PATH=$NEKSTAB_SOURCE_ROOT/bin:$NEK_SOURCE_ROOT/bin:$PATH

# Large cases need a large stack; ignore the error where the limit is fixed.
ulimit -s unlimited 2>/dev/null || true
# The README says the stack limit is removed: say so when it is not.
[ "$(ulimit -s)" = unlimited ] || echo "nekStab: stack limit is $(ulimit -s) kB, not unlimited; large cases may crash." >&2

# Tab completion of mks (bash only).
if [ -n "${BASH_VERSION:-}" ] && [ -f "$NEKSTAB_SOURCE_ROOT/bin/mks-completion.bash" ]; then
    . "$NEKSTAB_SOURCE_ROOT/bin/mks-completion.bash"
fi

echo "nekStab: NEKSTAB_SOURCE_ROOT=$NEKSTAB_SOURCE_ROOT"
echo "nekStab: NEK_SOURCE_ROOT=$NEK_SOURCE_ROOT"
unset _nekstab_file _nekstab_dir _nekstab_executed

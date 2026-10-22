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

if [ -n "${BASH_SOURCE:-}" ] && [ "${BASH_SOURCE[0]}" = "$0" ]; then
    echo "exports.sh: run it with 'source $_nekstab_dir/exports.sh', not 'bash exports.sh'."
    exit 1
fi

export NEKSTAB_SOURCE_ROOT=$_nekstab_dir
# Keep a Nek5000 that the user set before (an existing installation).
export NEK_SOURCE_ROOT=${NEK_SOURCE_ROOT:-$NEKSTAB_SOURCE_ROOT/Nek5000}
export PATH=$NEKSTAB_SOURCE_ROOT/bin:$NEK_SOURCE_ROOT/bin:$PATH

# Large cases need a large stack; ignore the error where the limit is fixed.
ulimit -s unlimited 2>/dev/null || true

# Tab completion of mks (bash only).
if [ -n "${BASH_VERSION:-}" ] && [ -f "$NEKSTAB_SOURCE_ROOT/bin/mks-completion.bash" ]; then
    . "$NEKSTAB_SOURCE_ROOT/bin/mks-completion.bash"
fi

echo "nekStab: NEKSTAB_SOURCE_ROOT=$NEKSTAB_SOURCE_ROOT"
echo "nekStab: NEK_SOURCE_ROOT=$NEK_SOURCE_ROOT"
unset _nekstab_file _nekstab_dir

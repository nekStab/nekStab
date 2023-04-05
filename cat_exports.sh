
#!/bin/bash

# Get the current working directory
nekstab_source_root=$(pwd)

# Define the exports as a single string
exports_string="
export NEKSTAB_SOURCE_ROOT=$nekstab_source_root
export NEK_SOURCE_ROOT=\"\$NEKSTAB_SOURCE_ROOT/Nek5000\"
export PATH=\$NEK_SOURCE_ROOT/bin:\$PATH
export PATH=\$NEKSTAB_SOURCE_ROOT/bin:\$PATH
"

# Determine the appropriate configuration file based on the current shell
if [ "$SHELL" = "/bin/bash" ]; then
    config_file=~/.bashrc
elif [ "$SHELL" = "/bin/zsh" ]; then
    config_file=~/.zshrc
elif [ "$SHELL" = "/bin/fish" ]; then
    config_file=~/.config/fish/config.fish
else
    echo "Unsupported shell: $SHELL"
    exit 1
fi

# Append the exports to the configuration file
echo "$exports_string" >> "$config_file"
echo "Exports added to $config_file"
echo " run: source $config_file"
echo " or restart the shell terminal window."

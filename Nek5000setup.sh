#!/bin/bash

# Welcome message
echo "Nek5000 Setup Script"
echo "--------------------"
echo "This script will perform the following actions upon your confirmation:"
echo "1. Remove existing Nek5000 directory (if found)."
echo "2. Install necessary dependencies."
echo "3. Clone the Nek5000 repository."
echo "4. Modify a line in core/prepost.f."
echo "5. Build genmap and genbox tools."
echo "6. Add environment variables to your shell configuration."
echo "--------------------"

# Detect the operating system
OS=$(uname)
should_clone="no"

# Check for existing Nek5000 directory
if [ -d "Nek5000" ]; then
    read -p "Directory Nek5000 exists. Do you want to remove it? (y/n): " confirm
    if [ "$confirm" == "y" ] || [ "$confirm" == "Y" ]; then
        echo "Removing Nek5000 directory..."
        rm -rf Nek5000
        should_clone="yes"
    else
        echo "Nek5000 directory will be retained."
    fi
else
    should_clone="yes"
fi

# Install dependencies
read -p "Do you want to install/update necessary packages? (y/n): " confirm
if [ "$confirm" == "y" ] || [ "$confirm" == "Y" ]; then
    echo "Installing dependencies..."
    if [ "$OS" == "Linux" ]; then
        sudo apt -y update
        sudo apt install build-essential gfortran libopenmpi-dev cmake libx11-dev libxt-dev
    elif [ "$OS" == "Darwin" ]; then
        brew update
        brew install gcc open-mpi cmake libx11 libxt
    else
        echo "Unsupported operating system. Exiting."
        exit 1
    fi
else
    echo "Skipping package installation."
fi

# Clone Nek5000 repository
if [ "$should_clone" == "yes" ]; then
    read -p "Do you want to clone the Nek5000 repository? (y/n): " confirm
    if [ "$confirm" == "y" ] || [ "$confirm" == "Y" ]; then
        echo "Cloning Nek5000 repository..."
        git clone https://github.com/Nek5000/Nek5000.git
        cd Nek5000
        git checkout master
    else
        echo "Skipping cloning."
    fi
else
    cd Nek5000
fi

# Modify core/prepost.f
read -p "Do you want to modify a line in core/prepost.f? (y/n): " confirm
if [ "$confirm" == "y" ] || [ "$confirm" == "Y" ]; then
    echo "Modifying core/prepost.f..."
    FILE_PATH="core/prepost.f"
    cp "$FILE_PATH" "${FILE_PATH}.original"
    if [ "$OS" == "Linux" ]; then
        sed -i 's/save    nopen/common \/RES_WANT\/ nopen/g' "$FILE_PATH"
    elif [ "$OS" == "Darwin" ]; then
        sed -i '' 's/save    nopen/common \/RES_WANT\/ nopen/g' "$FILE_PATH"
    fi
    if cmp -s "$FILE_PATH" "${FILE_PATH}.original"; then
        echo "No match found in $FILE_PATH. No replacement made."
    else
        echo "Line in $FILE_PATH successfully replaced."
    fi
    rm "${FILE_PATH}.original"
else
    echo "Skipping modification."
fi

# Build genmap and genbox tools
read -p "Do you want to build genmap and genbox tools? (y/n): " confirm
if [ "$confirm" == "y" ] || [ "$confirm" == "Y" ]; then
    echo "Building genmap and genbox tools..."
    cd tools
    ./maketools genmap genbox
    cd ../.. # return to nekStab root directory
else
    echo "Skipping tools building."
fi

# Add exports to shell configuration file
read -p "Do you want to add necessary environment variables to your shell configuration file? (y/n): " confirm
if [ "$confirm" == "y" ] || [ "$confirm" == "Y" ]; then
    # Get the current working directory
    nekstab_source_root=$(pwd)

    # Define the exports as a single string without leading spaces
    exports_string="# nekStab folder location"
    exports_string+="\nexport NEKSTAB_SOURCE_ROOT=$nekstab_source_root"
    exports_string+="\nexport PATH=\$NEKSTAB_SOURCE_ROOT/bin:\$PATH"
    exports_string+="\n# Nek5000 folder location"
    exports_string+="\nexport NEK_SOURCE_ROOT=\"\$NEKSTAB_SOURCE_ROOT/Nek5000\""
    exports_string+="\nexport PATH=\$NEK_SOURCE_ROOT/bin:\$PATH"

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
    printf "%b" "$exports_string" >> "$config_file"
    echo "Exports added to $config_file"
    echo " run: source $config_file"
    echo " or restart the shell terminal window."
else
    echo "Skipping export of environment variables."
fi

echo "Nek5000 setup complete."
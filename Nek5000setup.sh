#!/bin/bash

# Welcome message
echo "Nek5000 Setup Script"
echo "--------------------"
echo "This script will perform the following actions upon your confirmation:"
echo "1. Remove existing Nek5000 directory (if found)."
echo "2. Install necessary dependencies."
echo "3. Clone the Nek5000 repository."
echo "4. Build genmap and genbox tools."
echo "5. Add environment variables to your shell configuration."
echo "--------------------"

# Detect the operating system
OS=$(uname)
should_clone="no"

# Check for existing Nek5000 directory
if [ -d "Nek5000" ]; then
    read -p "Directory Nek5000 exists. Do you want to remove it? (y/n) [y]: " confirm
    confirm=${confirm:-y}
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
read -p "Do you want to install/update necessary packages? (y/n) [n]: " confirm
confirm=${confirm:-n}
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
    read -p "Do you want to clone the Nek5000 repository? (y/n) [y]: " confirm
    confirm=${confirm:-y}
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

# Build genmap and genbox tools
read -p "Do you want to build genmap and genbox tools? (y/n) [y]: " confirm
confirm=${confirm:-y}
if [ "$confirm" == "y" ] || [ "$confirm" == "Y" ]; then
    echo "Building genmap and genbox tools..."
    cd tools
    ./maketools genmap genbox
    cd ../.. # return to nekStab root directory
else
    echo "Skipping tools building."
fi

# Add exports to shell configuration file
read -p "Do you want to add necessary environment variables to your shell configuration file? (y/n) [n]: " confirm
confirm=${confirm:-n}
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
    else
        echo "Unsupported shell: $SHELL (only bash and zsh are supported)"
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
#!/bin/bash

echo "Sourcing Nek5000 and nekStab paths:"

NEKSTAB_SOURCE_ROOT="../.." #--> path to nekStab
NEK_SOURCE_ROOT="$NEKSTAB_SOURCE_ROOT/Nek5000"

echo "Setting NEKSTAB_SOURCE_ROOT to $NEKSTAB_SOURCE_ROOT"
echo "Setting NEK_SOURCE_ROOT to $NEK_SOURCE_ROOT"

export PATH=$NEK_SOURCE_ROOT/bin:$PATH
export PATH=$NEKSTAB_SOURCE_ROOT/bin:$PATH

echo "Adding the following paths to PATH:"
echo "$NEK_SOURCE_ROOT/bin"
echo "$NEKSTAB_SOURCE_ROOT/bin"

echo "Sourcing complete."

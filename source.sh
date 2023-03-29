#!/bin/bash
export NEKSTAB_SOURCE_ROOT=$(pwd) #--> path to nekStab
export NEK_SOURCE_ROOT="$NEKSTAB_SOURCE_ROOT/Nek5000"
export PATH=$NEK_SOURCE_ROOT/bin:$PATH
export PATH=$NEKSTAB_SOURCE_ROOT/bin:$PATH

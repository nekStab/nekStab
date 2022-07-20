#!/bin/bash
set -e

export CASE="1cyl" #--> case name goes here

#add to .bashrc
export NEKSTAB_SOURCE_ROOT="../../.." #--> path to nekStab
export NEK_SOURCE_ROOT="$NEKSTAB_SOURCE_ROOT/Nek5000"
export PATH=$NEK_SOURCE_ROOT/bin:$PATH

echo $NEKSTAB_SOURCE_ROOT
echo $NEK_SOURCE_ROOT

#uncomment the deisred compiler
source ${NEKSTAB_SOURCE_ROOT}/core/compiler.sh

if [ $#  -ne 1 ] ; then
    error_quit 'Wrong arguments number!'
fi
for il in "$@"
do
case $il in
        clean)
                ${SOURCE_ROOT}/bin/makenek clean
                shift
                ;;
        all)
                ${SOURCE_ROOT}/bin/makenek ${CASE}
                shift
                ;;
esac
done

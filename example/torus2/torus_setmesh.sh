#!/bin/bash 
GDIR="./geom/"
CNAME=torus

if [ $# -ne 1 ]; then echo -e "\nWrong argument. Abort.\n"; exit 1; fi
GNAME=$1
echo "Set links ... "
for suf in "re2" "ma2"; do
    oname="${GDIR}${GNAME}.${suf}"
    tname="${CNAME}.${suf}"
    if [ -e $oname ]; then
        if [ -e $tname ]; then
            echo "  remove $tname";
            rm $tname;
        fi
        echo "  $oname --> $tname";
        ln -s $oname $tname;
    else
        echo -e "\n  Mesh file $oname not found. Abort.\n"
        exit 1
    fi
done
echo "done."

if [ $# -ne 1 ]; then echo -e "\nWrong argument. Abort.\n"; exit 1; fi
#
#  Usage: bash mkmsh.sh torus_coarse2D
#
CNAME=$1
if [ ! -e ${CNAME}.geo ]; then echo -e "\nScript ${CNAME}.geo not found. Abort.\n"; exit 1; fi
TDIR="templates"
LOGFILE="mkmsh.log"
if [ -e $LOGFILE ]; then rm $LOGFILE; fi

echo "$CNAME:"
echo "  gmsh ..."
gmsh ${CNAME}.geo -2 -order 2 -o ${CNAME}.msh >> $LOGFILE

if grep -qi "Error" $LOGFILE; then echo -e "  Error! Check $LOGFILE for details. Abort \n"; exit 1; fi

echo "  gmsh2nek ..."
## gmsh2nek
sed "s/casename/${CNAME}/g" ${TDIR}/gmsh2nek_cmd_template.txt > gmsh2nek_cmd.txt;
gmsh2nek < gmsh2nek_cmd.txt >> $LOGFILE;
if grep -qi "Error" $LOGFILE; then echo -e "  Error! Check $LOGFILE for details. Abort \n"; exit 1; fi

echo "  re2torea ..."
## re2torea
cp -f torus_5.rea ${CNAME}.rea;
sed -i 's/3 D/2 D/g' ${CNAME}.rea;
nel=$(head -1 ${CNAME}.re2 | strings | head -1 | awk '{ print $2}');
sed -i "s/-3225  3        3225           NEL,NDIM/$nel 2        $nel           NEL,NDIM/g" ${CNAME}.rea;
sed "s/casename/${CNAME}/g" ${TDIR}/re2torea_cmd_template.txt > re2torea_cmd.txt;
re2torea < re2torea_cmd.txt >> $LOGFILE;
mv -f tmprea.rea ${CNAME}.rea;
if grep -qi "Error" $LOGFILE; then echo -e "  Error! Check $LOGFILE for details. Abort \n"; exit 1; fi

echo "  n2to3 ..."
## n2to3
sed "s/casename2D/${CNAME}/g" ${TDIR}/n2to3_cmd_template.txt > n2to3_cmd.txt;
sed -i "s/casename3D/${CNAME}_/g" n2to3_cmd.txt;
sed -i "s/2D_/3D/g" n2to3_cmd.txt;
n2to3 < n2to3_cmd.txt >> $LOGFILE;
if grep -qi "Error" $LOGFILE; then echo -e "  Error! Check $LOGFILE for details. Abort \n"; exit 1; fi

echo "  genmap ..."
## genmap
echo "${CNAME}_" > genmap_cmd.txt
echo "0.000001" >> genmap_cmd.txt
sed -i "s/2D_$/3D/g" genmap_cmd.txt
genmap < genmap_cmd.txt >> $LOGFILE;
if grep -qi "Error" $LOGFILE; then echo -e "  Error! Check $LOGFILE for details. Abort \n"; exit 1; fi

echo "done."
# cleanup
rm -f fort.99 gmsh2nek_cmd.txt n2to3_cmd.txt re2torea_cmd.txt genmap_cmd.txt

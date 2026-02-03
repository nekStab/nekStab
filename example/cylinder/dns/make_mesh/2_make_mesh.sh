if [ ! -f ~/.nekdefaults ]; then
    touch ~/.nekdefaults
fi

echo "-------------------- Running genbox --------------------"
genbox << EOF
mybox.box
EOF
echo "-------------------- Running pretex --------------------"
#
pretex << EOF
1cyl2
   1 READ PREVIOUS PARAMETERS 
box
   1 BUILD FROM FILE          
box
  10 IMPORT MESH              
import
y Would you like to displace existing elements in box?
   1 END    ELEMENTS          
   1 ACCEPT MATL,QVOL         
   1 SET BCs                  
  11 SET ENTIRE LEVEL         
   1 PERIODIC-AUTO            
n Is this a lattice for hex-close-packed spheres?:
1 END  LEVEL               
1 ACCEPT B.C.'s            
1 EXIT                     
EOF
rm pretex.jou
echo "-------------------- Running reatore2 --------------------"
echo "Latest file: $(ls -lhrt | tail -1)"
reatore2 << EOF
1cyl2
1cyl
EOF
rm 1cyl.rea 1cyl2.dra
echo "-------------------- Running genmap --------------------"
echo "Latest file: $(ls -t | head -n1)"
# Cross-platform script/genmap call
OS_TYPE=$(uname)
if [[ "$OS_TYPE" == "Darwin" ]]; then
    # macOS: run genmap with heredoc, capture output
    genmap > genmap.log <<EOF
1cyl
.1
EOF
else
    # Linux: script works as before
    script -q -c "genmap << EOF
1cyl
.1
EOF" genmap.log
fi
cat genmap.log
lelg_VALUE=$(grep "start rec_bisect:" genmap.log | awk '{print $3}')
ls -lhrt | tail -3
echo "Update SIZE lelg to $lelg_VALUE"
rm 1cyl.dra 1cyl2.rea genmap.log box.rea session.name fort.* 1cyl.rea 1cyl2.dra pretex.jou

# Create a temporary file for the operation
#awk -v val="$lelg_VALUE" '
#/parameter \(lelg=/{ 
#    printf "      parameter (lelg=%d)              ! max number of global elements\n", val
#    next
#}
#print }
#' SIZE > SIZE.tmp && mv SIZE.tmp SIZE
#echo "Updated SIZE file lelg line:"
#grep "lelg=" SIZE

mv 1cyl.ma2 ..
mv 1cyl.re2 ..

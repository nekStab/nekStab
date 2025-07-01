genbox << EOF
mybox.box
EOF
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
rm box.dra pretex.jou
echo "Latest file: $(ls -lhrt | tail -1)"
reatore2 << EOF
1cyl2
1cyl
EOF
rm 1cyl.rea 1cyl2.dra
echo "Latest file: $(ls -t | head -n1)"
script -q -c "genmap << EOF
1cyl
.1
EOF" genmap.log
BISECT_VALUE=$(grep "start rec_bisect:" genmap.log | awk '{print $3}')
rm 1cyl.dra 1cyl2.rea genmap.log box.rea session.name fort*
ls -lhrt | tail -3
echo "Update SIZE lelg to $BISECT_VALUE"

# Create a temporary file for the operation
#awk -v val="$BISECT_VALUE" '
#/parameter \(lelg=/ { 
#    printf "      parameter (lelg=%d)              ! max number of global elements\n", val
#    next
#}
#{ print }
#' SIZE > SIZE.tmp && mv SIZE.tmp SIZE
#echo "Updated SIZE file lelg line:"
#grep "lelg=" SIZE

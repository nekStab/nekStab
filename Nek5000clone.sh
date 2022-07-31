#sudo apt-get -y install libmpich-dev mpich libopenblas-dev build-essential cmake m4

git clone https://github.com/nekStab/Nek5000.git

cd Nek5000
#git checkout v19
git checkout master

cd tools
./maketools genmap genbox # n2to3 gmsh2nek nekamg_setup

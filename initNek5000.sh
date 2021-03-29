#sudo apt-get -y install libmpich-dev mpich libopenblas-dev build-essential cmake m4

git submodule update --init --recursive

cd Nek5000/tools

./maketools genmap genbox

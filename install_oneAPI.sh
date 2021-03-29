# following official Repo https://software.intel.com/content/www/us/en/develop/articles/installing-intel-oneapi-toolkits-via-apt.html

more /etc/lsb-release

sudo -E apt autoremove intel-hpckit intel-basekit

sudo apt install gnupg wget libigc-dev libigdfcl-dev libigfxcmrt-dev

cd $HOME
wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB

sudo add-apt-repository "deb https://apt.repos.intel.com/oneapi all main"

sudo apt update

sudo apt install intel-basekit intel-hpckit

sudo apt-cache pkgnames intel | grep kit$

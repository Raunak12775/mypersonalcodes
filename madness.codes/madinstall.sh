#!/bin/bash

echo "I hope you have read the readme.md before installing the script"
echo "before starting you should be logged in to node28 in rigi"
echo "to login do ssh rigi then ssh node28"

#if [`echo hostname` == "node28"]
#	then
#		continue
#else
#	echo "host is not node28"
#	break
#fi

#echo "before starting you should be logged in to node28 in rigi"


echo "creating a new directory and preparing everything for madness installation"
echo "libxc version used in this is 6.2.2"
mkdir madinstall
cd madinstall

mkdir -p install/libxc; mkdir install/madness 
mkdir -p build/libxc; mkdir build/madness
# clone madness
git clone https://github.com/m-a-d-n-e-s-s/madness.git

# download libxc, untar it
wget http://www.tddft.org/programs/libxc/down.php?file=6.2.2/libxc-6.2.2.tar.gz
mv *tar.gz libxc.tar.gz
tar -xvf libxc.tar.gz
rm *tar.gz

# building libxc :: something is wrong here check it tomorrow
cd build/libxc
cmake -DCMAKE_INSTALL_PREFIX=../../install/libxc/ ../../libxc-6.2.2/.
make -j 10 all
# installing libxc 
make -j 10 install
echo "libxc is installed at install/libxc"
echo "contents of the folder are:"
echo `ls ../../install/libxc`
cd ../../

# module loading
module load GNU-9
module load compiler
module load mkl
module load mpi

# building madness
cd build/madness
cmake -DENABLE_LIBXC=ON -DLIBXC_LIBRARIES=../../install/libxc/lib64/libxc.a -DLIBXC_INCLUDE_DIRS=../../install/libxc/include -DCMAKE_INSTALL_PREFIX=../../install/madness ../../madness

make -j 10 applications
# install madness
make -j 10 install

echo "installation complete!"
echo "madness is installed in install/madness"
echo `ls ../../install/madness`
cd ../../
echo "make your own module file and enjoy!"
echo "add this - export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:path/to/libxc/install/lib to your bashrc file"
echo "-- to prepend the path of LibXC libraries"
echo "or you can add this as a module file and load it in as madness environment"
echo "if you are unsure about the exporting stuff just uncomment and edit the path in the bottom line before running the script"
#echo LD_LIBRARY_PATH='/your/path/to/libxc/install/lib:LD_LIBRARY_PATH' >> ~/.bashrc && source ~/.bashrc


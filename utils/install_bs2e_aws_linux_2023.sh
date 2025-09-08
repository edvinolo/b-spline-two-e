#! /usr/bin/bash

# How many parallel jobs to use for make
N_jobs=1

# Specify the home_dir explicitly 
home_dir=/home/ec2-user

# CMake version variables
cmake_version=4.0.2
cmake_dir=cmake-$cmake_version
cmake_url=https://github.com/Kitware/CMake/releases/download/v$cmake_version/$cmake_dir.tar.gz

# Go to home directory
cd $home_dir

# Install GCC
sudo yum groupinstall "Development Tools"

# Install tmux
sudo yum install tmux

# Install htop
sudo yum install htop

# Install pip, and use it to install fypp
python3 -m ensurepip --upgrade
python3 -m pip install fypp
#python3 -m pip install cmake
#which cmake && cmake --version

# Install CMake
sudo yum install openssl-devel
wget $cmake_url
tar -xzvf $cmake_dir.tar.gz
rm -rf $cmake_dir.tar.gz
cd $cmake_dir
./bootstrap
make -j$N_jobs
sudo make install
cd $home_dir

# Install GSL
sudo yum install gsl-devel

# Install intel Fortran essentials
# Create repo file
tee > /tmp/oneAPI.repo << EOF
[oneAPI]
name=Intel® oneAPI repository
baseurl=https://yum.repos.intel.com/oneapi
enabled=1
gpgcheck=1
repo_gpgcheck=1
gpgkey=https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
EOF

# Move repo file
sudo mv /tmp/oneAPI.repo /etc/yum.repos.d

# Install
sudo yum install intel-fortran-essentials

# Set environment variables:
source /opt/intel/oneapi/setvars.sh
echo $MKLROOT

# Install Fortran stdlib
#export PATH=$PATH:/home/ec2-user/.local/bin:/home/ec2-user/bin:/home/ec2-user/.local/lib # Need to set some env variables to use fypp properly when calling this script as sudo
#export PYTHONPATH=/home/ec2-user/.local/lib/python3.9/site-packages 
cd $home_dir
git clone https://github.com/fortran-lang/stdlib.git
cd stdlib
cmake -B build -DCMAKE_Fortran_FLAGS="-O3 -march=native" -DCMAKE_BUILD_TYPE=NoConfig -DBLA_VENDOR=Intel10_64ilp -DWITH_ILP64=True
cmake --build build -j$N_jobs # Don't know if -j works here. Update: It does :)
# cmake --build build --target test # Some of them are currently failing, maybe due to ILP64?
sudo cmake --install build

# Install ARPACK
# First install openMPI (I will not be using PARPACK yet, so maybe don't need openMPI, but I don't think it hurts)
sudo yum install openmpi-devel
#sudo yum install hwloc-devel

# Update environment variables with paths to openMPI
echo 'export PATH=$PATH:/usr/lib64/openmpi/bin' >> $home_dir/.bashrc # This always adds these lines to .bashrc each time the script is run. Should maybe check if they are already present!
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64/openmpi/lib' >> $home_dir/.bashrc
source $home_dir/.bashrc

# Clone the ARPACK-NG repo:
cd $home_dir # Go back to home directory
git clone https://github.com/opencollab/arpack-ng.git
cd arpack-ng

# Build and Install ARPACK:
# to build with mkl blas/lapack might need to set BLA_VENDOR=Intel10_64ilp! But if MKL env variables set and no other BLAS/LAPACK installed it should find MKL I hope.
mkdir build
cd build
cmake -D EXAMPLES=ON -D BUILD_SHARED_LIBS=ON -DBLA_VENDOR=Intel10_64ilp -D INTERFACE64=ON -D ITF64SUFFIX="ILP64" ..
make -j$N_jobs
sudo make install

# Install b-spline-two-e
# Clone the repo:
cd $home_dir
git clone https://github.com/edvinolo/b-spline-two-e.git
cd b-spline-two-e

# Build
mkdir build
cd build
cmake ../
make -j$N_jobs

# Create suggested result directories
cd ..
mkdir basis_output quasi_output

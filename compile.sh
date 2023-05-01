#!/bin/bash
root_dir=`pwd`
binary_location=${root_dir}/executables
compiler='gfortran'

cd ${root_dir}/swap
   make clean
   make all
  
if [[ -d ${binary_location} ]]; then
    mv swap.x ${binary_location}
else
    mkdir ${binary_location}
    mv swap.x ${binary_location}
fi

make clean

cd $root_dir/avg
$compiler -g -O3 avg.f -o avg.x
mv avg.x ${binary_location}
cd $root_dir/

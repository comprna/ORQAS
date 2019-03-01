#!/bin/bash
cur_dir=`dirname $0`
pkg_dir=${cur_dir}/../pkg/
bin_dir=${cur_dir}/../bin/
lib_dir=${cur_dir}/../lib/
mkdir -p ${pkg_dir}
mkdir -p ${bin_dir}
mkdir -p ${lib_dir}
os=$1
if [ "$os" = linux ]; then
    salmon_url=https://github.com/kingsfordgroup/sailfish/releases/download/v0.2.7/Salmon-v0.2.7_Ubuntu-12.04.tar.gz
    salmon_dir=${pkg_dir}SalmonBeta-latest_ubuntu-12.04
elif [ "$os" = osx ]; then
    salmon_url=https://github.com/kingsfordgroup/sailfish/releases/download/v0.2.7/Salmon-v0.2.7_MacOSX-10.10.1.tar.gz
    salmon_dir=${pkg_dir}${salmon_url##*/}
    salmon_dir=${salmon_dir%.tar.gz}
else
    echo "Usage: ./include_prerequisites.sh os_type"
    echo "os_type = [ linux | osx ]"
    exit
fi
echo "downloading Saifish..."
salmon_tarball=${pkg_dir}${salmon_url##*/}
wget -P ${pkg_dir} --no-check-certificate -N ${salmon_url}
tar -zxvf ${salmon_tarball} -C ${pkg_dir}
echo "copying Sailfish executables and libraries to bin/ lib/"
cp ${salmon_dir}/bin/* ${bin_dir}
cp ${salmon_dir}/lib/* ${lib_dir}
echo "downloading Star..."
star_url=https://github.com/alexdobin/STAR/archive/STAR_2.4.0j.tar.gz
star_tarball=${pkg_dir}${star_url##*/}
wget -P ${pkg_dir} --no-check-certificate -N ${star_url}
tar -zxvf ${star_tarball} -C ${pkg_dir}
star_dir=${pkg_dir}STAR-${star_url##*/}
star_dir=${star_dir%.tar.gz}
if [ "$os" = linux ]; then
    echo "copying star executable to bin/"
    cp ${star_dir}/bin/Linux_x86_64/* ${bin_dir}
    echo "export PATH..."
    export LD_LIBRARY_PATH=${lib_dir}:$LD_LIBRARY_PATH
elif [ "$os" = osx ]; then
    echo "copying star executable to bin/"
    cp ${star_dir}/bin/MacOSX_x86_64/* ${bin_dir}
    echo "export PATH..."
    export DYLD_FALLBACK_LIBRARY_PATH=${lib_dir}:$DYLD_FALLBACK_LIBRARY_PATH
else
    echo "os_type '$os' not recoginzed!"
    exit
fi
export PATH=${bin_dir}:$PATH
echo "done"


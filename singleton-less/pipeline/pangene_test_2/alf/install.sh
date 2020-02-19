#!/bin/bash

machine=$(uname -m)
os=$(uname -s)

DARWIN_BINARY="alfdarwin.linux32"

install_prefix="/usr/local"
if [ $# -eq 1 ]
then
    if [ $1 = "--help" -o $1 = "-h" -o $1 = "-help" ]
    then
        echo "usage: ./install.sh [install_prefix]"
        exit 0
    fi
    install_prefix=$(dirname $1)/$(basename $1)
fi
mkdir -p $install_prefix/bin

if [ $os = "Linux" ]
then
    if [ $machine = "x86_64" ]
    then
        DARWIN_BINARY="alfdarwin.linux64"
    fi
elif [ $os = "Darwin" ]
then
    if [ $machine = "x86_64" ]
    then
        DARWIN_BINARY="alfdarwin.mac64"
    else
        DARWIN_BINARY="alfdarwin.mac32"
    fi
else
    echo "Operating system not supported!"
fi

echo "installing darwin binary..."
cp bin/alfdarwin bin/$DARWIN_BINARY $install_prefix/bin/
echo "installing ALF..."
cp bin/alfsim $install_prefix/bin/
echo "instaling library..."
mkdir -p $install_prefix/share/alfdarwin/
cp -rf lib $install_prefix/share/alfdarwin/lib
echo "installation complete."
echo "Make sure $install_prefix/bin is in your PATH."

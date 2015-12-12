#!/bin/bash

MODE=$1
if [ -z "$MODE" ]; then
    MODE="quick"
fi

case $MODE in
    "quick")
        rm -rf `find . -name CMakeFiles`
        rm -rf `find . -name CMakeCache.txt`
        rm -rf `find . -name cmake_install.cmake`
        rm -rf `find . -name Makefile`
        rm -rf `find . -name '*~'`
        rm -rf `find . -name '*.o'`
    ;;
    "deep")
        rm -rf `find . -name CMakeFiles`
        rm -rf `find . -name CMakeCache.txt`
        rm -rf `find . -name cmake_install.cmake`
        rm -rf `find . -name Makefile`
        rm -rf `find . -name '*~'`
        rm -rf `find . -name '*.o'`
        rm -rf `find . -name '*.a'`
        rm -rf `find . -name '*.so'`
        rm -rf `find . -name '*.so.*'`
        rm -rf `find . -name '*.mod'`
        rm -rf bin/*
    ;;
esac


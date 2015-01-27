#!/bin/bash

if [ -f $1.c ]
then
    a=$1.c
else
    if [ -f $1.cpp ]
    then
	a=$1.cpp
    else
	echo "Error, neither file $1.c and $1.cpp present!"
	exit
    fi
fi

g++ -o $1 $a -I /home/francesco/Prace/sunpho/Analysis_programs/src -Wall `~/bin/rootlib` -I $(dirname $1) -L/usr/lib/x86_64-linux-gnu/root5.34 -O0 -ggdb3 -llapack -lblas -lfftw3 -std=c++1y $2

#g++ -o $1 $a -I /Users/francesco/Prace/sunpho/Analysis_programs/src -Wall `rootlib` -I $(dirname $1) -O0 -ggdb3

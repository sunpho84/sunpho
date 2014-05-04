#!/bin/bash

if [ "$1" == "" ]
then
    echo Use $0 file
    exit
fi

HPPFILE=$(echo $1|sed 's|cpp|hpp|')
MACRO=_$(echo $(basename $HPPFILE)|tr [:lower:] [:upper:]|sed 's|\.|_|')

echo "#ifndef $MACRO"
echo "#define $MACRO"
if [ -f "$HPPFILE" ]
then
    grep include $HPPFILE
fi

egrep '[[:alnum:]_]+ [[:alnum:]_]+\([[:alnum:]_, &*]+\)' $1|grep -v "//"|grep -v "^{"|sort|uniq|awk -F { '{print $1";"}' 
echo "#endif"

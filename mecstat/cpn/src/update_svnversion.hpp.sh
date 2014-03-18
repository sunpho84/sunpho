#!/bin/bash

if [ "$1" == "" ]
then
    echo "Error, use $0 dir"
    exit 1
fi

#take old version
old_svnversion=$(if [ -f src/svnversion.hpp ];then sed 's|\"||g' src/svnversion.hpp|awk '{print $NF}';else echo 0;fi)

#take new version
cd $1
new_svnversion=$(svnversion 2>&1)
cd $OLDPWD

#if svn unavailable or not svn downloaded reset to old
if  [ "$new_svnversion" == "" ] || \
    [ "$new_svnversion" == exported ] || \
    [ "$new_svnversion" == "Unversioned directory" ]
then
    new_svnversion="$old_svnversion"
fi

#compare and decide if to update
echo \#define SVN_VERSION \""$new_svnversion"\" > test
if [ ! -f src/svnversion.hpp ]
then
    mv test src/svnversion.hpp
else
    diff test src/svnversion.hpp > /dev/null
    if [ "$?" != "0" ]
    then
	mv test src/svnversion.hpp
    else
	rm test
    fi
fi

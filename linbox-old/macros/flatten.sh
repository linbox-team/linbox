#!/bin/sh
## In current directory, create subdir flat/ whose contents are
## hard links to all code files in the tree at the current directory.
## -bds 7/01

LINK="ln -s"

if [ ! -d flat ] ; then mkdir flat ; fi

tmp=`find . -name flat -prune -o -name ARCHIVES -prune -o \( -name "*.h" -a -type f -print \)` 
if [ "x$tmp" != "x" ] ; then 
 	for file in $tmp ; do	
		$LINK $PWD/$file flat/ 
	done
fi

tmp=`find . -name flat -prune -o -name ARCHIVES -prune -o \( -name "*.inl" -a -type f -print \)` 
if [ "x$tmp" != "x" ] ; then 
 	for file in $tmp ; do	
		$LINK $PWD/$file flat/ 
	done
fi

tmp=`find . -name flat -prune -o -name ARCHIVES -prune -o \( -name "*.c" -a -type f -print \)` 
if [ "x$tmp" != "x" ] ; then 
 	for file in $tmp ; do	
		$LINK $PWD/$file flat/ 
	done
fi

tmp=`find . -name flat -prune -o -name ARCHIVES -prune -o \( -name "*.C" -a -type f -print \)` 
if [ "x$tmp" != "x" ] ; then 
 	for file in $tmp ; do	
		$LINK $PWD/$file flat/ 
	done
fi

tmp=`find . -name flat -prune -o -name ARCHIVES -prune -o \( -name "*.cpp" -a -type f -print \)` 
if [ "x$tmp" != "x" ] ; then 
 	for file in $tmp ; do	
		$LINK $PWD/$file flat/ 
	done
fi

tmp=`find . -name flat -prune -o -name ARCHIVES -prune -o \( -name "*.cc" -a -type f -print \)` 
if [ "x$tmp" != "x" ] ; then 
 	for file in $tmp ; do	
		$LINK $PWD/$file flat/ 
	done
fi

# Check for Maple, for maple examples
# Rich Seagraves, 2002-10-31
# a gutted and modified version of ntl-check.m4 by Bradford Hovinen, 2001-06-13
# which was in turn inspired by gnome-bonobo-check.m4 by Miguel de Icaza, 99-04-12
# Stolen from Chris Lahey
# stolen from Manish Singh again
# stolen back from Frank Belew
# stolen from Manish Singh
# Shamelessly stolen from Owen Taylor
# . . . and David begat Solomon who begat . . . (sorry, couldn't help it)

dnl LB_CHECK_MAPLE()

AC_DEFUN([LB_CHECK_MAPLE],
[

AC_ARG_WITH(maple-prefix,[ --with-maple-prefix=PFX],
[maple_prefix="$withval"],[maple_prefix=""])

AC_MSG_CHECKING(for maple support)

if test "x${maple_prefix}" != x; then
	MAPLE_BINPATHIS=${maple_prefix}/`ls $maple_prefix | grep "bin." `
	if test ${MAPLE_BINPATHIS} = ${maple_prefix}; then
		AC_MSG_RESULT(not found)
		cp -f interfaces/maple/Makefile.in.1 interfaces/maple/Makefile.in	
	else
		AC_MSG_RESULT(found)	
		LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${MAPLE_BIN_PATHIS}"
		LD_RUN_PATH="${LD_RUN_PATH}:${MAPLE_BIN_PATHIS}"
		export LD_LIBRARY_PATH
		export LD_RUN_PATH
		MAPLE_LIBS="-L${MAPLE_BINPATHIS} -lmaplec"
		MAPLE_INCLUDES="-I${maple_prefix}/extern/include"
		MAPLE_BUILD_LOC="`pwd`/interfaces/maple"
		cp -f interfaces/maple/lbmaple.mpl.head interfaces/maple/lbmaple.mpl
		echo "`pwd`/intefaces/maple/" >> interfaces/maple/lbmaple.mpl
		cat interfaces/maple/lbmaple.mpl.tail >> interfaces/maple/lbmaple.mpl
       		AC_SUBST(MAPLE_LIBS)
		AC_SUBST(MAPLE_INCLUDES)
		AC_SUBST(MAPLE_BUILD_LOC)
		cp -f interfaces/maple/Makefile.in.2 interfaces/maple/Makefile.in
	fi
else
	AC_MSG_RESULT(not found)
	cp -f interfaces/maple/Makefile.in.1 interfaces/maple/Makefile.in	
fi
])

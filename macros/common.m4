dnl Copyright(c)'2019 LinBox
dnl
dnl ========LICENCE========
dnl This file is part of the library LinBox.
dnl
dnl LinBox is free software: you can redistribute it and/or modify
dnl it under the terms of the  GNU Lesser General Public
dnl License as published by the Free Software Foundation; either
dnl version 2.1 of the License, or (at your option) any later version.
dnl
dnl This library is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
dnl Lesser General Public License for more details.
dnl
dnl You should have received a copy of the GNU Lesser General Public
dnl License along with this library; if not, write to the Free Software
dnl Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
dnl ========LICENCE========
dnl/

dnl Set OPTIM_FLAGS, DEBUG_FLAGS depending on compiler (in CCNAM) and command
dnl line arguments of ./configure (--enable-debug, --enable-warnings and
dnl --enable-profile)
AC_DEFUN([SET_FLAGS],[
    dnl require AC_DEBUG so USE_DEBUG is defined
    AC_REQUIRE([AC_DEBUG])
    dnl require AC_PROFILE so PROF is defined
    AC_REQUIRE([AC_PROFILE])
    dnl require AC_WARNINGS so WARN is defined
    AC_REQUIRE([AC_WARNINGS])
    dnl require AC_COMPILER_NAME so CCNAM is defined
    AC_REQUIRE([AC_COMPILER_NAME])
    
    # --enable-debug ?
    # __LINBOX_DEBUG will be (un)set in linbox/config.h
    #TODO use -fast for icc, -ipa for eko...
    AS_IF([test "x$USE_DEBUG" = "xyes"],
          [OPTIM_FLAGS="-O0"
            DEBUG_FLAGS="-Wall -g -UNDEBUG -DDEBUG"],
          [OPTIM_FLAGS="-O2"
            DEBUG_FLAGS="-Wall -DNDEBUG -UDEBUG"]
          )

    # --enable-profile ?
    AS_IF([test "x$PROF" = "xyes"], [ DEBUG_FLAGS+=" -pg" ])

    # --enable-warnings ?
    AS_IF([test "x$WARN" = "xyes" -o "x$WARN" = "xfull"],
          [AS_CASE([$CCNAM],
                   [eko], [],
                   [gcc*|icc*|clang*], [ DEBUG_FLAGS+=" -Wextra" ],
                   [AS_BOX([Unsupported compiler ($CCNAM). Please file a bug.],[*])]
                   )
          ])

    AS_IF([test "x$WARN" = "xfull"],
          [AS_CASE([$CCNAM],
                   [eko],
                    [],
                   [gcc*|icc*|clang*],
                    [
                      DEBUG_FLAGS+=" -Wuninitialized -Wconversion -Wcast-qual "
                      DEBUG_FLAGS+=" -pedantic -Wshadow -Wpointer-arith "
                      DEBUG_FLAGS+=" -Wwrite-strings -Wno-long-long"
                      AS_CASE([$CCNAM],
                        [icc],
                            [ DEBUG_FLAGS+=" -Wcheck -ansi" ],
                        [gcc*],
                            [
                             DEBUG_FLAGS+=" -Wno-vla"
                             DEBUG_FLAGS+=" -Wcast-align -Wno-variadic-macros"
                            ],
                        [clang*],
                            [
                             DEBUG_FLAGS+=" -Wno-vla-extension -D__STRICT_ANSI__"
                             DEBUG_FLAGS+=" -Wcast-align -Wno-variadic-macros"
                            ])
                    ],
                   [AS_BOX([Unsupported compiler ($CCNAM). Please file a bug.],[*])]
                   )
          ])
    ])

AC_DEFUN([AX_CHECK_COMPILE_FLAG],
[dnl
AS_VAR_PUSHDEF([CACHEVAR],[ax_cv_check_[]_AC_LANG_ABBREV[]flags_$4_$1])dnl
AC_CACHE_CHECK([whether _AC_LANG compiler accepts $1], CACHEVAR, [
  ax_check_save_flags=$[]_AC_LANG_PREFIX[]FLAGS
  _AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS $4 $1"
  AC_COMPILE_IFELSE([AC_LANG_SOURCE(m4_default([$5],[AC_LANG_PROGRAM()]))],
    [AS_VAR_SET(CACHEVAR,[yes])],
    [AS_VAR_SET(CACHEVAR,[no])])
  _AC_LANG_PREFIX[]FLAGS=$ax_check_save_flags])
AS_VAR_IF(CACHEVAR,yes,
  [m4_default([$2], :)],
  [m4_default([$3], :)])
AS_VAR_POPDEF([CACHEVAR])dnl
])dnl AX_CHECK_COMPILE_FLAGS


dnl Append -march=native or -mcpu=native (if recognized by the compiler) to
dnl OPTIM_FLAGS if not present in CXXFLAGS and not cross-compiling and
dnl --without-archnative is not set
AC_DEFUN([ARCH_FLAGS],[
    AC_ARG_WITH(archnative, [AC_HELP_STRING([--without-archnative],
        [do not use -march=native or -mcpu=native (default is to use it if not already present in CXXFLAGS)])])

    AX_CHECK_COMPILE_FLAG([-march=native], [
        AS_CASE([$CXXFLAGS],
                [*-march=*], [], # do nothing if already set in CXXFLAGS
                [AS_IF([test "x${with_archnative}" == "xno"],
                    [], # do nothing if option is set to no
                    [AS_IF([test "${host}" != "${build}" -o "${host}" != "${target}"],
                        [AC_MSG_NOTICE("For efficiency you may want to add a '-march=...' flag in CXXFLAGS")],
                        [AC_MSG_NOTICE("Adding '-march=native' to OPTIM_FLAGS")
                        OPTIM_FLAGS+=" -march=native"])])])
    ], [
        AX_CHECK_COMPILE_FLAG([-mcpu=native], [
            AS_CASE([$CXXFLAGS],
                    [*-cpu=*], [], # do nothing if already set in CXXFLAGS
                    [AS_IF([test "x${with_archnative}" == "xno"],
                        [], # do nothing if option is set to no
                        [AS_IF([test "${host}" != "${build}" -o "${host}" != "${target}"],
                            [AC_MSG_NOTICE("For efficiency you may want to add a '-cpu=...' flag in CXXFLAGS")],
                            [AC_MSG_NOTICE("Adding '-mcpu=native' to OPTIM_FLAGS")
                            OPTIM_FLAGS+=" -mcpu=native"])])])
        ], [], [], [])
    ], [], [])
])

dnl Append -mfpmath=sse to OPTIM_FLAGS on i386 and i686 architecture with SSE
AC_DEFUN([FPMATH_FLAGS],[
    AC_REQUIRE([ARCH_FLAGS])

    BACKUP_CXXFLAGS="${CXXFLAGS}"
    CXXFLAGS="${OPTIM_FLAGS} ${CXXFLAGS}"
    AS_CASE([$target],
            [*i386*|*i686*],
                [AC_RUN_IFELSE([AC_LANG_PROGRAM([[]], [[#ifdef __SSE__
                                                          return 0;
                                                        #else
                                                          return 1;
                                                        #endif
                                                        ]])],
                    [AC_MSG_NOTICE("Adding '-mfpmath=sse' to OPTIM_FLAGS")
                     OPTIM_FLAGS+=" -mfpmath=sse"],
                    [], # either the flag is not recognized by the compiler or
                        # SSE is not avail => do nothing
                    [AC_MSG_NOTICE("If available you may want to add
                     '-mfpmath=sse' to flags")])] # cross-compilation case
                []) # not on i386 nor i686 => do nothing
    CXXFLAGS="${BACKUP_CXXFLAGS}"
    ])

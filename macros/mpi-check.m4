dnl Check for MPI support
dnl Copyright (c) the LinBox group
dnl This file is part of LinBox

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
 dnl

dnl LB_CHECK_MPI()
dnl
dnl Check for MPI support

AC_DEFUN([LB_CHECK_MPI], [
    AC_MSG_CHECKING(for MPI)

    AC_ARG_WITH(mpi,
        [
            AC_HELP_STRING([--with-mpi=yes|no], [
                Use MPI for specific tests or benchmarks.
            ])
        ],
        [],
        []
    )

    AS_IF([ test "x$with_mpi" = "xyes" ],
        [
            BACKUP_CXX=${CXX}
            CXX="mpicxx"

            AC_TRY_RUN(
                [
                    #include <mpi.h>
                    int main(void) {
                        MPI_Init(NULL, NULL);
                        MPI_Finalize();
                        return 0;
                    }
                ],
                [ mpi_found="yes" ],
                [ mpi_found="no" ],
                [ mpi_found="no" ] # Cross compiling
            )

            AS_IF([ test "x$mpi_found" = "xyes" ],
                [
                    AC_DEFINE(HAVE_MPI, 1, [Define if MPI is available])
                    AC_MSG_RESULT(yes)
                    HAVE_MPI=yes
                ],
                [
                    AC_MSG_RESULT(no)
                    CXX=${BACKUP_CXX}
                ]
            )
        ],
        [ AC_MSG_RESULT(not found) ]
    )

    AM_CONDITIONAL(LINBOX_HAVE_MPI, test "x$HAVE_MPI" = "xyes")
])

AC_INIT(mapgd, 0.4.6, matthew.s.ackerman@gmail.com, mapgd.tgz, https://github.com/LynchLab/MAPGD/)
AC_PROG_CXX
AC_LANG(C++)

AC_SUBST(HAVE_GSL, [no])
AC_SUBST(HAVE_OMP, [no])

AC_CHECK_LIB(sqlite3,sqlite3_open,AC_SUBST(HAVE_SQL, [yes]), AC_SUBST(HAVE_SQL, [no]) )
AC_CHECK_PROG(HAVE_MAN, [man], [yes], [no])

AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

AC_OPENMP

if test "$ac_cv_lib_gsl_gsl_blas_dgemm" == yes
then
	AC_SUBST(HAVE_GSL, yes)
else
	AC_MSG_WARN([GSL (gnu scientific libraries) not found. Unable to use relatedness command.])
fi

if test "$ac_cv_prog_cxx_openmp" == "unsupported"
then
	AC_MSG_WARN([No openmp support detected, will not be multi-threaded.])
else
	AC_SUBST(HAVE_OMP, yes)
fi

AC_OUTPUT(Makefile)
#cp confdefs.h my_config.h

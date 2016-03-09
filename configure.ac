AC_INIT(mapgd, 0.4.6, matthew.s.ackerman@gmail.com, mapgd.tgz, https://github.com/LynchLab/MAPGD/)
AC_PROG_CXX
AC_LANG(C++)

#AX_LIB_SQLITE3(v3)
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

#AX_CXX_COMPILE_STDCXX_11()
AC_OPENMP

AC_SUBST(HAVE_SQL, [ ])
AC_SUBST(HAVE_GSL, [ ])
AC_SUBST(HAVE_OMP, [ ])

if test "$ac_cv_lib_gsl_gsl_blas_dgemm" == yes
then
	AC_SUBST(HAVE_GSL, [HAVE_GSL=true])
fi

if test "x$ac_cv_prog_cxx_openmp" == "x"
then
	AC_MSG_WARN([No openmp support detected])
else
	AC_SUBST(HAVE_OMP, [HAVE_OMP=true])
fi

AC_OUTPUT(Makefile)
cp confdefs.h my_config.h

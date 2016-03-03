AC_INIT(myconfig, version-0.1)
AC_PROG_CXX
AC_LANG(C++)
AX_CXX_COMPILE_STDCXX_11()
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])
#AC_OPENMP()
cp confdefs.h my_config.h

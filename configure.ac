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

AM_GNU_GETTEXT

AC_CHECK_LIB([lzo],[lzo1x_1_15_compress],[LZO_LIBS="-llzo"],[HAVE_LZO="no"])

AC_CONFIG_MACRO_DIR([m4])

AC_OPENMP
AX_CXX_COMPILE_STDCXX([11], [noext], [mandatory])
#AX_CXX_COMPILE_STDCXX_11
AX_WITH_HTSLIB
AX_MPI

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

if test "$ac_cv_func_MPI_Init" == "yes"
then
	AC_SUBST(HAVE_MPI, yes)
	AC_SUBST(MPICXX, mpicxx)
else
	AC_MSG_WARN([No mpi support detected, cannot run on multiple nodes.])
fi

if test "$ax_cv_htslib" == yes
then
	AC_SUBST(HAVE_HTS, yes)
else
	AC_SUBST(HAVE_HTS, no)
	AC_MSG_WARN([Htslib not detected, will be unable to create vcf files.])
fi

if test "$ax_cv_htslib_which" = source; then
  Hsource=
  Hinstall='#'
else
  Hsource='#'
  Hinstall=
fi

#if test "$ax_cv_htslib" == yes

AC_SUBST([Hsource])
AC_SUBST([Hinstall])
AC_OUTPUT(Makefile)
#cp confdefs.h my_config.h

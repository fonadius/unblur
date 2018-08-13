#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
#

# Check that Autoconf is of required version
AC_PREREQ([2.62])

# Get started
AC_INIT([Unblur], [1.0.2],[grantt@janelia.hhmi.org])

# We use some macros from http://www.gnu.org/software/autoconf-archive/
m4_include([m4/ax_f90_header.m4])
m4_include([m4/ax_f90_internal_headmod.m4])
m4_include([m4/ax_f90_module_extension.m4])
m4_include([m4/ax_f90_module_flag.m4])
m4_include([m4/ax_f90_module.m4])
m4_include([m4/ax_f90_library.m4])
m4_include([m4/ax_f90_library_setup.m4])

#
AM_INIT_AUTOMAKE([-Wall -Werror] [subdir-objects])

#Set the language to Fortran
AC_LANG([Fortran])

# Determine Fortran compiler to use
AC_PROG_F77()
AC_PROG_FC()

# The macros AC_PROG_FC and AC_PROG_F77 by default set FCFLAGS=-g, which we don't want
FCFLAGS=""
F77FLAGS=""
FFLAGS=""

# Determine C compiler to use
AC_PROG_CC
AC_PROG_CPP
# Get the Fortran compiler to accept the .f90 extension. If it is necessary to use specific flags for this purpouse, do so!
AC_FC_SRCEXT(f90)
AC_FC_SRCEXT(f)

# Make sure free format source code is allowed by compiler, even if .f extension is found
AC_FC_FREEFORM()

# Enable the building of C files
AM_PROG_CC_C_O()

# Check if user specified location for fftw library
AC_ARG_WITH(fftw-libdir,
[[  --with-fftw-libdir=DIR  declare the directory where the fftw3 library is
                          to be found, if its current location is not
                          system accessible ]],
[
    if test "$withval" != "yes" -a "$withval" != ""; then
        LDFLAGS="$LDFLAGS -L$withval"
    fi
])

# Check if user specified location for gsl/cblas library
AC_ARG_WITH(gsl-libdir,
[[  --with-gsl-libdir=DIR  declare the directory where the GSL & GSL CBLAS
                          libraries are to be found, if their current location
                          are not system accessible ]],
[
    if test "$withval" != "yes" -a "$withval" != ""; then
        LDFLAGS="$LDFLAGS -L$withval"
    fi
])



# Check for C libraries we need
AC_LANG_PUSH([C])
AC_SEARCH_LIBS([cblas_sswap],[gslcblas],[],[AC_MSG_ERROR(Could not find your installation of the GSL CBLAS library)])
AC_SEARCH_LIBS([gsl_log1p],[gsl],[],[AC_MSG_ERROR(Could not find your installation of the GSL library)],[-lgslcblas])
AC_SEARCH_LIBS([deflate],[z],[],[AC_MSG_ERROR(Could not find your installation of the zlib library)],[-lz])
AC_SEARCH_LIBS([jbg_enc_init],[jbig],[],[AC_MSG_ERROR(Could not find your installation of the jbig library)],[-ljbig])
AC_SEARCH_LIBS([jpeg_read_scanlines],[jpeg],[],[AC_MSG_ERROR(Could not find your installation of the JPEG library)],[-ljpeg])
AC_SEARCH_LIBS([TIFFOpen],[tiff],[],[AC_MSG_ERROR(Could not find your installation of the TIFF library)],[-ltiff])
AC_LANG_POP([C])

# Check for FGSL. This involves checking for the fgsl module
AX_F90_MODULE_FLAG()
AX_F90_MODULE_EXTENSION()
if test x$av_cv_f90_modext = xunknown; then
 AC_MSG_ERROR([unable to find f90 modules extension])
fi

# Check whether user specified static linking option
static_default="no"
AC_ARG_ENABLE(static, [  --enable-static=[no/yes] turn on static linking
                       [default=$static_default]],, enable_static=$static_default)
if test "x$enable_static" = "xyes"; then
	static_link=true
	AC_MSG_NOTICE(Static linking enabled)
else
	static_link=false
fi
AM_CONDITIONAL([STATIC_LINKING_AM],[test x$static_link = xtrue])

# Check whether user specified the openmp option // add -tcheck to use with intel's thread checker
# preprocessor symbol FMC_EVAL_THREAD means that the scoring function will be threaded
openmp_default="no"
AC_ARG_ENABLE(openmp, [  --enable-openmp=[no/yes] turn on openmp
                       [default=$openmp_default]],, enable_openmp=$openmp_default)
AM_CONDITIONAL([OPENMP_AM],[test "x$enable_openmp" = "xyes"])
if test "x$enable_openmp" = "xyes"; then
	AC_MSG_NOTICE(OpenMP enabled)
else
	AC_MSG_NOTICE(OpenMP disabled)
fi

# Check whether Intel's MKL is available. If not, link against FFTW.
AC_LANG_PUSH([C])
AC_CHECK_LIB([mkl_rt],[fftwf_malloc],[HAVE_MKL="yes"],[HAVE_MKL="no"])
AC_LANG_POP([C])
if test "x$HAVE_MKL" = "xyes"; then
if test "x$FC" = "xifort"; then
	FCFLAGS="$FCFLAGS -mkl"
fi
if test "x$enable_static" = "xyes"; then
	if test "x$FC" = "xnagfor"; then
		AC_MSG_NOTICE(Static linking against Intel's MKL is not supported for the NAG compiler)
		HAVE_MKL="no"
	else
		if test "x$enable_openmp" = "xyes"; then
			LIBS="$LIBS  -Wl,--start-group \$(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \$(MKLROOT)/lib/intel64/libmkl_intel_thread.a \$(MKLROOT)/lib/intel64/libmkl_core.a \$(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm"
		else
			LIBS="$LIBS  -Wl,--start-group \$(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \$(MKLROOT)/lib/intel64/libmkl_core.a \$(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm"
		fi
	fi
else
if test "x$FC" = "xifort"; then
	LIBS="$LIBS -mkl"
else
	LIBS="$LIBS -lmkl_rt"
	LDFLAGS="$LDFLAGS -L\$(MKLROOT)/lib/intel64"
fi
fi
fi

if test "x$HAVE_MKL" = "xno"; then
	AC_LANG_PUSH([C])
	AC_SEARCH_LIBS([fftwf_malloc],[fftw3f],[],[AC_MSG_ERROR(could not find your installation of the FFTW library (single-precision). FFTW >=3.3 needed.)])
	AC_LANG_POP([C])
	AC_MSG_NOTICE(Using FFTW rather than MKL for fast Fourier transforms)
else
	AC_MSG_NOTICE(Will use Intel MKL rather than FFTW for fast Fourier transforms)
fi


# Set a conditional that automake can use to add compiler-specific flags
AM_CONDITIONAL([FC_IS_GFORTRAN],[test $FC = gfortran])
AM_CONDITIONAL([FC_IS_INTEL],[test $FC = ifort])
AM_CONDITIONAL([FC_IS_PORTLAND],[test $FC = pgfortran])
AM_CONDITIONAL([FC_IS_NAG],[test $FC = nagfor])

# Set a conditional so that automake can know about debugging
AM_CONDITIONAL([DEBUG],[test "x$enable_debug" = "xyes"])

# This needed with the library manip in automake
AC_PROG_RANLIB

# Make sure OpenMP flags are set. This does not work, and I'm not sure why.
# AC_OPENMP()


# Check whether user specified debugging option
debug_default="no"
AC_ARG_ENABLE(debug, [  --enable-debug=[no/yes] turn on debugging
                       [default=$debug_default]],, enable_debug=$debug_default)
if test "x$enable_debug" = "xyes"; then
	AC_MSG_NOTICE(Debugging enabled)
else
	AC_MSG_NOTICE(Debugging disabled)
fi

# Check whether user specified optimisation option
optimisations_default="yes"
AC_ARG_ENABLE(optimisations, [  --enable-optimisations=[no/yes] turn on optimisations
                       [default=$optimisations_default]],, enable_optimisations=$optimisations_default)
if test "x$enable_optimisations" = "xyes"; then
if test "x$FC" = "xnagfor"; then
	FCFLAGS="$FCFLAGS -O4" 
	FFLAGS="$FFFLAGS -O4"
elif test "x$FC" = "xpgfortran"; then
	FCFLAGS="$FCFLAGS -O4"
	FFLAGS="$FFLAGS -O4"
elif test "x$FC" = "xifort"; then
	FCFLAGS="$FCFLAGS -O3 -no-prec-div"
	FFLAGS="$FFLAGS -O3 -no-prec-div"
else
	FCFLAGS="$FCFLAGS -O3" 
	FFLAGS="$FFFLAGS -O3" 
fi
	AC_MSG_NOTICE(Optimisations enabled)
else
	FCFLAGS="$FCFLAGS -O0"
    FFLAGS=" -O0"
	AC_MSG_NOTICE(Optimisations disabled)
fi

# Check whether user specified instruction set option
latest_instructions_default="yes"
AC_ARG_ENABLE(latest_instructions, [  --enable-latest-instruction-set=[no/yes] use the latest available CPU instruction sets on the current processor
                       [default=$latest_instructionss_default]],, enable_latest_instruction_set=$latest_instructrions_default)
if test "x$enable_latest_instruction_set" = "xyes"; then
if test "x$FC" = "xnagfor"; then
	FCFLAGS="$FCFLAGS" 
	FFLAGS="$FFLAGS"
elif test "x$FC" = "xpgfortran"; then
	FCFLAGS="$FCFLAGS"
	FFLAGS="$FFLAGS"
elif test "x$FC" = "xifort"; then
	FCFLAGS="$FCFLAGS -xHost"
	FFLAGS="$FFLAGS -xHost"
else
	FCFLAGS="$FCFLAGS -march=native" 
	FFLAGS="$FFLAGS -march=native" 
fi
	AC_MSG_NOTICE(Latest instruction sets enabled)
else
	FCFLAGS="$FCFLAGS"
    FFLAGS="$FFLAGS"
	AC_MSG_NOTICE(Latest instruction sets disabled)
fi



# Check whether user specified profiling option
profiling_default="no"
AC_ARG_ENABLE(profiling, [  --enable-profiling=[no/yes] turn on gprof profiling
                       [default=$profiling_default]],, enable_profiling=$profiling_default)
AM_CONDITIONAL([PROFILING_AM],[test "x$enable_profiling" = "xyes"])
if test "x$enable_profiling" = "xyes"; then
	AC_MSG_NOTICE(Profiling enabled)
else
	AC_MSG_NOTICE(Profiling disabled)
fi




# Check whether user specified the openmp profiling option
openmp_profiling_default="no"
AC_ARG_ENABLE(openmp-profiling, [  --enable-openmp-profiling=[no/yes] turn on openmp profiling
                      	 [default=$openmp_profiling_default]],, enable_openmp_profiling=$openmp_default)
if test "x$enable_openmp_profiling" = "xyes"; then
	FC="kinst-ompp "$FC
	F77="kinst-ompp "$F77
	F90="kinst-ompp "$F90
	AC_MSG_NOTICE(OpenMP profiling enabled)
else
	AC_MSG_NOTICE(OpenMP profiling disabled)
fi

# Check whether we want to leave the run-time checks programmed into frealix
AC_MSG_CHECKING(--enable-runtime-checks argument)
AC_ARG_ENABLE(runtime_checks,
    [  --enable-runtime-checks         Enable runtime checks, especially useful during development, but may slow down execution.],
    [enable_runtime_checks=$enableval],
    [enable_runtime_checks="yes"])
AC_MSG_RESULT($enable_runtime_checks)
if test "$enable_runtime_checks" = "no"; then
  FCFLAGS="$FCFLAGS -DSKIP_RUNTIME_CHECKS"
fi

# Check whether user specified the coarray option
if test "x$FC" = "xifort"; then
coarray_default="no"
else
coarray_default="no"
fi
AC_ARG_ENABLE(coarray, [  --enable-coarray=[no/yes] turn on coarray support
                       [default=$coarray_default]],, enable_coarray=$coarray_default)
AM_CONDITIONAL([COARRAY_AM],[test "x$enable_coarray" = "xyes"])
if test "x$enable_coarray" = "xyes"; then
	AC_MSG_NOTICE(Coarrays enabled)
else
	AC_MSG_NOTICE(Coarrays disabled)
fi


# Get SVN information from the shell command
define([svnversion], esyscmd([svnversion -n]))
AC_MSG_NOTICE(svn version = svnversion)
FCFLAGS="$FCFLAGS"


# Debug
AC_MSG_NOTICE(FCFLAGS at the end: $FCFLAGS)

AC_CONFIG_FILES([Makefile])

# Checks for library functions.
AC_OUTPUT

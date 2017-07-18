# Set up conditional compilation so that one Makefile can handle
# multiple platforms/compilers. Branch according to an environmental
# variable F90. I wish someone would show me a better way of doing this.
#

# Shell
SHELL		= /bin/bash

# enumlib library paths
LBD1 = ../../../symlib/src

FOUND = false
ifeq (${F90},gfortran)  # gfortran compiler
  ifeq (${DEBUG},false)
     FFLAGS = -O3 -ffree-line-length-none -fno-underscoring -I${LBD1} 
     FOUND = true
  else
     FFLAGS = -fPIC -g -fbounds-check -Wall -ffree-line-length-none -fno-underscoring -I${LBD1} 
     FOUND = true
  endif
endif

ifeq (${F90},ifc)  # Intel compiler
  FFLAGS = -fPIC -g -error-limit 7 -traceback -check bounds -warn  -I${LBD1}
  FOUND = true
endif

ifeq (${F90},ifort)  # Intel compiler
  ifeq (${DEBUG},false)
     FFLAGS =  -fPIC -O3 -I${LBD1} 
     FOUND = true
  else
#  F90 =  /opt/intel/fc/10.0.016/bin/ifort
     FFLAGS =  -fPIC -g -debug -error-limit 7 -heap-arrays -traceback -check bounds -warn -I${LBD1}
     FOUND = true
#-prof-use -prof-dir .
   endif
endif


ifeq (${F90},xlf90) # IBM compiler
  FFLAGS = -g -C -fPIC -qsuffix=f=f90  -I${LBD1}
  FOUND = true
endif
ifeq (${F90},f95) # Absoft PPC compiler
#  FFLAGS =  -profile -p ${LBDR} 
  FFLAGS = -g -Rb -Rc  -fPIC -nodefaultmod -p ${LBD1} #-ea
# B80  show entry in subprograms ; Rb bounds; Rc array conformance;
# z2 warning level
# -O3 optimization
# -ea stop after one error 
  FOUND = true
endif
ifeq (${F90},)  # If the variable isn't set, make sure compilation bombs
error:
	echo Error environmental variable "F90" not set!
endif
ifneq (${FOUND},true) # If it's set but no matching flags then bomb
error:	
	echo Error: makefile doesn\'t have flags for this compiler
endif

# if type "python" > /dev/null; then
#     python irevision.py derivative_structure_generator.f90 version.enum	
# fi

SRC = sp_hnfs.f90

OBJS = ${SRC:.f90=.o} 
LIBS =  $(LBD1)/libcomparestructs.a $(LBD1)/libutils.a $(LBD1)/libsym.a \
         $(LBD1)/librational.a $(LBD1)/libcombinatorics.a

.SUFFIXES :  
.SUFFIXES :  .f .f90 .f95 .o


libsphnf.a: ${OBJS}
	ar ru $@ $?
	ranlib  $@

all: libsphnf.a

.f95.o : 
	${F90} ${FFLAGS} -c $<
.f90.o : 
	${F90} ${FFLAGS} -c $<
.f.o : 
	${F90} -c $<



CLEAN  = *.o *.mod *.a *.x svn-commit.*
clean :
	rm -f ${CLEAN}
clobber : 
	rm -f  ${CLEAN}  *~ \#*
	make

SHELL = /bin/bash
MAIN = SVMC

FC       = ifort
#INCPATH  = /xnilu_wrk/flex_wrk/bin64/grib_api/include
#INCPATH1 = /usr/include
#INCPATH_NCDF =  /homevip/flexpart/gcc-5.4.0/include
INCPATH1 = ${EBROOTNETCDFMINFORTRAN}/include
INCPATH2 = /cluster/software/ecCodes/2.9.2-intel-2018b/include
#/usr/include

#LIBPATH1 = /xnilu_wrk/flex_wrk/bin64/grib_api/lib/
#LIBPATH1 = /homevip/flexpart/gcc-5.4.0/lib
#/usr/lib/x86_64-linux-gnu/ 
#/homevip/flexpart/gcc-5.4.0/lib
#LIBPATH2 =   /usr/lib/x86_64-linux-gnu/
LIBPATH1 = /cluster/software/netCDF-Fortran/4.4.4-intel-2018b/lib
LIBPATH2 = /cluster/software/ecCodes/2.9.2-intel-2018b/lib
LIBS = -leccodes_f90 -leccodes -lm -ljasper -lnetcdff -DUSE_NCF

FFLAGS   = -I$(INCPATH1) -I$(INCPATH2) -O2 -g -cpp -m64  -mcmodel=medium -lnetcdff -traceback
#FFLAGS   =   -I$(INCPATH1) -I$(INCPATH_NCDF) -g  -fbacktrace -mcmodel=medium -fbounds-check -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan
LDFLAGS  = $(FFLAGS)  -L$(LIBPATH1) -L$(LIBPATH2) $(LIBS)  #-Wl,-rpath,$(LIBPATH1) $(LIBPATH2) $(LIBS) 
#LDFLAGS  = $(FFLAGS)  -L$(LIBPATH1) -lgrib_api_f90 -lgrib_api -lm -ljasper -lnetcdff
#-lnetcdff
#-L$(LIBPATH2)


MODOBJS = \
readpara_mod.o \
readclim_mod.o \
phydro_mod.o \
alloc_mod.o

OBJECTS = \
svmc.o \
yassofortran.o


$(MAIN): $(MODOBJS) $(OBJECTS)
        $(FC) *.o -o $(MAIN) $(LDFLAGS)

$(OBJECTS): $(MODOBJS)

%.o: %.f90
        $(FC) -c $(FFLAGS) $<

clean:
        rm *.o *.mod

cleanall:
        rm *.o *.mod $(MAIN)
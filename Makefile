# Makefile
# Bernd Bruegmann, 12/99, 5/02,  Wolfgang Tichy 4/2005
# Builds the sgrid executable based on the file MyConfig
# See http://www.gnu.org/software/make/manual for the manual of GNU make

# where am I
UNAME := $(shell uname)
TOP   := $(shell pwd)
# to shorten the output: about to be removed for non-standard directories
TOP   := ../../..

# name the fruit of our labor
EXEC = sgrid
EXECDIR = $(TOP)/exe


# variables common to all setups
CC = gcc	# gcc or icc

INCS = -I$(TOP)/src/main/main
LIBS = -L$(TOP)/lib
SPECIALLIBS =
libsys = -lm

DFLAGS =
OFLAGS =
WARN = # -Wall


# --------------------------------------------------------------------------
# different machines and environments

# Linux
OFLAGS = -O2
#CC = icc:  quite different compared to Cactus: icc gives slower results!    


# --------------------------------------------------------------------------
# some libraries are currently required
libpaths = src/main/MemoryMan 
libpaths += src/utility/output src/utility/evolve
libpaths += src/utility/Spectral
libpaths += src/utility/Coordinates
#libpaths += src/utility/NumericUtils

# --------------------------------------------------------------------------
# the user choses the libraries and some options in the file MyConfig
include MyConfig

# --------------------------------------------------------------------------
# manage how the sgrid sources are compiled

# note that the order matters, e.g.
# main has to go last since it has to be compiled last
libpaths += src/main/main

# extract the list of directory names
libdirs = $(dir $(libpaths))

# extract list of names
libnames = $(notdir $(libpaths))

# make the list of libraries
liblist := $(foreach libname,$(libnames),-l$(libname))

# remove -lmain from that list
liblist := $(subst -lmain,,$(liblist))

# make final list of libraries that is passed to the linker
# uses standard hack to resolve interdependencies by repeating the libraries
# system libraries go in the end
LIBS += $(MPIDIRL) $(liblist) $(liblist)
LIBS += $(SPECIALLIBS) $(libsys) $(MPILIBS)


# make the list of include files that will be automatically included for each
# module
libincludes := $(foreach libpath,$(libpaths),\
	$(libpath)/sgrid_$(notdir $(libpath)).h)

# define the automatic configuration files
autoinclude = src/main/main/sgrid_automatic_include.h
autoinitial = src/main/main/sgrid_automatic_initialize.c
autotext    = \/\* automatically generated from MyConfig \*\/



# --------------------------------------------------------------------------
# some of the above variables are meant to be global, so pass them on
# to the shell 
CFLAGS = $(DFLAGS) $(OFLAGS) $(INCS) $(MPIDIRI) $(HDF5DIRI) $(WARN)
export


# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# default target
sgrid: $(autoinclude) $(autoinitial)
	for X in $(libnames); do mkdir -p lib/obj/$$X; done
	for X in $(libpaths); do $(MAKE) -C $$X; done


# --------------------------------------------------------------------------
# other targets 

# if there is no MyConfig file, use the example provided in doc
MyConfig:
	-if test ! -f MyConfig; then cp doc/MyConfig.example MyConfig; fi 

# automatic configuration files
$(autoinclude): MyConfig
	echo $(autotext) > $(autoinclude) 
	for X in $(libincludes); do \
	  echo \#include \"$(TOP)/$$X\" >> $(autoinclude); \
	done

$(autoinitial): MyConfig
	echo $(autotext) > $(autoinitial) 
	for X in $(libnames); do \
	  echo sgrid\_$$X\(\)\; >> $(autoinitial); \
	done


# create tar file
tar:
	cd ..; tar czf sgrid.tgz --exclude lib --exclude exe ./sgrid

# take a fresh look at things
clean:
	-rm -r lib 
	-rm $(autoinclude)
	-rm $(autoinitial)

# remove emacs backup files
cleantilde:
	find . -name "*~" -exec rm {} \;


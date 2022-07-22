# Makefile
# Bernd Bruegmann, 12/99, 5/02,  Wolfgang Tichy 4/2005
# Builds the sgrid executable based on the file MyConfig
# See http://www.gnu.org/software/make/manual for the manual of GNU make

# where am I
UNAME := $(shell uname)
TOP   := $(shell pwd)
## to shorten the output: about to be removed for non-standard directories
#TOP   := ../../..

# name the fruit of our labor
EXEC = sgrid
EXECDIR = $(TOP)/exe
RELAPROJECTDIR = src/Projects
PROJECTDIR = $(TOP)/$(RELAPROJECTDIR)

# variables common to all setups
CC = gcc	# gcc or icc
CXX =		# g++ or icc
CLINKER =	# will be used only in src/main/main/Makefile for linking
AR = ar		# ar command we use to build library from object files

INCS = -I$(TOP)/src/main/main
LIBS = -L$(TOP)/lib
SPECIALINCS =
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
libpaths = src/main/MemoryMan src/utility/ParManipulator
libpaths += src/utility/Spectral
libpaths += src/utility/boundary
libpaths += src/utility/output src/utility/evolve
libpaths += src/utility/NumberChecker
libpaths += src/utility/checkpoint

# --------------------------------------------------------------------------
# the user choses the libraries and some options in the file MyConfig
projects =#
include MyConfig

# --------------------------------------------------------------------------
# set variable projectnames from git targets in projects. NOTE: there must
# be a / just before the actual projectname (e.g. mars.fau.edu:/DNSdata)
projectbasenames = $(basename $(projects))
projectnames = $(notdir $(projectbasenames))

# set projectpaths and add them to libpaths
projectpaths = $(addprefix $(RELAPROJECTDIR)/,$(projectnames))
libpaths += $(projectpaths)

# --------------------------------------------------------------------------
# some more libraries are currently required, those need to be last
libpaths += src/utility/Coordinates
libpaths += src/utility/NumericUtils
## check if we have old NumericUtils, if not add new numerics
#ifneq ($(findstring NumericUtils,$(libpaths)),NumericUtils)
#libpaths += src/utility/numerics
#endif

# --------------------------------------------------------------------------
# set CXX and CLINKER to CC if they are not set in MyConfig
ifeq ($(CXX),)
CXX = $(CC)
endif

ifeq ($(CLINKER),)
CLINKER = $(CC)
endif

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
CFLAGS = $(DFLAGS) $(OFLAGS) $(INCS) $(MPIDIRI) $(HDF5DIRI) $(SPECIALINCS) $(WARN)
export


# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# default target
sgrid: $(autoinclude) $(autoinitial)
	@echo CC=$(CC)
	@echo CXX=$(CXX)
	@echo CLINKER=$(CLINKER)
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
	  echo int sgrid\_$$X\(\)\; >> $(autoinitial); \
	done
	echo "/* call sgrid initialization functions: */" >> $(autoinitial); \
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

# print some vars
printvars:
	@echo projects=$(projects)
	@echo PROJECTDIR=$(PROJECTDIR)
	@echo projectpaths=$(projectpaths)
	@echo projectnames=$(projectnames)
	@echo libpaths=$(libpaths)


# targets to get git projects
git_clone:
	@echo ==================== Cloning sgrid projects ====================
	-for X in $(projects); do N=$$(basename $$X .git); printf "==== %s ====\n" $$N; git clone $$X $(PROJECTDIR)/$$N; done
#	@$(MAKE) install_git_hooks

#git_pull: install_git_hooks
git_pull:
	@echo ====================== main part of sgrid ======================
	git pull
	@echo ======================== sgrid projects ========================
	@for X in $(projects); do N=$$(basename $$X .git); if [ -d "$(PROJECTDIR)/$$N" ]; then printf "==== %s ====\n" $$N; cd $(PROJECTDIR)/$$N; git pull; fi done

git_status:
	@echo ====================== main part of sgrid ======================
	git status -uno
	@echo ======================== sgrid projects ========================
	@for X in $(projects); do N=$$(basename $$X .git); if [ -d "$(PROJECTDIR)/$$N" ]; then printf "==== %s ====\n" $$N; cd $(PROJECTDIR)/$$N; git status -uno; fi done


# remove code that is not needed once the corresponding libs have been built
rm_MemoryMan_code:
	find src/main/MemoryMan/ -name "*.c*" -print -exec rm -rf '{}' \;
	find src/main/MemoryMan/ -name "*.m*" -print -exec rm -rf '{}' \;
	echo -e "ls:\n\tls *.h" > donothing_Makefile
	find src/main/MemoryMan/ -name "Makefile" -print -exec cp -f donothing_Makefile '{}' \;
	rm -f donothing_Makefile

rm_Math_code:
	rm -rf src/Math

rm_utility_code:
	find src/utility/ -name "*.c*" -print -exec rm -rf '{}' \;
	find src/utility/ -name "*.m*" -print -exec rm -rf '{}' \;
	echo -e "ls:\n\tls *.h" > donothing_Makefile
	find src/utility/ -name "Makefile" -print -exec cp -f donothing_Makefile '{}' \;
	rm -f donothing_Makefile

rm_physics_code:
	find src/physics/ -name "*.c*" -print -exec rm -rf '{}' \;
	find src/physics/ -name "*.m*" -print -exec rm -rf '{}' \;
	echo -e "ls:\n\tls *.h" > donothing_Makefile
	find src/physics/ -name "Makefile" -print -exec cp -f donothing_Makefile '{}' \;
	rm -f donothing_Makefile

rm_code:
	make rm_MemoryMan_code
	make rm_Math_code
	make rm_utility_code
	make rm_physics_code

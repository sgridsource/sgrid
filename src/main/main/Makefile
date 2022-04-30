# main/Makefile

NAME := main

OBJS := sgrid_main.o main.o parameters.o skeleton.o tensors.o utilities.o \
	variables.o endianIO.o sgrid_MPI.o


# main is special, so we do not use Makefile.subdirs here
# By default we build an executable. Only if CFLAGS contains -DLIBSGRID we
# build a library

# directory of the libraries
LIBDIR  = $(TOP)/lib

# directory for the object files of this library
OBJDIR  = $(LIBDIR)/obj/$(NAME)
OBJSWITHPATH = $(addprefix $(OBJDIR)/,$(OBJS))


# default case where we build an executable:
ifneq ($(findstring -DLIBSGRID,$(CFLAGS)),-DLIBSGRID)
# main target is the executable
# note that $@ is the matched target
$(EXECDIR)/$(EXEC): $(OBJSWITHPATH) $(wildcard $(LIBDIR)/*)
	$(CLINKER) $(CFLAGS) -o $@ $(OBJSWITHPATH) $(LIBS)

# special case where we build the library libsgrid:
else
$(LIBDIR)/libsgrid.a: $(OBJSWITHPATH) $(wildcard $(LIBDIR)/*)
	-$(eval ALLOBJS = $(foreach X,$(libnames),$(wildcard $(LIBDIR)/obj/$(X)/*.o)))
	$(AR) rcs $@ $(ALLOBJS)
endif


# secondary target: build object files
# note that $@ is the matched target, 
# $* is the stem of the target defined by %
$(OBJDIR)/%.o: %.c sgrid.h sgrid_$(NAME).h $(NAME).h Makefile
	$(CC) $(CFLAGS) -o $@ -c $*.c


# dependencies for automatic files
sgrid.h: sgrid_automatic_include.h
main.o: sgrid_automatic_initialize.c

# ScalarOnKerr/Makefile

NAME := ScalarOnKerr
OBJS := sgrid_$(NAME).o ScalarOnKerr.o Kerr.o Kerr3d.o

include $(TOP)/Makefile.subdirs

Kerr.c: Kerr.m
	math < Kerr.m > /dev/null

Kerr3d.c: Kerr3d.m
	math < Kerr3d.m > /dev/null
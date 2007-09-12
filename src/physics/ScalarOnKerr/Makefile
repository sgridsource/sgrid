# ScalarOnKerr/Makefile

NAME := ScalarOnKerr
OBJS := sgrid_$(NAME).o ScalarOnKerr.o Kerr.o

include $(TOP)/Makefile.subdirs

Kerr.c: Kerr.m
	math < Kerr.m > /dev/null

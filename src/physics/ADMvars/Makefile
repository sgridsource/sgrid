# ADMvars/Makefile

NAME := ADMvars
OBJS := sgrid_$(NAME).o ADMvars.o ADMconstraints.o utility.o CheckIfFinite.o

include $(TOP)/Makefile.subdirs

ADMconstraints.c: ADMconstraints.m
	math < ADMconstraints.m > /dev/null

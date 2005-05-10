# SingleBH/Makefile

NAME := SingleBH
OBJS := sgrid_$(NAME).o SingleBH.o SingleBHKS.o SingleBHisotropic.o

include $(TOP)/Makefile.subdirs

SingleBHKS.c: SingleBHKS.m
	math < SingleBHKS.m > /dev/null


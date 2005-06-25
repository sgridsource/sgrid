# Z4secondO/Makefile
# Wolfgang Tichy 24.6.2005

NAME := Z4secondO
OBJS := sgrid_$(NAME).o Z4secondO.o Z4secondO_rhs.o

include $(TOP)/Makefile.subdirs

Z4secondO_rhs.c: Z4secondO_rhs.m
	math < Z4secondO_rhs.m > /dev/null

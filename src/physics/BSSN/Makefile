# BSSN/Makefile

NAME := BSSN
OBJS := sgrid_$(NAME).o BSSN.o BSSN_rhs.o BSSN_init.o

include $(TOP)/Makefile.subdirs

BSSN_init.c: BSSN_init.m
	math < BSSN_init.m > /dev/null

BSSN_rhs.c: BSSN_rhs.m
	math < BSSN_rhs.m > /dev/null

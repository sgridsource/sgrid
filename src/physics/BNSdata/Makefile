# BNSdata/Makefile

NAME := BNSdata
OBJS := sgrid_$(NAME).o BNSdata.o BNS_CTS.o setADMvars.o TOVstar.o BNSgrid.o

include $(TOP)/Makefile.subdirs

BNS_CTS.c: BNS_CTS.m
	math < BNS_CTS.m > /dev/null

setADMvars.c: setADMvars.m
	math < setADMvars.m > /dev/null

# BNSdata/Makefile

NAME := BNSdata
OBJS := sgrid_$(NAME).o BNSdata.o BNS_CTS.o setADMvars.o TOVstar.o BNSgrid.o BNS_compute_new_q.o BNS_compute_new_q_atXYZ.o BNS_set_restmassintegrand.o

include $(TOP)/Makefile.subdirs

BNS_CTS.c: BNS_CTS.m
	math < BNS_CTS.m > /dev/null

setADMvars.c: setADMvars.m
	math < setADMvars.m > /dev/null

BNS_compute_new_q.c: BNS_compute_new_q.m
	math < BNS_compute_new_q.m > /dev/null

BNS_compute_new_q_atXYZ.c: BNS_compute_new_q_atXYZ.m
	math < BNS_compute_new_q_atXYZ.m > /dev/null

BNS_set_restmassintegrand.c: BNS_set_restmassintegrand.m
	math< BNS_set_restmassintegrand.m  > /dev/null

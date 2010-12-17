# BNSdata/Makefile

NAME := BNSdata
OBJS := sgrid_$(NAME).o BNSdata.o BNS_CTS.o setADMvars.o TOVstar.o \
	BNSgrid.o BNS_compute_new_q.o BNS_compute_new_q_atXYZ.o \
 	BNS_set_restmassintegrand.o BNS_set_J_ADM_VolInt_integrand.o \
	BNS_set_M_ADM_VolInt_integrand.o BNS_Interpolate_ADMvars.o \
	set_BNSdata_Sigma_BCs.o BNS_BCs.o

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
	math < BNS_set_restmassintegrand.m > /dev/null

BNS_set_J_ADM_VolInt_integrand.c: BNS_set_J_ADM_VolInt_integrand.m
	math < BNS_set_J_ADM_VolInt_integrand.m > /dev/null

BNS_set_M_ADM_VolInt_integrand.c: BNS_set_M_ADM_VolInt_integrand.m
	math < BNS_set_M_ADM_VolInt_integrand.m > /dev/null

set_BNSdata_Sigma_BCs.c: set_BNSdata_Sigma_BCs.m
	math < set_BNSdata_Sigma_BCs.m > /dev/null

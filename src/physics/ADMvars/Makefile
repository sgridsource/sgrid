# ADMvars/Makefile

NAME := ADMvars
OBJS := sgrid_$(NAME).o ADMvars.o ADMconstraints.o utility.o CheckIfFinite.o ADMenergy_spheric_intergrand.o

include $(TOP)/Makefile.subdirs

ADMconstraints.c: | ADMconstraints.m
	math < ADMconstraints.m > /dev/null

ADMenergy_spheric_intergrand.c: | ADMenergy_spheric_intergrand.m
	math < ADMenergy_spheric_intergrand.m > /dev/null

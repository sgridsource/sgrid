# Coordinates/Makefile

NAME := Coordinates

OBJS := sgrid_$(NAME).o $(NAME).o cartesianDerivs.o doubleCovering.o\
	coordFilters.o ComplexFunctions.o transfDerivs.o findXYZ_of_xyz.o\
	coordtrans_AnsorgNS0.o

include $(TOP)/Makefile.subdirs

coordtrans_AnsorgNS0.c: coordtrans_AnsorgNS0.m
	math < coordtrans_AnsorgNS0.m > /dev/null

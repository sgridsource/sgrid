# Coordinates/Makefile

NAME := Coordinates

OBJS := sgrid_$(NAME).o $(NAME).o cartesianDerivs.o doubleCovering.o\
	coordFilters.o ComplexFunctions.o transfDerivs.o findXYZ_of_xyz.o\
	Coordinates_AnsorgNS.o coordtrans_AnsorgXRphi.o\
	coordtrans_AnsorgNS0.o coordtrans_AnsorgNS1.o\
	coordtrans_AnsorgNS2.o coordtrans_AnsorgNS3.o\
	coordtrans_Bfunc.o

include $(TOP)/Makefile.subdirs

coordtrans_AnsorgXRphi.c: coordtrans_AnsorgXRphi.m
	math < coordtrans_AnsorgXRphi.m > /dev/null

coordtrans_AnsorgNS0.c: coordtrans_AnsorgNS0.m
	math < coordtrans_AnsorgNS0.m > /dev/null

coordtrans_AnsorgNS1.c: coordtrans_AnsorgNS1.m
	math < coordtrans_AnsorgNS1.m > /dev/null

coordtrans_AnsorgNS2.c: coordtrans_AnsorgNS2.m
	math < coordtrans_AnsorgNS2.m > /dev/null

coordtrans_AnsorgNS3.c: coordtrans_AnsorgNS3.m
	math < coordtrans_AnsorgNS3.m > /dev/null

coordtrans_Bfunc.c: coordtrans_Bfunc.m
	math < coordtrans_Bfunc.m > /dev/null

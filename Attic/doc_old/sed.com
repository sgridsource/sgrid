Command1:
cp $1 l ; sed s/W.\ T.\*/"Wolfgang Tichy, April 2005  \&  "'&'/ l > $1

forCommand1:
main/
main/MemoryMan/
main/MemoryMan/Makefile
main/MemoryMan/MemoryMan.h
main/MemoryMan/grid.c
main/MemoryMan/print.c
main/MemoryMan/sgrid_MemoryMan.c
main/MemoryMan/sgrid_MemoryMan.h
main/MemoryMan/sgrid_MemoryMan_loops.h
main/MemoryMan/storage.c
main/main/
main/main/Makefile
main/main/main.c
main/main/main.h
main/main/parameters.c
main/main/sgrid.h
main/main/sgrid_automatic_include.h
main/main/sgrid_automatic_initialize.c
main/main/sgrid_main.c
main/main/sgrid_main.h
main/main/skeleton.c
main/main/tensors.c
main/main/utilities.c
main/main/variables.c
physics/
physics/testwave/
physics/testwave/Makefile
physics/testwave/sgrid_testwave.c
physics/testwave/sgrid_testwave.h
physics/testwave/testwave.c
physics/testwave/testwave.h
physics/testwave/testwave.par
physics/testwave/wave.c
utility/
utility/BAMoutput/
utility/BAMoutput/CVS/
utility/BAMoutput/CVS/Entries
utility/BAMoutput/CVS/Repository
utility/BAMoutput/CVS/Root
utility/BAMoutput/HDF5Bugs.txt
utility/BAMoutput/HDF5Readme.txt
utility/BAMoutput/Makefile
utility/BAMoutput/bam_output.c
utility/BAMoutput/bam_output.h
utility/BAMoutput/output.c
utility/BAMoutput/output.h
utility/BAMoutput/output0d.c
utility/BAMoutput/output1d.c
utility/BAMoutput/output2d.c
utility/BAMoutput/output3d.c
utility/BAMoutput/outputFMR.c
utility/BAMoutput/outputHDF5.c
utility/BAMoutput/prgnuplot.c
utility/BAMoutput/prstdio.c
utility/BAMoutput/prxgraph.c
utility/BAMps/
utility/BAMps/CVS/
utility/BAMps/CVS/Entries
utility/BAMps/CVS/Repository
utility/BAMps/CVS/Root
utility/BAMps/Makefile
utility/BAMps/bam_ps.c
utility/BAMps/bam_ps.h
utility/BAMps/diffdirect.c
utility/BAMps/diffmatrix.c
utility/BAMps/ft.c
utility/BAMps/ft.par
utility/BAMps/ps.c
utility/BAMps/ps.h
utility/BAMps/ps_bssn_rhs.m
utility/BAMps/pswave.par
utility/BAMps/pswave1d.par
utility/BAMps/pswavf1d.par
utility/BAMps/pswavr1d.par
utility/BAMps/pswavsph.par
utility/BAMps/wave.c
utility/Coordinates/
utility/Coordinates/.ispell.joe.bak
utility/Coordinates/Coordinates.c
utility/Coordinates/Coordinates.h
utility/Coordinates/Makefile
utility/Coordinates/cartesianDerivs.c
utility/Coordinates/l
utility/Coordinates/sgrid_Coordinates.c
utility/Coordinates/sgrid_Coordinates.h
utility/NumericUtils/
utility/NumericUtils/Makefile
utility/NumericUtils/NumericUtils.h
utility/NumericUtils/attenuation.c
utility/NumericUtils/brent.c
utility/NumericUtils/dpythag.c
utility/NumericUtils/dsvbksb.c
utility/NumericUtils/dsvdcmp.c
utility/NumericUtils/dsvdfit.c
utility/NumericUtils/f_to_minimize.c
utility/NumericUtils/fd_jacobian.c
utility/NumericUtils/fit.c
utility/NumericUtils/golden.c
utility/NumericUtils/integral.c
utility/NumericUtils/integral3D.c
utility/NumericUtils/linmin.c
utility/NumericUtils/lnsrch_double.c
utility/NumericUtils/lubksb.c
utility/NumericUtils/ludcmp.c
utility/NumericUtils/minimization.c
utility/NumericUtils/mnbrak.c
utility/NumericUtils/newton_lnsrch.c
utility/NumericUtils/nrutil.c
utility/NumericUtils/nrutil.h
utility/NumericUtils/powell.c
utility/NumericUtils/sgrid_NumericUtils.c
utility/NumericUtils/sgrid_NumericUtils.h
utility/Spectral/
utility/Spectral/.ispell.joe.bak
utility/Spectral/Makefile
utility/Spectral/Spec/
utility/Spectral/Spec/a.out*
utility/Spectral/Spec/adv.c
utility/Spectral/Spec/dnl16
utility/Spectral/Spec/dnl32
utility/Spectral/Spec/fd192
utility/Spectral/Spec/fd192d
utility/Spectral/Spec/fd384
utility/Spectral/Spec/fd48
utility/Spectral/Spec/fd96
utility/Spectral/Spec/l
utility/Spectral/Spec/l10
utility/Spectral/Spec/l16
utility/Spectral/Spec/l20
utility/Spectral/Spec/l24
utility/Spectral/Spec/l32
utility/Spectral/Spec/l48
utility/Spectral/Spec/n16
utility/Spectral/Spec/n24
utility/Spectral/Spec/n32
utility/Spectral/Spec/n48
utility/Spectral/Spec/nl16
utility/Spectral/Spec/nl24
utility/Spectral/Spec/nl32
utility/Spectral/Spec/nl48
utility/Spectral/Spec/nl8
utility/Spectral/Spec/sl24
utility/Spectral/Spec/sl32
utility/Spectral/Spec/sl48
utility/Spectral/Spec/slow_Cheb_trafos.c
utility/Spectral/Spec/snl24
utility/Spectral/Spec/snl48
utility/Spectral/Spec/wave.c
utility/Spectral/Spec/wave_compactified.c
utility/Spectral/Spec/wave_fd.c
utility/Spectral/Spec/wave_overlap.c
utility/Spectral/Spectral.h
utility/Spectral/a.out*
utility/Spectral/derivs.c
utility/Spectral/diffmatrices.c
utility/Spectral/explicit_Cheb_trafos.c
utility/Spectral/explicit_Cheb_trafos.c_bak
utility/Spectral/explicit_Cheb_trafos.c_flipsign
utility/Spectral/l
utility/Spectral/l2
utility/Spectral/sgrid_Spectral.c
utility/Spectral/sgrid_Spectral.h
utility/Spectral/slow_Cheb_trafos.c
utility/Spectral/slow_Cheb_trafos2.c
utility/Spectral/wave.c
utility/evolve/
utility/evolve/.tree.joe
utility/evolve/Makefile
utility/evolve/evolve.c
utility/evolve/evolve.h
utility/evolve/evolve_test.c
utility/evolve/evolve_test.par
utility/evolve/icn.c
utility/evolve/l
utility/evolve/rk.c
utility/evolve/sed.com
utility/evolve/sed.sh
utility/evolve/sgrid_evolve.c
utility/evolve/sgrid_evolve.h
utility/evolve/tags
utility/output/
utility/output/.ispell.joe.bak
utility/output/Makefile
utility/output/gnuplot_out2d.c
utility/output/output.c
utility/output/output.h
utility/output/sgrid_output.c
utility/output/sgrid_output.h
utility/output/xgraph_out1d.c

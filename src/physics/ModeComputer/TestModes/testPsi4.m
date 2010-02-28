<< ToFile.m
<< SpinWeightedSphericalHarmonics.m

Psi4[t_,p_] = Sum[ Sum[(l*100+m + I*m)*YY[-2,l,m,t,p], {m,-l,l}], {l,2,10}];


SphericalDFData["rpsi4.r3.l3", Re[Psi4[t,p]], t,p, 94,46];
SphericalDFData["ipsi4.r3.l3", Im[Psi4[t,p]], t,p, 94,46];

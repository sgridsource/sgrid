(* Get second derivs in domain 0 for Ansorgs coord in Eq. 24 *)

(* Set some assuptions *)

$Assumptions = A>=0 && A<=1  &&  B>=0 && B<=1  &&  phi>=0 && phi<2Pi

$Assumptions = $Assumptions && lep<0 $$ Ap[A_]>=0 && Ap[A_]<=1 && 
                                        Element[Ap[A_], Reals]

(* Ap = Sinh[A lep]/Sinh[lep] *)


$Assumptions = $Assumptions && sigp[B_,phi_]>0 
$Assumptions = $Assumptions && Element[sigp[B_,phi_], Reals] 

(* e.g.:  sigp[B_,phi_]= 0.8 + 0.15*Cos[B*2*Pi] + 0.15*Sin[phi] *)

(* Cp[B_,phi_] = Sqrt[ Tanh[ (sigp[B,phi] + I Pi B)/4] ] *)

(* Get all parts of Cp *)
Abstanh[x_,y_] = Sqrt[ (Cosh[2 x]-Cos[2 y])/(Cos[2 y] + Cosh[2 x]) ]
Argtanh[x_,y_] = ArcTan[Sin[2y]/Sinh[2x]]/2  (* <--is this right ??? *)

AbsCp[B_,phi_] = Sqrt[ Abstanh[sigp[B,phi]/4, Pi*B/4] ]
ArgCp[B_,phi_] = Argtanh[sigp[B,phi]/4, Pi*B/4]/2
 ReCp[B_,phi_] = AbsCp[B,phi] * Cos[ArgCp[B,phi]]
 ImCp[B_,phi_] = AbsCp[B,phi] * Sin[ArgCp[B,phi]]



(* Eq. 24 *)

X[A_,B_,Phi_] = (1-Ap[A])*(ReCp[B,phi] - B*ReCp[1,phi]) +
                B*Cos[(1-Ap[A])*ArgCp[1,phi]] + Ap[A]*(1-B)

R[A_,B_,Phi_] = (1-Ap[A])*(ImCp[B,phi] - B*ImCp[1,phi]) +
                B*Sin[(1-Ap[A])*ArgCp[1,phi]]


(* Define x,y,z in terms of X,R,phi *)
$Assumptions = X[A_,B_,Phi_]>=0  &&  R[A_,B_,Phi_]>=0 &&  phi>=0 && phi<2Pi

XRtoXRofABphi = {R->R[A,B,Phi],X->X[A,B,Phi]}
x[A_,B_,Phi_] = b (1/(R^2 + X^2)^2 + 1)(X^2 - R^2)/2 /. XRtoXRofABphi
y[A_,B_,Phi_] = b (1/(R^2 + X^2)^2 - 1) R X Cos[phi] /. XRtoXRofABphi
z[A_,B_,Phi_] = b (1/(R^2 + X^2)^2 - 1) R X Sin[phi] /. XRtoXRofABphi


(* Now take second derivs and replace Derivative *)

ReplaceDerivs= {Derivative[1, 0][sigp][B_,phi_]->dsigpdB[B,phi],
		Derivative[0, 1][sigp][B_,phi_]->dsigpdphi[B,phi],
		Derivative[2, 0][sigp][B_,phi_]->ddsigpddB[B,phi],
		Derivative[1, 1][sigp][B_,phi_]->ddsigpdBdphi[B,phi],
		Derivative[0, 2][sigp][B_,phi_]->ddsigpddphi[B,phi],
		Derivative[1][Ap][A_]->dApdA[A],
		Derivative[2][Ap][A_]->ddApddA[A]}

Replace1phi = { dsigpdB[1,phi_]     ->dsigpdB1[phi],
                dsigpdphi[1,phi_]   ->dsigpdphi1[phi],
               ddsigpddB[1,phi_]   ->ddsigpddB1[phi],
               ddsigpdBdphi[1,phi_]->ddsigpdBdphi1[phi],
               ddsigpddphi[1,phi_] ->ddsigpddphi1[phi]}

CFunctionFile = "ddxyzddABphi_AnsorgNS0.c"

Print["Writing to ", CFunctionFile, "\n"];
DeleteFile[CFunctionFile];
filepointer = OpenAppend[CFunctionFile];
pr[x_] := Module[{s},
  WriteString[filepointer, x];
];

pr["/* "<>CFunctionFile<>" */\n"];
pr["/* Copyright (C) 2007 Wolfgang Tichy, "<>
       ToString[Date[][[3]]]<>"."<>
       ToString[Date[][[2]]]<>"."<>
       ToString[Date[][[1]]]<>" */\n"];
pr["\n"];

pr["ddxyzddABphi[1][1][1] = "<>ToString[CForm[D[ x[A,B,phi], A,A] /. ReplaceDerivs /. Replace1phi]]<>";\n"];
pr["ddxyzddABphi[1][1][2] = "<>ToString[CForm[D[ x[A,B,phi], A,B] /. ReplaceDerivs /. Replace1phi]]<>";\n"];
pr["ddxyzddABphi[1][1][3] = "<>ToString[CForm[D[ x[A,B,phi], A,phi] /. ReplaceDerivs /. Replace1phi]]<>";\n"];
pr["ddxyzddABphi[1][2][2] = "<>ToString[CForm[D[ x[A,B,phi], B,B] /. ReplaceDerivs /. Replace1phi]]<>";\n"];
pr["ddxyzddABphi[1][2][3] = "<>ToString[CForm[D[ x[A,B,phi], B,phi] /. ReplaceDerivs /. Replace1phi]]<>";\n"];
pr["ddxyzddABphi[1][3][3] = "<>ToString[CForm[D[ x[A,B,phi], phi,phi] /. ReplaceDerivs /. Replace1phi]]<>";\n"];

pr["ddxyzddABphi[2][1][1] = "<>ToString[CForm[D[ y[A,B,phi], A,A] /. ReplaceDerivs /. Replace1phi]]<>";\n"];
pr["ddxyzddABphi[2][1][2] = "<>ToString[CForm[D[ y[A,B,phi], A,B] /. ReplaceDerivs /. Replace1phi]]<>";\n"];
pr["ddxyzddABphi[2][1][3] = "<>ToString[CForm[D[ y[A,B,phi], A,phi] /. ReplaceDerivs /. Replace1phi]]<>";\n"];
pr["ddxyzddABphi[2][2][2] = "<>ToString[CForm[D[ y[A,B,phi], B,B] /. ReplaceDerivs /. Replace1phi]]<>";\n"];
pr["ddxyzddABphi[2][2][3] = "<>ToString[CForm[D[ y[A,B,phi], B,phi] /. ReplaceDerivs /. Replace1phi]]<>";\n"];
pr["ddxyzddABphi[2][3][3] = "<>ToString[CForm[D[ y[A,B,phi], phi,phi] /. ReplaceDerivs /. Replace1phi]]<>";\n"];

pr["ddxyzddABphi[3][1][1] = "<>ToString[CForm[D[ z[A,B,phi], A,A] /. ReplaceDerivs /. Replace1phi]]<>";\n"];
pr["ddxyzddABphi[3][1][2] = "<>ToString[CForm[D[ z[A,B,phi], A,B] /. ReplaceDerivs /. Replace1phi]]<>";\n"];
pr["ddxyzddABphi[3][1][3] = "<>ToString[CForm[D[ z[A,B,phi], A,phi] /. ReplaceDerivs /. Replace1phi]]<>";\n"];
pr["ddxyzddABphi[3][2][2] = "<>ToString[CForm[D[ z[A,B,phi], B,B] /. ReplaceDerivs /. Replace1phi]]<>";\n"];
pr["ddxyzddABphi[3][2][3] = "<>ToString[CForm[D[ z[A,B,phi], B,phi] /. ReplaceDerivs /. Replace1phi]]<>";\n"];
pr["ddxyzddABphi[3][3][3] = "<>ToString[CForm[D[ z[A,B,phi], phi,phi] /. ReplaceDerivs /. Replace1phi]]<>";\n"];

Close[filepointer];

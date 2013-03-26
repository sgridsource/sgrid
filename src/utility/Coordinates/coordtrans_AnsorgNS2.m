(* Mathematica utility to make C-functions for coord transform *)

CoordTransfName = "AnsorgNS2";
FuncArgs = "tBox *box, int ind, double A, double B, double phi";

CCodeAtBeginning =
"double sigOfBphi = Coordinates_AnsorgNS_sigmam(box, ind, B, phi);
double sigOf1phi = Coordinates_AnsorgNS_sigmam(box, -1, 1.0, phi);
double d1sigOfBphi = Coordinates_AnsorgNS_dsigmam_dB(box, ind, B, phi);
double d2sigOfBphi = Coordinates_AnsorgNS_dsigmam_dphi(box, ind, B, phi);
double d1d1sigOfBphi = Coordinates_AnsorgNS_ddsigmam_dBdB(box, ind, B, phi);
double d1d2sigOfBphi = Coordinates_AnsorgNS_ddsigmam_dBdphi(box, ind, B, phi);
double d2d2sigOfBphi = Coordinates_AnsorgNS_ddsigmam_dphidphi(box, ind, B, phi);
double AbstanhWsBphi = Abstanh(0.25*sigOfBphi, 0.25*PI*B);
double d1AbstanhWsBphi = dAbstanhdx(0.25*sigOfBphi, 0.25*PI*B);
double d2AbstanhWsBphi = dAbstanhdy(0.25*sigOfBphi, 0.25*PI*B);
double d1d1AbstanhWsBphi = ddAbstanhdxdx(0.25*sigOfBphi, 0.25*PI*B);
double d1d2AbstanhWsBphi = ddAbstanhdxdy(0.25*sigOfBphi, 0.25*PI*B);
double d2d2AbstanhWsBphi = ddAbstanhdydy(0.25*sigOfBphi, 0.25*PI*B);
double ArgtanhWsBphi = Argtanh(0.25*sigOfBphi, 0.25*PI*B);
double d1ArgtanhWsBphi = dArgtanhdx(0.25*sigOfBphi, 0.25*PI*B);
double d2ArgtanhWsBphi = dArgtanhdy(0.25*sigOfBphi, 0.25*PI*B);
double d1d1ArgtanhWsBphi = ddArgtanhdxdx(0.25*sigOfBphi, 0.25*PI*B);
double d1d2ArgtanhWsBphi = ddArgtanhdxdy(0.25*sigOfBphi, 0.25*PI*B);
double d2d2ArgtanhWsBphi = ddArgtanhdydy(0.25*sigOfBphi, 0.25*PI*B);
double AbstanhWs1phi = Abstanh(0.25*sigOf1phi, 0.25*PI);
double ArgtanhWs1phi = Argtanh(0.25*sigOf1phi, 0.25*PI);
";

CCodeBeforeReturnValues =
"/* end var defs */\n";

(* some rules *)
POW2Rule    = Power[XXX_,2]  -> pow2[XXX];
POW2invRule = Power[XXX_,-2] -> pow2inv[XXX];
SomeRules = { POW2Rule, POW2invRule, POW2Rule, POW2invRule };

(* our coordtrafo is of the form: x = x(X) *)
x[1] = X;
x[2] = R;
x[3] = phi;
X[1] = A;
X[2] = B;
X[3] = phi;

(* how X,R,phi relate to our new coords *)
Ap = A;
(* we assume that Cp[1,phi] is independent of phi!!! *)
ToXYZrule[1] = x[1] -> (1-Ap)*(ReCpOfBphi[B,phi] - B*ReCpOf1phi) +
                        B*Cos[PIq*Ap + (1-Ap)*ArgCpOf1phi];
ToXYZrule[2] = x[2] -> (1-Ap)*(ImCpOfBphi[B,phi] - B*ImCpOf1phi) + 
                        B*Sin[PIq*Ap + (1-Ap)*ArgCpOf1phi];
ToXYZrule[3] = x[3] -> phi;


(* general rules to replace derivs of funcs *)
dXFofXrule = Derivative[1][F_][X_] :> ToExpression[
 ToString[d]<>ToString[F]];
dXFofXYrule = Derivative[1,0][F_][X_,Y_] :> ToExpression[
 ToString[d1]<>ToString[F]];
dYFofXYrule = Derivative[0,1][F_][X_,Y_] :> ToExpression[
 ToString[d2]<>ToString[F]];
dXdXFofXrule = Derivative[2][F_][X_] :> ToExpression[
 ToString[dd]<>ToString[F]];
dXdXFofXYrule = Derivative[2,0][F_][X_,Y_] :> ToExpression[
 ToString[d1d1]<>ToString[F]];
dXdYFofXYrule = Derivative[1,1][F_][X_,Y_] :> ToExpression[
 ToString[d1d2]<>ToString[F]];
dYdYFofXYrule = Derivative[0,2][F_][X_,Y_] :> ToExpression[
 ToString[d2d2]<>ToString[F]];
replaceDerivs =
 {dXFofXrule,dXFofXYrule,dYFofXYrule,
  dXdXFofXrule,dXdXFofXYrule,dXdYFofXYrule,dYdYFofXYrule};

(* rules to remove B and phi function arguments *)
removeBphiFuncArgs = { F_[phi] :> F , F_[B,phi] :> F };


removeArgsFromWsBphiFuncs = {AbstanhWsBphi[X_,Y_] :> AbstanhWsBphi,
                             ArgtanhWsBphi[X_,Y_] :> ArgtanhWsBphi,
                             AbstanhWs1phi[X_,Y_] :> AbstanhWs1phi,
                             ArgtanhWs1phi[X_,Y_] :> ArgtanhWs1phi}

(* lists of needed functions *)
FuncList = {
  AbsCpOfBphi[B,phi],
  ArgCpOfBphi[B,phi],
  ReCpOfBphi[B,phi],
  ImCpOfBphi[B,phi],
  AbsCpOf1phi, (*  AbsCpOf1phi[phi], <-- does not depend on phi *)
  ArgCpOf1phi, (*  ArgCpOf1phi[phi], *)
  ReCpOf1phi,  (*  ReCpOf1phi[phi],  *)
  ImCpOf1phi   (*  ImCpOf1phi[phi]   *)
};
FuncDefs = {
  Sqrt[ AbstanhWsBphi[sigOfBphi[B,phi]/4, PI*B/4] ],
  ArgtanhWsBphi[sigOfBphi[B,phi]/4, PI*B/4]/2,
  AbsCpOfBphi[B,phi] * Cos[ArgCpOfBphi[B,phi]],
  AbsCpOfBphi[B,phi] * Sin[ArgCpOfBphi[B,phi]],
  Sqrt[ AbstanhWs1phi[sigOf1phi[phi]/4, PI/4] ],
  ArgtanhWs1phi[sigOf1phi[phi]/4, PI/4]/2,
  AbsCpOf1phi[phi] * Cos[ArgCpOf1phi[phi]],
  AbsCpOf1phi[phi] * Sin[ArgCpOf1phi[phi]]
};

(* take derivs of fun lists *)
DFuncs[FuncList_]  := Flatten[{D[FuncList,B], D[FuncList,phi]}];
DDFuncs[FuncList_] := Flatten[{D[FuncList,B,B],
                      D[FuncList,B,phi], D[FuncList,phi,phi]}];



(* module to make specialVar array and their RHS *)
SetSpecialVars[specialVarList_,specialVarRHSList_,nstart_] := Module[{k,n},
 n=nstart;
 For[k = 1, k <= Length[specialVarList], k++,
   If[specialVarList[[k]]==0,
     n = n,
     n = n+1;
     specialVar[n] = specialVarList[[k]] /. removeBphiFuncArgs /.
                      removeArgsFromWsBphiFuncs /. SomeRules;
     specialVarRHS[n] = specialVarRHSList[[k]] /. removeBphiFuncArgs /.
                         removeArgsFromWsBphiFuncs /. SomeRules,
     n = n+1;
     specialVar[n] = specialVarList[[k]] /. removeBphiFuncArgs /.
                      removeArgsFromWsBphiFuncs /. SomeRules;
     specialVarRHS[n] = specialVarRHSList[[k]] /. removeBphiFuncArgs /.
                         removeArgsFromWsBphiFuncs /. SomeRules
   ];
 ];
 n
];


(* special Vars *)
NspecialVarsxOfX = SetSpecialVars[FuncList, FuncDefs, 0];
NspecialVarsdxdX  = SetSpecialVars[DFuncs[FuncList] /. replaceDerivs,
                                   DFuncs[FuncDefs]  /. replaceDerivs,
                                   NspecialVarsxOfX];
NspecialVarsddxdXdX = SetSpecialVars[DDFuncs[FuncList] /. replaceDerivs,
                                     DDFuncs[FuncDefs]  /. replaceDerivs,
                                     NspecialVarsdxdX];

(* special rules *)
NspecialRulesxOfX = 3;
NspecialRulesdxdX  = 3;
NspecialRulesddxdXdX = 3;
specialRule[1] = replaceDerivs;
specialRule[2] = removeBphiFuncArgs;
specialRule[3] = removeArgsFromWsBphiFuncs;


(* ***************************************************************** *)


CFunctionFile = "coordtrans_"<>CoordTransfName<>".c";

IncludeAndDefine[] := Module[{},

  pr["#include \"sgrid.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Log(x)     (log((double) (x)))\n"];
  pr["#define Sqrt(x)    (sqrt((double) (x)))\n"];
  pr["#define Sin(x)     (sin((double) (x)))\n"];
  pr["#define Cos(x)     (cos((double) (x)))\n"];
  pr["#define Csc(x)     (1.0/sin((double) (x)))\n"];
  pr["#define Cot(x)     (1.0/tan((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n"];
  pr["\n\n"];
];

GlobalVars[] := Module[{},

  pr["extern double (*Coordinates_AnsorgNS_sigmam)(tBox *box, int ind, double B, double phi);\n"];
  pr["extern double (*Coordinates_AnsorgNS_dsigmam_dB)(tBox *box, int ind, double B, double phi);\n"];
  pr["extern double (*Coordinates_AnsorgNS_dsigmam_dphi)(tBox *box, int ind, double B, double phi);\n"];
  pr["extern double (*Coordinates_AnsorgNS_ddsigmam_dBdB)(tBox *box, int ind, double B, double phi);\n"];
  pr["extern double (*Coordinates_AnsorgNS_ddsigmam_dBdphi)(tBox *box, int ind, double B, double phi);\n"];
  pr["extern double (*Coordinates_AnsorgNS_ddsigmam_dphidphi)(tBox *box, int ind, double B, double phi);\n"];
  pr["\n\n"];
];


(* ***************************************************************** *)

Print["Writing to ", CFunctionFile, "\n"];
DeleteFile[CFunctionFile];
filepointer = OpenAppend[CFunctionFile];
pr[x_] := Module[{s},
  WriteString[filepointer, x];
];

pr["/* "<>CFunctionFile<>" */\n"];
pr["/* Copyright (C) 2013 Wolfgang Tichy, "<>
       ToString[Date[][[3]]]<>"."<>
       ToString[Date[][[2]]]<>"."<>
       ToString[Date[][[1]]]<>" */\n"];
pr["/* Produced with Mathematica */\n"];
pr["\n"];


(* ***************************************************************** *)



(* module to print special Vars *)
prspecialVars[N_] := Module[{ss},
 For[k = 1, k <= N, k++,
   ss = StringForm["double `` = ``;\n",
                   CForm[specialVar[k]], CForm[specialVarRHS[k]]];
   pr[ss];
 ];
];

(* module for coord transform *)
prxOfX := Module[{xOfXYZ, fs, rs},
  pr["/* coord transforms */\n"];
  For[i = 1, i <= 3, i++,
    fs = StringForm["double ``_of``(``)\n",
                    x[i], CoordTransfName, FuncArgs];
    xOfXYZ = x[i] /. ToXYZrule[i];
    For[k = 1, k <= NspecialRulesxOfX, k++,
      xOfXYZ = xOfXYZ /. specialRule[k];
    ];
    xOfXYZ = xOfXYZ /. SomeRules;
    xOfXYZ = xOfXYZ /. Power[XXX_,-2] -> pow2inv[XXX];
    rs = StringForm["return ``;\n", CForm[xOfXYZ]];
    pr[fs];
    pr["{\n"];
    pr[CCodeAtBeginning];
    prspecialVars[NspecialVarsxOfX];
    pr[CCodeBeforeReturnValues];
    pr[rs];
    pr["}\n\n"];
    Print[fs];
  ];
  pr["\n\n\n"];
];

(* module for 1st derivs *)
prdxdX := Module[{dxdX, fs, rs},
  pr["/* 1st derivs */\n"];
  For[i = 1, i <= 3, i++,
    For[j = 1, j <= 3, j++,
      fs = StringForm["double d``_``_d``(``)\n",
                      x[i], CoordTransfName, X[j], FuncArgs];
      dxdX[i,j] = D[ x[i] /. ToXYZrule[i], X[j] ];
      For[k = 1, k <= NspecialRulesdxdX, k++,
        dxdX[i,j] = dxdX[i,j] /. specialRule[k];
      ];
      dxdX[i,j] = dxdX[i,j] /. SomeRules;
      dxdX[i,j] = Simplify[PowerExpand[Simplify[ dxdX[i,j] ]]];
      dxdX[i,j] = Simplify[ dxdX[i,j] ];
      dxdX[i,j] = N[ dxdX[i,j], 20 ] /. -1.*XXX_ :> -XXX;

      rs = StringForm["return ``;\n", CForm[dxdX[i,j]]];
      pr[fs];
      pr["{\n"];
      pr[CCodeAtBeginning];
      prspecialVars[NspecialVarsdxdX];
      pr[CCodeBeforeReturnValues];
      pr[rs];
      pr["}\n\n"];
      Print[fs];
    ];
  ];
  pr["\n\n\n"];
];

(* module for 2nd derivs *)
prddxdXdX := Module[{ddxdXdX, fs, rs},
  pr["/* 2nd derivs */\n"];
  For[i = 1, i <= 3, i++,
    For[j = 1, j <= 3, j++,
      For[k = j, k <= 3, k++,
        fs = StringForm["double dd``_``_d``d``(``)\n",
                        x[i], CoordTransfName, X[j], X[k], FuncArgs];
        ddxdXdX[i,j,k] = D[ x[i] /. ToXYZrule[i], X[j],X[k] ];
        For[n = 1, n <= NspecialRulesddxdXdX, n++,
          ddxdXdX[i,j,k] = ddxdXdX[i,j,k] /. specialRule[n];
        ];
        ddxdXdX[i,j,k] = ddxdXdX[i,j,k] /. SomeRules;
        ddxdXdX[i,j,k] = Simplify[PowerExpand[Simplify[ ddxdXdX[i,j,k] ]]];
        ddxdXdX[i,j,k] = Simplify[ ddxdXdX[i,j,k] ];
        ddxdXdX[i,j,k] = N[ ddxdXdX[i,j,k], 20 ] /. -1.*XXX_ :> -XXX;
        rs = StringForm["return ``;\n", CForm[ddxdXdX[i,j,k]]];
        pr[fs];
        pr["{\n"];
        pr[CCodeAtBeginning];
        prspecialVars[NspecialVarsddxdXdX];
        pr[CCodeBeforeReturnValues];
        pr[rs];
        pr["}\n\n"];
        Print[fs];
      ];
    ];
  ];
  pr["\n\n\n"];
];

(* module for one combned function with all derivs *)
prxyzdxdXddxdXdX := Module[{ddxdXdX, fs, rs},
  pr["/* func with trafo and 1st and 2nd derivs */\n"];
  fs = StringForm["void ``````_d``````_dd``````_of_``_``````(``, double ``````[4],double d``````[4][4],double dd``````[4][4][4])\n",
                  x[1],x[2],x[3], x[1],x[2],x[3], x[1],x[2],x[3], 
                  CoordTransfName, X[1],X[2],X[3], 
                  FuncArgs, x[1],x[2],x[3], x[1],x[2],x[3], x[1],x[2],x[3]];
  Print[fs];
  pr[fs];
  pr["{\n"];
  pr[CCodeAtBeginning];
  prspecialVars[NspecialVarsddxdXdX];
  pr[CCodeBeforeReturnValues];
  pr["\n"];
  pr["/* coord transforms */\n"];
  For[i = 1, i <= 3, i++,
    xOfXYZ = x[i] /. ToXYZrule[i];
    For[k = 1, k <= NspecialRulesxOfX, k++,
      xOfXYZ = xOfXYZ /. specialRule[k];
    ];
    xOfXYZ = xOfXYZ /. SomeRules;
    xOfXYZ = xOfXYZ /. Power[XXX_,-2] -> pow2inv[XXX];
    xOfXYZ = N[xOfXYZ] /. -1.*XXX_ :> -XXX;
    rs = StringForm["``````[``] = ``;\n", x[1],x[2],x[3], i, CForm[xOfXYZ]];
    pr[rs];
  ];
  pr["\n"];
  pr["/* 1st derivs */\n"];
  For[i = 1, i <= 3, i++,
    For[j = 1, j <= 3, j++,
      fs = StringForm["double d``_``_d``(``)\n",
                      x[i], CoordTransfName, X[j], FuncArgs];
      dxdX[i,j] = D[ x[i] /. ToXYZrule[i], X[j] ];
      For[k = 1, k <= NspecialRulesdxdX, k++,
        dxdX[i,j] = dxdX[i,j] /. specialRule[k];
      ];
      dxdX[i,j] = dxdX[i,j] /. SomeRules;
      dxdX[i,j] = Simplify[PowerExpand[Simplify[ dxdX[i,j] ]]];
      dxdX[i,j] = Simplify[ dxdX[i,j] ];
      dxdX[i,j] = N[ dxdX[i,j], 20 ] /. -1.*XXX_ :> -XXX;
      rs = StringForm["d``````[``][``] = ``;\n", 
                      x[1],x[2],x[3], i,j, CForm[dxdX[i,j]]];
      pr[rs];
    ];
  ];
  pr["\n"];
  pr["/* 2nd derivs */\n"];
  For[i = 1, i <= 3, i++,
    For[j = 1, j <= 3, j++,
      For[k = j, k <= 3, k++,
        ddxdXdX[i,j,k] = D[ x[i] /. ToXYZrule[i], X[j],X[k] ];
        For[n = 1, n <= NspecialRulesddxdXdX, n++,
          ddxdXdX[i,j,k] = ddxdXdX[i,j,k] /. specialRule[n];
        ];
        ddxdXdX[i,j,k] = ddxdXdX[i,j,k] /. SomeRules;
        ddxdXdX[i,j,k] = Simplify[PowerExpand[Simplify[ ddxdXdX[i,j,k] ]]];
        ddxdXdX[i,j,k] = Simplify[ ddxdXdX[i,j,k] ];
        ddxdXdX[i,j,k] = N[ ddxdXdX[i,j,k], 20 ] /. -1.*XXX_ :> -XXX;
        rs = StringForm["dd``````[``][``][``] = ``;\n", x[1],x[2],x[3], i,j,k,
                        CForm[ddxdXdX[i,j,k]]];
        pr[rs];
      ];
    ];
  ];
  pr["}\n\n"];
];




(***************************************************************************)
(* write functions *)
IncludeAndDefine[];
GlobalVars[];
(*
prxOfX;
prdxdX;
prddxdXdX;
*)
prxyzdxdXddxdXdX;

(***************************************************************************)
(* count operations *)

char = ReadList[CFunctionFile, Character];
nmul = Count[char, "*"];
ndiv = Count[char, "/"];
nsum = Count[char, "+"] + Count[char, "-"];
count = StringForm[
  "/* nvars = ``, n* = ``,  n/ = ``,  n+ = ``, n = ``, O = `` */\n", 
  Length[vars], nmul, ndiv, nsum, nmul+ndiv+nsum, 
  If [optimizeflag, 1, 0, 0]
];
Print[""];
Print[CFunctionFile];
Print[count];
pr["/* end of: "<>CFunctionFile<>" */\n"];
pr[count];

Close[filepointer];

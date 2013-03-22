(* Mathematica utility to make C-functions for coord transform *)

CoordTransfName = "AnsorgXRphi";
FuncArgs = "void *aux, int ind, double X, double R, double phi";

CCodeBeforeReturn =
"/* tBox *box = (tBox *) aux; */
if(XsqrPlusRsqr==0.0) XsqrPlusRsqr=dequaleps*dequaleps*dequaleps;\n";

(* some rules *)
POW2Rule    = Power[XXX_,2]  -> pow2[XXX];
POW2invRule = Power[XXX_,-2] -> pow2inv[XXX];
SomeRules = { POW2Rule, POW2invRule, POW2Rule, POW2invRule };

(* our coordtrafo is of the form: x = x(X) *)
x[1] = x;
x[2] = y;
x[3] = z;
X[1] = X;
X[2] = R;
X[3] = phi;

(* how x,y,z relate to our new coords *)
ToXYZrule[1] = x[1] -> b (1/(R^2 + X^2)^2 + 1)(X^2 - R^2)/2;
ToXYZrule[2] = x[2] -> b (1/(R^2 + X^2)^2 - 1) R X Cos[phi];
ToXYZrule[3] = x[3] -> b (1/(R^2 + X^2)^2 - 1) R X Sin[phi];

(* make it smarter *)
ToXYZrule[1] = x[1] -> b (1/(XRfunc[X,R])^2 + 1)(X^2 - R^2)/2;
ToXYZrule[2] = x[2] -> b (1/(XRfunc[X,R])^2 - 1) R X Cos[phi];
ToXYZrule[3] = x[3] -> b (1/(XRfunc[X,R])^2 - 1) R X Sin[phi];

(* special Vars *)
NspecialVarsxOfX = 1;
NspecialVarsdxdX  = 3;
NspecialVarsddxdXdX = 5;
specialVar[1] = XsqrPlusRsqr;
specialVar[2] = DXsqrPlusRsqrDX;
specialVar[3] = DXsqrPlusRsqrDR;
specialVar[4] = DDXsqrPlusRsqrDXDX;
specialVar[6] = DDXsqrPlusRsqrDXDR;
specialVar[5] = DDXsqrPlusRsqrDRDR;

(* special VarRHSs *)
specialVarRHS[1] = X^2 + R^2 /. SomeRules;
specialVarRHS[2] = D[X^2 + R^2, X] /. SomeRules;
specialVarRHS[3] = D[X^2 + R^2, R] /. SomeRules;
specialVarRHS[4] = D[X^2 + R^2, X,X] /. SomeRules;
specialVarRHS[6] = D[X^2 + R^2, X,R] /. SomeRules;
specialVarRHS[5] = D[X^2 + R^2, R,R] /. SomeRules;

(* special rules *)
NspecialRulesxOfX = 1;
NspecialRulesdxdX  = 3;
NspecialRulesddxdXdX = 6;
specialRule[1] = XRfunc[X,R] -> specialVar[1];
specialRule[2] = D[XRfunc[X,R],X] -> specialVar[2];
specialRule[3] = D[XRfunc[X,R],R] -> specialVar[3];
specialRule[4] = D[XRfunc[X,R],X,X] -> specialVar[4];
specialRule[6] = D[XRfunc[X,R],X,R] -> 0;
specialRule[5] = D[XRfunc[X,R],R,R] -> specialVar[5];


CFunctionFile = "coordtrans_"<>CoordTransfName<>".c";

IncludeAndDefine[] := Module[{},

  pr["#include \"sgrid.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Log(x)     (log((double) (x)))\n"];
  pr["#define Sin(x)     (sin((double) (x)))\n"];
  pr["#define Cos(x)     (cos((double) (x)))\n"];
  pr["#define Csc(x)     (1.0/sin((double) (x)))\n"];
  pr["#define Cot(x)     (1.0/tan((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];
  pr["\n\n\n"];
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
    prspecialVars[NspecialVarsxOfX];
    pr[CCodeBeforeReturn];
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
      dxdX[i,j] = N[ dxdX[i,j], 20 ] /. -1. XXX_ -> -XXX;
      rs = StringForm["return ``;\n", CForm[dxdX[i,j]]];
      pr[fs];
      pr["{\n"];
      prspecialVars[NspecialVarsdxdX];
      pr[CCodeBeforeReturn];
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
        ddxdXdX[i,j,k] = N[ ddxdXdX[i,j,k], 20 ] /. -1. XXX_ -> -XXX;
        rs = StringForm["return ``;\n", CForm[ddxdXdX[i,j,k]]];
        pr[fs];
        pr["{\n"];
        prspecialVars[NspecialVarsddxdXdX];
        pr[CCodeBeforeReturn];
        pr[rs];
        pr["}\n\n"];
        Print[fs];
      ];
    ];
  ];
  pr["\n\n\n"];
];

(***************************************************************************)
(* write functions *)
IncludeAndDefine[];
prxOfX;
prdxdX;
prddxdXdX;

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

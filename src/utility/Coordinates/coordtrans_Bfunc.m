(* Mathematica utility to make C-functions for coord transform *)

CoordTransfName = "Bfunc";
FuncArgs = "tBox *box, int ind, double Af, double Bf, double Phi";

CCodeAtBeginning =
"double ss = 0.5;
double aa = 2.0*ss-2.0;
double bb = 3.0-3.0*ss;
double es = 0.05;
";

CCodeBeforeReturnValues = "";

(* some rules *)
POW2Rule    = Power[XXX_,2]  -> pow2[XXX];
POW2invRule = Power[XXX_,-2] -> pow2inv[XXX];
SomeRules = { POW2Rule, POW2invRule, POW2Rule, POW2invRule };

(* our coordtrafo is of the form: x = x(X) *)
x[1] = A;
x[2] = B;
x[3] = phi;
X[1] = Af;
X[2] = Bf;
X[3] = Phi;

(* how x,y,z relate to our new coords *)
(* B=B[Bf] redistributes the points away or toward the axis *)
ToXYZrule[1] = x[1] -> Af;
ToXYZrule[2] = x[2] -> aa*Bf*Bf*Bf + bb*Bf*Bf + ss*Bf;
(*
ToXYZrule[2] = x[2] -> 1/2 + ( Sqrt[Bf+es] - Sqrt[1-Bf+es] )/
                             (2 (Sqrt[1+es]-Sqrt[es]) );
*)
ToXYZrule[3] = x[3] -> Phi;

(* special Vars *)
NspecialVarsxOfX = 0;
NspecialVarsdxdX  = 0;
NspecialVarsddxdXdX = 0;
specialVar[1] = XsqrPlusRsqr;

(* special VarRHSs *)
specialVarRHS[1] = X^2 + R^2 /. SomeRules;

(* special rules *)
NspecialRulesxOfX = 1;
NspecialRulesdxdX  = 1;
NspecialRulesddxdXdX = 1;
specialRule[1] = SomeRules;


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

  pr["extern int Coordinates_AnsorgNS_b_ParIndex;\n"];
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
      dxdX[i,j] = N[ dxdX[i,j], 20 ] /. -1. XXX_ -> -XXX;
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
        ddxdXdX[i,j,k] = N[ ddxdXdX[i,j,k], 20 ] /. -1. XXX_ -> -XXX;
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

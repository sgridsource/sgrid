(* Mathematica utility to make C-functions for coord transform *)

CoordTransfName = "SphericalDF";
FuncArgs = "void *aux, int ind, double r, double thm, double phi";

CCodeBeforeReturn =
"tBox *box = (tBox *) aux;
 double N = box->n2;
 double theta = thm + PI/N;\n";

X[1] = r;
X[2] = thm;
X[3] = phi;

ToXYZrule[1] = x -> r Cos[phi] Sin[theta];
ToXYZrule[2] = y -> r Sin[phi] Sin[theta];
ToXYZrule[3] = z -> r Cos[theta];

Toxyzrule[1] = r -> Sqrt[x^2 + y^2 + z^2];
Toxyzrule[2] = thm ->ArcCos[z/Sqrt[x^2 + y^2 + z^2]] - PI/N;
Toxyzrule[3] = phi -> ArcTan[y/x];

CFunctionFile = "coordinates_"<>CoordTransfName<>".c";

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
pr["/* Copyright (C) 2005 Wolfgang Tichy, "<>
       ToString[Date[][[3]]]<>"."<>
       ToString[Date[][[2]]]<>"."<>
       ToString[Date[][[1]]]<>" */\n"];
pr["/* Produced with Mathematica */\n"];
pr["\n"];


(* ***************************************************************** *)
x[1] = x;
x[2] = y;
x[3] = z;

ToXYZ = { ToXYZrule[1], ToXYZrule[2], ToXYZrule[3] };

(* module for coord transform *)
prxOfX := Module[{xOfXYZ, fs, rs},
  pr["/* coord transforms */\n"];
  For[i = 1, i <= 3, i++,
    xOfXYZ = x[i] /. ToXYZrule[i];
    fs = StringForm["double ``_of``(``)\n",
                    x[i], CoordTransfName, FuncArgs];
    rs = StringForm["return ``;\n", CForm[xOfXYZ]];
    pr[fs];
    pr["{\n"];
    pr[CCodeBeforeReturn];
    pr[rs];
    pr["}\n\n"];
    Print[fs];
  ];
  pr["\n\n\n"];
];

(* module for 1st derivs *)
prdXdx := Module[{dXdx, fs, rs},
  pr["/* 1st derivs */\n"];
  For[i = 1, i <= 3, i++,
    For[j = 1, j <= 3, j++,
      fs = StringForm["double d``_``_d``(``)\n",
                      X[i], CoordTransfName, x[j], FuncArgs];
      dXdx[i,j] = D[ X[i] /. Toxyzrule[i], x[j] ];
      dXdx[i,j] = dXdx[i,j] /. ToXYZ;
      dXdx[i,j] = dXdx[i,j] /. ToXYZ;
      dXdx[i,j] = dXdx[i,j] /. ToXYZ;
      dXdx[i,j] = Simplify[PowerExpand[Simplify[ dXdx[i,j] ]]];
      dXdx[i,j] = Simplify[ dXdx[i,j] ];
      dXdx[i,j] = N[ dXdx[i,j], 20 ] /. -1. xxx_ -> -xxx;
      rs = StringForm["return ``;\n", CForm[dXdx[i,j]]];
      pr[fs];
      pr["{\n"];
      pr[CCodeBeforeReturn];
      pr[rs];
      pr["}\n\n"];
      Print[fs];
    ];
  ];
  pr["\n\n\n"];
];

(* module for 2nd derivs *)
prddXdxdx := Module[{ddXdxdx, fs, rs},
  pr["/* 2nd derivs */\n"];
  For[i = 1, i <= 3, i++,
    For[j = 1, j <= 3, j++,
      For[k = j, k <= 3, k++,
        fs = StringForm["double dd``_``_d``d``(``)\n",
                        X[i], CoordTransfName, x[j], x[k], FuncArgs];
        ddXdxdx[i,j,k] = D[ X[i] /. Toxyzrule[i], x[j],x[k] ];
        ddXdxdx[i,j,k] = ddXdxdx[i,j,k] /. ToXYZ;
        ddXdxdx[i,j,k] = ddXdxdx[i,j,k] /. ToXYZ;
        ddXdxdx[i,j,k] = ddXdxdx[i,j,k] /. ToXYZ;
        ddXdxdx[i,j,k] = Simplify[PowerExpand[Simplify[ ddXdxdx[i,j,k] ]]];
        ddXdxdx[i,j,k] = Simplify[ ddXdxdx[i,j,k] ];
        ddXdxdx[i,j,k] = N[ ddXdxdx[i,j,k], 20 ] /. -1. xxx_ -> -xxx;
        rs = StringForm["return ``;\n", CForm[ddXdxdx[i,j,k]]];
        pr[fs];
        pr["{\n"];
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
prdXdx;
prddXdxdx;

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

(* BSSN_init.m 
   Wolfgang Tichy 4/2004, Bernd Bruegmann 10/98, 10/02  *)

(* initialize BSSN variables from ADM variables *)


(* variables *)
variables = {gb[a,b], K[a,b], psi, 
             gt[a,b], At[a,b], G[a], K, phi, alpha, alphaDensity,
             dgt[a,b,c]}


(* compute in this order *)
tocompute = {

  (* loop of all points *)
  Cinstruction == "forallpoints(box, ijk) {",

  (* conformal factors *)
  detgb == matrixdet[gb],
  p == detgb^(-1/3),
  p0 == Cal[usepsi,psi^(-4),1],

  (* basic rescaling *)
  gt[a,b] == p gb[a,b],
  Kt[a,b] == p p0 K[a,b], 
  alphaDensity == alpha (p p0)^(1.5 alphaDensityWeight),

  (* inverse BSSN metric *)
  detgtinv == 1/matrixdet[gt],
  gtinv[a,b] == detgtinv matrixinvdet[gt,a,b],

  (* additional variables *)
  phi == - log[p] / 4,
  K == gtinv[a,b] Kt[a,b],
  At[a,b] == Kt[a,b] - gt[a,b] K / 3,

  (* new loop *)
  Cinstruction == "} /* end 1st forallpoints loop */",
  Cinstruction == "FirstDerivsOf_Sab(box, Ind(\"BSSN_gxx\"), 
                   Ind(\"ADMvars_dgxxx\"));",
  Cinstruction == "forallpoints(box, ijk) {",

  (* inverse BSSN metric *)
  detgtinv == 1/matrixdet[gt],
  gtinv[a,b] == detgtinv matrixinvdet[gt,a,b],

  (* G *)
  Gtdown[a] == gtinv[b,c] dgt[a,b,c],
  G[a] == gtinv[a,b] Gtdown[b]
}


(* symmetries *)
g[a_,b_] := g[b,a] /; !OrderedQ[{a,b}]
gb[a_,b_] := gb[b,a] /; !OrderedQ[{a,b}]
gt[a_,b_] := gt[b,a] /; !OrderedQ[{a,b}]
gtinv[a_,b_] := gtinv[b,a] /; !OrderedQ[{a,b}]
K[a_,b_] := K[b,a] /; !OrderedQ[{a,b}]
Kt[a_,b_] := Kt[b,a] /; !OrderedQ[{a,b}]
At[a_,b_] := At[b,a] /; !OrderedQ[{a,b}]

dgt[a_,b_,c_] := dgt[b,a,c] /; !OrderedQ[{a,b}]




(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be in Cfunctionfile 
*)

CFunctionFile = "BSSN_init.c"

(* the head of the function *)
BeginCFunction[] := Module[{},

  pr["#include \"sgrid.h\"\n"];
  pr["#include \"BSSN.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Log(x)     log((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];

  pr["\n\n\n"];

  pr["void BSSN_init(tGrid *grid, int igb, int iK, int ipsi, int igt, int iAt, int iG, int itrK, int iphi, int i_alpha, int i_alphaDensity)\n"];
  pr["{\n"];
  pr["int bi;\n"];
  pr["\n"];

  pr["double alphaDensityWeight = Getd(\"BSSN_alphaDensityWeight\");\n"];
  pr["int usepsi = 1;\n"];

  pr["for(bi = 0; bi < grid->nboxes; bi++)\n"];
  pr["{\n"];
  pr["tBox *box = grid->box[bi];\n"];
  pr["int ijk;\n\n"];
  pr["\n"];
];

(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  prdecvar[{gb[a,b]},    "igb"];
  prdecvar[{K[a,b]},     "iK"];
  prdecvar[{psi},        "ipsi"];
  prdecvar[{gt[a,b]},    "igt"];
  prdecvar[{At[a,b]},    "iAt"];
  prdecvar[{G[a]},       "iG"];
  prdecvar[{K},          "itrK"];
  prdecvar[{phi},        "iphi"];
  prdecvar[{alpha},      "i_alpha"];
  prdecvar[{alphaDensity},"i_alphaDensity"];

  prdecvarname[{dgt[a,b,c]}, "ADMvars_dgxxx"];
  pr["\n"];
];    
(* auxillary variables are automatically inserted here *)

(* custom C-code which is placed directly after the variable declarations *)
InitializationCommands[] := Module[{},

  pr["/* Jetzt geht's los! */\n"];

];

(* the end or tail of the function (we need at least a }) *)
EndCFunction[] := Module[{},

  pr["} /* end of points */\n"];
  pr["} /* end of boxes */\n"];
  pr["\n\n"];
  pr["}  /* end of function */\n\n"];
];


(* to turn off optimization for development set optimizeflag = False *)
optimizeflag = True;


(************************************************************************)
(* now we are ready to go *)

<< "../../Math/MathToC/TensorEquations3dToC.m"

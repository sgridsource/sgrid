(* example4d.m
   Wolfgang Tichy 2.2.2005 *)


(* Variables *)
variables = { psi, g[a,b], K[a,b],  AR[a,b], ATT[a,b], tau,  U, LW[a,b] }

(* compute in this order *)
tocompute = {

  Cinstruction == "forallpoints(level,ijk) {",

  (* NOTE: g^{physical}_{ij} = Psi^4 g_{ij} *)

  (* Psi = psi + U *)
  Psi == psi + U,

  Psi2 == Psi * Psi,
  Psi4 == Psi2 * Psi2,
  ooPsi2 == 1.0/Psi2,

  g[a,b] == g[a,b] (Psi/psi)^4,

  psqrPN1 == mu^2 (4 + r12m)/r12m^2,

  W[a,b,c] == Del[a, g[b,c] ],
  V[a,b,c] == func[a, g[b,c] ],

  Bup[a,b] == AR[a,b] + ATT[a,b] + LW[a,b],
  Bdown[a,b] == g[a,c] g[b,d] Bup[c,d],
  K[a,b] == ooPsi2 (Bdown[a,b]) + Psi4 g[a,b] tau/3.0,

  Cinstruction == "}  /* end loop */"
}


(* symmetries *)
g[a_,b_]     := g[b,a]     /; !OrderedQ[{a,b}]
ginv[a_,b_]  := ginv[b,a]  /; !OrderedQ[{a,b}]

K[a_,b_]     := K[b,a]     /; !OrderedQ[{a,b}]
AR[a_,b_]    := AR[b,a]    /; !OrderedQ[{a,b}]
ATT[a_,b_]   := ATT[b,a]   /; !OrderedQ[{a,b}]
LW[a_,b_]    := LW[b,a]    /; !OrderedQ[{a,b}]
Bdown[a_,b_] := Bdown[b,a] /; !OrderedQ[{a,b}]
Bup[a_,b_]   := Bup[b,a]   /; !OrderedQ[{a,b}]


(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be in Cfunctionfile 
*)

CFunctionFile = "example4d.c"

(* the head of the function *)
BeginCFunction[] := Module[{},

  pr["#include \"bam.h\"\n"];
  pr["#include \"RealisticBBH.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Log(x)     log((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];

  pr["#define DEL(x,g)         ( partial((x), ijk, (g), ijk, x,y,z) )\n"];
  pr["#define DELDEL(x,y ,g)   ( partial((x), (y), (g), ijk, x,y,z) )\n"];
  pr["#define DEL2(x, g)       (  DELDEL((x), (x) ,(g)) )\n"];

  pr["\n\n\n"];

  pr["void example4d(tL *level, tVarList *vlMetric, tVarList *vlCurv)\n"];
  pr["{\n"];
  pr["int ijk=0;\n\n"];

];

(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  prdecvar[{psi},        "Ind(\"psi\")"]; 
  prdecvar[{g[a,b]},     "Ind(\"gxx\")"];
  prdecvar[{K[a,b]},     "Ind(\"Kxx\")"];
  prdecvar[{AR[a,b]},    "Ind(\"RealisticBBH_ARxx\")"];
  prdecvar[{ATT[a,b]},   "Ind(\"RealisticBBH_ATTxx\")"];
  prdecvar[{tau},        "Ind(\"RealisticBBH_tau\")"];
  prdecvar[{U},		 "Ind(\"RealisticBBH_U\")"];
  prdecvar[{LW[a,b]},    "Ind(\"RealisticBBH_LWxx\")"];
  pr["\n"];
];    
(* auxillary variables are automatically inserted here *)

(* custom C-code which is placed directly after the variable declarations *)
InitializationCommands[] := Module[{},

  pr["/* Jetzt geht's los! */\n"];
];


(* the end or tail of the function (we need at least a }) *)
EndCFunction[] := Module[{},

  pr["\n\n"];
  pr["}  /* end of function */\n\n"];
];


(* to turn off optimization for development set optimizeflag = False *)
optimizeflag = True;

(* use 4d tensors the default is 3 *)
TensorEquationsDim = 4;

(************************************************************************)
(* now we are ready to go *)

<< "TensorEquationsToC.m"

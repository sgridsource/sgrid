(* example.m
   Wolfgang Tichy 7.5.2005 *)


(* Variables *)
variables = { psi, g[a,b], K[a,b], LW[a,b] }

(* compute in this order *)
tocompute = {


  Cinstruction == "forallpoints(level,ijk) {",

  (* NOTE: g^{physical}_{ij} = Psi^4 g_{ij} *)

  (* Psi = psi + U *)
  Psi == psi + U,

  Psi2 == Psi * Psi,
  Psi4 == Psi2 * Psi2,

  g[a,b] == g[a,b] (Psi/psi)^4,

  psqrPN1 == mu^2 (4 + r12m)/r12m^2,


  W[a,b,c] == del[a, g[b,c] ],
  Z[a,b,c,d] == deldel[a,b, g[c,d] ],
  V[a,b,c] == func[a, g[b,c] ],

  Bup[a,b] == LW[a,b],
  Bdown[a,b] == g[a,c] g[b,d] Bup[c,d],
  K[a,b] == Bdown[a,b],

  Cinstruction == "}  /* end loop */"
}


(* symmetries *)
g[a_,b_]     := g[b,a]     /; !OrderedQ[{a,b}]

K[a_,b_]     := K[b,a]     /; !OrderedQ[{a,b}]
LW[a_,b_]    := LW[b,a]    /; !OrderedQ[{a,b}]
Bdown[a_,b_] := Bdown[b,a] /; !OrderedQ[{a,b}]
Bup[a_,b_]   := Bup[b,a]   /; !OrderedQ[{a,b}]


(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be in Cfunctionfile 
*)

CFunctionFile = "example.c"

(* the head of the function *)
BeginCFunction[] := Module[{},

  pr["#include \"sgrid.h\"\n"];
  pr["#include \"RealisticBBH.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Log(x)     log((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];

  pr["#define DEL(x,u)       dxu)\n"];
  pr["#define DELDEL(x,y,u)  dxyu)\n"];
  pr["#define DEL2(x,u)      dxxu)\n"];
(*
  pr["#define DEL(x,u)  partial(box, (x), (u), temp1, temp2, temp3)\n"];
  pr["#define DELDEL(x,y ,u) partial2((x), (y), (u), ijk, x,y,z)\n"];
  pr["#define DEL2(x, u)     DELDEL((x), (x) ,(u))\n"];
*)
  pr["\n\n\n"];

  pr["void example(tL *level, tVarList *vlMetric, tVarList *vlCurv)\n"];
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

(************************************************************************)
(* now we are ready to go *)

<< "TensorEquations3dToC.m"

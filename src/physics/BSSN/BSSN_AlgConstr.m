(* BSSN_AlgConstr.m 
   Wolfgang Tichy 9/2008, Bernd Bruegmann 4/2006  *)

(* compute trA and detg, and also enforce algebraic constraints *)


(* variables *)
variables = {g[a,b], A[a,b], trAstore, detgstore}

(* compute in this order *)
tocompute = {

  Cinstruction == "forallpoints(box, ijk) {",

  (* determinant of metric *)
  detg == matrixdet[g],

  (* enforce detg = 1 *)
  Cif == normalizedetg,
    aux == detg^(-oothree),
    g[a,b] == aux g[a,b],
    detg == matrixdet[g],
  Cif == end,

  (* trace of extrinsic curvature using new metric if rescaled *)
  detginv == 1/detg,
  trA == detginv matrixinvdet[g,a,b] A[a,b],

  (* enforce trA = 0 *)
  Cif == subtractA,
    aux == -trA oothree,
    A[a,b] == A[a,b] + aux g[a,b],
  Cif == end,

  (* store result if storage is enabled *)
  Cinstruction == "/* no storage at the moment */"
  (* Cinstruction == "if (detgstore) detgstore[ccc] = detg;\n",
     Cinstruction == "if (trAstore)  trAstore[ccc] = trA;\n" *)
}



(* symmetries *)
g[a_,b_]    := g[b,a] /; !OrderedQ[{a,b}]
A[a_,b_]    := A[b,a] /; !OrderedQ[{a,b}]
ginv[a_,b_] := ginv[b,a] /; !OrderedQ[{a,b}]
Ainv[a_,b_] := Ainv[b,a] /; !OrderedQ[{a,b}]


(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be in Cfunctionfile 
*)

CFunctionFile = "BSSN_AlgConstr.c"

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

  pr["void BSSN_AlgConstr(tVarList *ucur)\n"];
  pr["{\n"];
  pr["tGrid *grid = ucur->grid;\n"];
  pr["int bi;\n"];
  pr["\n"];

  pr["const int subtractA     = Getv(\"BSSN_subtractA\", \"yes\");\n"];
  pr["const int normalizedetg = Getv(\"BSSN_normalizedetg\", \"yes\");\n"];
  pr["const double oothree = 1.0/3.0;\n"];

  pr["\n"];

  pr["forallboxes(grid,bi)\n"];
  pr["{\n"];
  pr["tBox *box = grid->box[bi];\n"];
  pr["int ijk;\n\n"];
  pr["\n"];
  pr["\n"];
];

(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  (* prdecvl[{g[a,b], A[a,b], trAstore, detgstore}, "ucur"]; *)
  prdecvl[{g[a,b], A[a,b]}, "ucur"];

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

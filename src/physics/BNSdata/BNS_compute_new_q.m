(* BNS_compute_new_q.m 
   Wolfgang Tichy  2/2008       *)

(* compute the updated q *)


(* variables *)
variables = {Psi, B[a], alphaP, Sigma, dSigma[a], q, wB[a], x, y}

constvariables = {OmegaCrossR[a]}

(* compute in this order *)
tocompute = {

  Cinstruction == "FirstDerivsOf_S(box,index_BNSdata_Sigma,
                                   Ind(\"BNSdata_Sigmax\"));",

  (* loop of all points *)
  Cinstruction == "forallpoints(box, ijk) {",

  (* Omega \times r term *)
  OmegaCrossR1 == - Omega y,
  OmegaCrossR2 == + Omega (x-xCM),
  OmegaCrossR3 == 0,

  (* shift in rotating frame (in the inertial frame beta^i = B^i) *)
  beta[a] == B[a] + OmegaCrossR[a],

  (* some abbreviations *)
  alpha  == alphaP/Psi,
  alpha2 == alpha*alpha,
  Psi2   == Psi*Psi,
  Psi4   == Psi2*Psi2,

  (**************)
  (* corotation *)
  (**************)
  Cif == ( ((bi<=1 || bi==5) && corot1) || ((bi>=2 && bi<=4) && corot2) ),
    (* vR[a] is zero for corotation *)
    vR[a] == 0,

    (* compute u^0 in rotating frame *)
    oouzerosqr == alpha2 - Psi4 delta[b,c] (beta[b] + vR[b]) (beta[c] + vR[c]),
    Cif == (oouzerosqr==0),
      oouzerosqr == -1,
    Cif == end,
    Cif == (oouzerosqr<0),
      uzero == -1, (* -Sqrt[-1/oouzerosqr], *)
    Cif == else,
      uzero == Sqrt[1/oouzerosqr],
    Cif == end,

    (* killing vec xi^i in rotating frame, xi^a = (xi^0, xi^i) *)
    xi[a] == 0,

    (* constant F *)
    F == -(alpha2 - Psi4 delta[b,c] (beta[b] + xi[b]) (beta[c] + vR[c]) ) *
	uzero,

    (* set q *)
    Cif == ((bi==0) || (bi==1) || (bi==5)),  
      q == (C1/F - 1.0)/(n+1.0),
    Cif == else,
      q == (C2/F - 1.0)/(n+1.0),
    Cif == end,

  (****************)
  (* general case *)
  (****************)
  Cif == else,
    (* which const to use ? *)
    Cif == ((bi<=1) || (bi==5)),  
      CC == C1,
    Cif == else,
      CC == C2,
    Cif == end,

    Psim4 == 1/Psi4,
    Psim6 == Psim4/Psi2,
    DSigmaUp[a] == Psim4 dSigma[a],
    w[a] == Psim6 wB[a],
    twoalpha2wdSigma == 2 alpha2 w[c] dSigma[c],
    betadSigmaMinusCC == beta[c] dSigma[c] - CC, 
    bb == betadSigmaMinusCC^2 - twoalpha2wdSigma,
    L2 == (bb + sqrt[bb*bb - twoalpha2wdSigma^2])/(2 alpha2),
    h == sqrt[L2 - dSigma[a] DSigmaUp[a]],

    (* h == (n+1) q + 1, *)
    q == (h - 1.0)/(n + 1.0),
  Cif == end,

  Cinstruction == "} /* end of points loop */\n"
}


(* symmetries *)
g[a_,b_] :=  g[b,a] /; !OrderedQ[{a,b}]

(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be in Cfunctionfile 
*)

CFunctionFile = "BNS_compute_new_q.c"

(* the head of the function *)
BeginCFunction[] := Module[{},

  pr["#include \"sgrid.h\"\n"];
  pr["#include \"BNSdata.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Sqrt(x)    (sqrt((double) (x)))\n"];
  pr["#define Log(x)     (log((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];

  pr["\n\n\n"];

  pr["void BNS_compute_new_q(tGrid *grid)\n"];
  pr["{\n"];

  pr["int corot1 = Getv(\"BNSdata_rotationstate1\",\"corotation\");\n"];
  pr["int corot2 = Getv(\"BNSdata_rotationstate2\",\"corotation\");\n"];
  pr["double n = Getd(\"BNSdata_n\");\n"];
  pr["double C1 = Getd(\"BNSdata_C1\");\n"];
  pr["double C2 = Getd(\"BNSdata_C2\");\n"];
  pr["double kappa = Getd(\"BNSdata_kappa\");\n"];
  pr["double Omega = Getd(\"BNSdata_Omega\");\n"];
  pr["double xCM = Getd(\"BNSdata_x_CM\");\n"];
  pr["\n"];

  pr["int bi;\n"];
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

  prdecvarname[{x},      "x"];
  prdecvarname[{y},      "y"];

  (* prdecvarname[{g[a,b]},    "gxx"]; prdecvarname[{K[a,b]}, "Kxx"]; *)
  prdecvarname[{Psi},       "BNSdata_Psi"];
  prdecvarname[{alphaP},    "BNSdata_alphaP"];
  prdecvarname[{B[a]},      "BNSdata_Bx"];
  prdecvarname[{q},         "BNSdata_q"];
  prdecvarname[{wB[a]},    "BNSdata_wBx"];
  prdecvarname[{Sigma},     "BNSdata_Sigma"];
  prdecvarname[{dSigma[a]}, "BNSdata_Sigmax"];

  pr["\n"];
];    
(* auxillary variables are automatically inserted here *)

(* custom C-code which is placed directly after the variable declarations *)
InitializationCommands[] := Module[{},

  pr["/* Jetzt geht's los! */\n"];

];

(* the end or tail of the function (we need at least a }) *)
EndCFunction[] := Module[{},

  pr["} /* end of boxes */\n"];
  pr["\n\n"];
  pr["}  /* end of function */\n\n"];
];


(* to turn off optimization for development set optimizeflag = False *)
optimizeflag = True;

(* use 3d tensors the default is 3 *)
TensorEquationsDim = 3;

(************************************************************************)
(* now we are ready to go *)

<< "../../Math/MathToC/TensorEquationsToC.m"

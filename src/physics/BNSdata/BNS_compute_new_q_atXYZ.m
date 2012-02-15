(* BNS_compute_new_q_atXYZ.m 
   Wolfgang Tichy  2/2008       *)

(* compute the updated q *)


(* variables *)
variables = {BNSPsi, BNSB[a], BNSalphaP, BNSSigma, BNSdSigma[a], BNSwB[a],
             BNSqgold, temp4}

constvariables = {Psi, B[a], alphaP, Sigma, dSigma[a], wB[a], x,y, qgold,
                  OmegaCrossR[a]}

(* compute in this order *)
tocompute = {

  (* deriv of Sigma *)
  Cinstruction == "FirstDerivsOf_S(box,index_BNSdata_Sigma,
                                   Ind(\"BNSdata_Sigmax\"));",
  Cif == (bi==1),
    Cinstruction == "\n
    copy_Var_at_i0_from_Box1_Box2(grid, Ind(\"BNSdata_Sigmax\"), 0,1);\n
    copy_Var_at_i0_from_Box1_Box2(grid, Ind(\"BNSdata_Sigmay\"), 0,1);\n
    copy_Var_at_i0_from_Box1_Box2(grid, Ind(\"BNSdata_Sigmaz\"), 0,1);\n
    copy_Var_at_i0_from_Box1_Box2(grid, Ind(\"BNSdata_wBx\"), 0,1);\n
    copy_Var_at_i0_from_Box1_Box2(grid, Ind(\"BNSdata_wBy\"), 0,1);\n
    copy_Var_at_i0_from_Box1_Box2(grid, Ind(\"BNSdata_wBz\"), 0,1);",
  Cif == end,
  Cif == (bi==2),
    Cinstruction == "\n
    copy_Var_at_i0_from_Box1_Box2(grid, Ind(\"BNSdata_Sigmax\"), 3,2);\n
    copy_Var_at_i0_from_Box1_Box2(grid, Ind(\"BNSdata_Sigmay\"), 3,2);\n
    copy_Var_at_i0_from_Box1_Box2(grid, Ind(\"BNSdata_Sigmaz\"), 3,2);\n
    copy_Var_at_i0_from_Box1_Box2(grid, Ind(\"BNSdata_wBx\"), 3,2);\n
    copy_Var_at_i0_from_Box1_Box2(grid, Ind(\"BNSdata_wBy\"), 3,2);\n
    copy_Var_at_i0_from_Box1_Box2(grid, Ind(\"BNSdata_wBz\"), 3,2);",
  Cif == end,

  (* x,y *)
  Cinstruction == "if(box->x_of_X[1] != NULL) {",
  Cinstruction == "x = box->x_of_X[1]((void *) box, -1, X,Y,Z);",
  Cinstruction == "y = box->x_of_X[2]((void *) box, -1, X,Y,Z);",
  Cinstruction == "} else {",
  Cinstruction == "x = X;",
  Cinstruction == "y = Y;",
  Cinstruction == "}",

  (* Psi, B[a], alphaP, Sigma, dSigma[a], wB[a] by interpolation *)
  Cinstruction == "spec_Coeffs(box, BNSPsi, temp4);",
  Cinstruction == "Psi = spec_interpolate(box, temp4, X,Y,Z);",

  Cinstruction == "spec_Coeffs(box, BNSB1, temp4);",
  Cinstruction == "B1 = spec_interpolate(box, temp4, X,Y,Z);",

  Cinstruction == "spec_Coeffs(box, BNSB2, temp4);",
  Cinstruction == "B2 = spec_interpolate(box, temp4, X,Y,Z);",

  Cinstruction == "spec_Coeffs(box, BNSB3, temp4);",
  Cinstruction == "B3 = spec_interpolate(box, temp4, X,Y,Z);",

  Cinstruction == "spec_Coeffs(box, BNSalphaP, temp4);",
  Cinstruction == "alphaP = spec_interpolate(box, temp4, X,Y,Z);",

  Cinstruction == "spec_Coeffs(box, BNSSigma, temp4);",
  Cinstruction == "Sigma = spec_interpolate(box, temp4, X,Y,Z);",

  Cinstruction == "spec_Coeffs(box, BNSdSigma1, temp4);",
  Cinstruction == "dSigma1 = spec_interpolate(box, temp4, X,Y,Z);",

  Cinstruction == "spec_Coeffs(box, BNSdSigma2, temp4);",
  Cinstruction == "dSigma2 = spec_interpolate(box, temp4, X,Y,Z);",

  Cinstruction == "spec_Coeffs(box, BNSdSigma3, temp4);",
  Cinstruction == "dSigma3 = spec_interpolate(box, temp4, X,Y,Z);",

  Cinstruction == "spec_Coeffs(box, BNSwB1, temp4);",
  Cinstruction == "wB1 = spec_interpolate(box, temp4, X,Y,Z);",

  Cinstruction == "spec_Coeffs(box, BNSwB2, temp4);",
  Cinstruction == "wB2 = spec_interpolate(box, temp4, X,Y,Z);",

  Cinstruction == "spec_Coeffs(box, BNSwB3, temp4);",
  Cinstruction == "wB3 = spec_interpolate(box, temp4, X,Y,Z);",


  (**************)
  (* compute q: *)

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
    Psim2 == Psim4 Psi2,
    Psim6 == Psim4 Psim2,
    DSigmaUp[a] == Psim4 dSigma[a],
    w[a] == Psim6 wB[a],
    wBDown[a] == wB[a],
    wDown[a] == Psim2 wBDown[a],
    Cif == qFromFields,
      twoalpha2wdSigmapw == 2 alpha2 w[c] (dSigma[c]+wDown[c]),
      betadSigmaMinusCC == beta[c] dSigma[c] - CC, 
      bb == betadSigmaMinusCC^2 + twoalpha2wdSigmapw,
      L2 == (bb + Sqrt[Abs[bb*bb + twoalpha2wdSigmapw^2]])/(2 alpha2),
      h == Sqrt[Abs[L2 - (dSigma[a]+wDown[a]) (DSigmaUp[a]+w[a])]],
    Cif == else,
      Cinstruction == "spec_Coeffs(box, BNSqgold, temp4);",
      Cinstruction == "qgold = spec_interpolate(box, temp4, X,Y,Z);",
      (* h == (n+1) q + 1, *)
      hOLD == (n+1.0) qgold + 1.0,
      hOLD2 == hOLD*hOLD,
      LOLD2 == hOLD2 + (wDown[c] + dSigma[c]) (w[c] + DSigmaUp[c]),
      uzerosqr == LOLD2/(alpha2 hOLD2),
      uzero == sqrt[uzerosqr],
      vR[a] == (w[a] + DSigmaUp[a])/(uzero*hOLD) - beta[a],
      h == -( CC + dSigma[a] vR[a] ) uzero,
    Cif == end,

    (* h == (n+1) q + 1, *)
    q == (h - 1.0)/(n + 1.0),
  Cif == end,

  Cinstruction == "/* end of computation */\n"
}


(* symmetries *)
g[a_,b_] :=  g[b,a] /; !OrderedQ[{a,b}]

(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be in Cfunctionfile 
*)

CFunctionFile = "BNS_compute_new_q_atXYZ.c"

(* the head of the function *)
BeginCFunction[] := Module[{},

  pr["#include \"sgrid.h\"\n"];
  pr["#include \"BNSdata.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Sqrt(x)    (sqrt((double) (x)))\n"];
  pr["#define Abs(x)     (fabs((double) (x)))\n"];
  pr["#define Log(x)     (log((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];

  pr["\n\n\n"];

  pr["double BNS_compute_new_q_atXYZ(tGrid *grid, int bi, double X, double Y, double Z)\n"];
  pr["{\n"];

  pr["int corot1 = Getv(\"BNSdata_rotationstate1\",\"corotation\");\n"];
  pr["int corot2 = Getv(\"BNSdata_rotationstate2\",\"corotation\");\n"];
  pr["int qFromFields = Getv(\"BNSdata_new_q\",\"FromFields\");\n"];
  pr["double n = Getd(\"BNSdata_n\");\n"];
  pr["double C1 = Getd(\"BNSdata_C1\");\n"];
  pr["double C2 = Getd(\"BNSdata_C2\");\n"];
  pr["double kappa = Getd(\"BNSdata_kappa\");\n"];
  pr["double Omega = Getd(\"BNSdata_Omega\");\n"];
  pr["double xCM = Getd(\"BNSdata_x_CM\");\n"];
  pr["\n"];

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

  (* prdecvarname[{g[a,b]},    "gxx"]; prdecvarname[{K[a,b]}, "Kxx"]; *)
  prdecvarname[{BNSPsi},       "BNSdata_Psi"];
  prdecvarname[{BNSalphaP},    "BNSdata_alphaP"];
  prdecvarname[{BNSB[a]},      "BNSdata_Bx"];
  prdecvarname[{BNSq},         "BNSdata_q"];
  prdecvarname[{BNSwB[a]},     "BNSdata_wBx"];
  prdecvarname[{BNSSigma},     "BNSdata_Sigma"];
  prdecvarname[{BNSdSigma[a]}, "BNSdata_Sigmax"];
  prdecvarname[{BNSqgold},     "BNSdata_qgold"];
  prdecvarname[{temp4},        "BNSdata_temp4"];

  pr["double Psi, B1,B2,B3, alphaP;\n"];
  pr["double Sigma, dSigma1,dSigma2,dSigma3, wB1,wB2,wB3, x,y, qgold;\n"];

  pr["\n"];
];    
(* auxillary variables are automatically inserted here *)

(* custom C-code which is placed directly after the variable declarations *)
InitializationCommands[] := Module[{},

  pr["/* Jetzt geht's los! */\n"];

];

(* the end or tail of the function (we need at least a }) *)
EndCFunction[] := Module[{},

  pr["return q;\n"];
  pr["}  /* end of function */\n\n"];
];


(* to turn off optimization for development set optimizeflag = False *)
optimizeflag = True;

(* use 3d tensors the default is 3 *)
TensorEquationsDim = 3;

(************************************************************************)
(* now we are ready to go *)

<< "../../Math/MathToC/TensorEquationsToC.m"

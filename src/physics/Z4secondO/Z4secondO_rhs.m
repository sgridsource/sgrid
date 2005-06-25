(* Z4secondO_rhs.m 
   Bernd Bruegmann 10/98, 10/02
   Wolfgang Tichy  4/2004       *)

(* compute right hand side of Z4secondO equations *)


(* variables *)
variables = {g[a,b],  K[a,b],  Theta,  Z[a],   alpha,  beta[a],  B[a], 
             ng[a,b], nK[a,b], nTheta, nZ[a], nalpha, nbeta[a], nB[a],
             pg[a,b], pK[a,b], pTheta, pZ[a], palpha, pbeta[a], pB[a],
	     psi, dpop[a], ddpop[a,b], K0,
             dg[a,b,c], ddg[a,b,c,d], dK[a,b,c], dTheta[a],  
             dZ[a,b], dalp[a], ddalp[a,b], dbeta[a,b], ddbeta[a,b,c] }

(* compute in this order *)
tocompute = {

  Cif == useDD,

  Cinstruction == "allDerivsOf_Sab(box, index_g11, \
                    Ind(\"ADMvars_dgxxx\"), Ind(\"ADMvars_ddgxxxx\"));",
  Cinstruction == "FirstDerivsOf_Sab(box, index_K11, Ind(\"ADMvars_dKxxx\"));",

  Cinstruction == "FirstDerivsOf_Sa(box, index_Z1, Ind(\"Z4secondO_dZxx\"));",
  Cinstruction == "FirstDerivsOf_S(box, index_Theta, Ind(\"Z4secondO_dThetax\"));",

  Cinstruction == "allDerivsOf_S(box, index_alpha, \
                    Ind(\"Z4secondO_dalpx\"), Ind(\"Z4secondO_ddalpxx\"));",
  Cinstruction == "allDerivsOf_Sa(box, index_beta1, \
                    Ind(\"Z4secondO_dbetaxx\"), Ind(\"Z4secondO_ddbetaxxx\"));",

  Cif == else,

  Cinstruction == "FirstAndSecondDerivsOf_Sab(box, index_g11, \
                    Ind(\"ADMvars_dgxxx\"), Ind(\"ADMvars_ddgxxxx\"));",
  Cinstruction == "FirstDerivsOf_Sab(box, index_K11, Ind(\"ADMvars_dKxxx\"));",

  Cinstruction == "FirstDerivsOf_Sa(box, index_Z1, Ind(\"Z4secondO_dZxx\"));",
  Cinstruction == "FirstDerivsOf_S(box, index_Theta, Ind(\"Z4secondO_dThetax\"));",

  Cinstruction == "FirstAndSecondDerivsOf_S(box, index_alpha, \
                    Ind(\"Z4secondO_dalpx\"), Ind(\"Z4secondO_ddalpxx\"));",
  Cinstruction == "FirstAndSecondDerivsOf_Sa(box, index_beta1, \
                    Ind(\"Z4secondO_dbetaxx\"), Ind(\"Z4secondO_ddbetaxxx\"));",

  Cif == end,

  (* loop of all points *)
  Cinstruction == "forallpoints(box, ijk) {",

  (* former shift adv terms *)
  betadg[a,b] == beta[d] dg[a,b,d],
  betadK[a,b] == beta[d] dK[a,b,d],
  betadTheta  == beta[d] dTheta[d],
  betadZ[a]   == beta[d] dZ[a,d],
  betadalp    == beta[d] dalp[d],

  (* make transition to physical metric *)
  (*
     f == psi^4,
     df[a] == 4 f dpop[a],
     ddf[a,b] == 4 f (ddpop[a,b] + 3 dpop[a] dpop[b]), 
     ddg[c,d,a,b] == f ddg[c,d,a,b] + ddf[a,b] g[c,d] +
                       df[a] dg[c,d,b] + df[b] dg[c,d,a],  
     dg[a,b,c] == f dg[a,b,c] + g[a,b] df[c],
  *)

  (* inverse physical metric *)
  detginv == 1/( matrixdet[g]),
  ginv[a,b] == detginv matrixinvdet[g,a,b],

  (* connection of physical metric *) 
  gammado[c,a,b] == 1/2 (dg[b,c,a] + dg[a,c,b] - dg[a,b,c]),
  gamma[c,a,b] == ginv[c,d] gammado[d,a,b], 

  (* curvature of physical metric *)
  R[a,b] == ginv[c,d] ( 1/2 (
    - ddg[a,b,c,d] - ddg[c,d,a,b] + 
      ddg[b,d,a,c] + ddg[a,d,b,c]) +
      gamma[e,a,c] gammado[e,b,d] - 
      gamma[e,a,b] gammado[e,c,d]),
  R == ginv[a,b] R[a,b],

  (* extrinsic curvature terms*)
  KK[a,b] == ginv[c,d] K[a,c] K[d,b],
  trK  == ginv[a,b] K[a,b],
  trKK == ginv[a,b] KK[a,b],

  (* Hamiltonian constraint *)
  hamil == R + trK trK - trKK,

  (* subtract hamil from R  (In Bernds bssn: RtoRminusHfactor=1) *)
  R == R - RtoRminusHfactor * hamil,

  (* derivatives of inverse metric *)
  dginv[b,c,a] == - ginv[b,d] ginv[c,e] dg[d,e,a],
  
  (* derivative of the trace of the extrinsic curvature *)
  dK[a] == dginv[b,c,a] K[b,c] + ginv[b,c] dK[b,c,a],

  (* covariant derivatives of the lapse *)
  cddalp[a,b] == ddalp[a,b] - gamma[c,a,b] dalp[c],

  (* covariant derivatives of Z *)
  cdZ[a,b] == dZ[a,b] - gamma[c,a,b] Z[c],

  (* covariant derivatives of the extrinsic curvature *)
  cdK[a,b,c] == dK[a,b,c] - gamma[d,a,b] K[d,c] - gamma[d,a,c] K[b,d],

  (* Lie derivs *)
  lieg[a,b] == betadg[a,b] + dbeta[c,a] g[b,c] + dbeta[c,b] g[a,c],
  lieK[a,b] == betadK[a,b] + dbeta[c,a] K[b,c] + dbeta[c,b] K[a,c],
  lieTheta  == betadTheta,
  lieZ[a]   == betadZ[a] + dbeta[c,a] Z[c],

  (* right hand sides *)

  rg[a,b] == - 2 alpha K[a,b] + lieg[a,b],

  rK[a,b] == - cddalp[a,b] + alpha (R[a,b] + cdZ[a,b] + cdZ[b,a] -
             2 KK[a,b] + K[a,b](trK - 2 Theta)) + lieK[a,b],

  rTheta == alpha/2 ( R + 2 cdZ[a,b] ginv[a,b] + 
                      trK(trK - 2 Theta) - trKK ) - 
            ginv[a,b] Z[a] dalp[b] + lieTheta,

  rZ[a] == alpha ( cdK[b,a,c] ginv[b,c] - dK[a] + dTheta[a] -
           2 K[a,b] Z[c] ginv[b,c] ) - Theta dalp[a] + lieZ[a],

  (* gauge conditions, avoid unnecessary if statements, use flags *)
  liealpha == betadalp,
  ralpha0 == -alpha (trK - subtractK0 K0 - lapseharmonicm Theta) *
               psi^lapsepsipower,
  ralpha  == nonconstantlapse * ( withshift liealpha +
              (oploglapse lapseharmonicf  + harmoniclapse alpha) ralpha0 ),

  rbeta[a] == 0,
  rB[a] == 0,
(*
  ralpha0 == (6 withshift liephi - 
              alpha (trK - subtractK0 K0)) psi^lapsepsipower,
  ralpha  == ralpha0 * nonconstantlapse *
             (oploglapse lapseharmonicf +
              oploglapse2 8/3/(3-alpha) + harmoniclapse alpha),

   betaF == shiftgammacoeff alpha^shiftalphapower / psi^shiftpsipower,
   rbeta[a] == evolveshift (gamma2factor + gamma0factor betaF) B[a],
   rB[a]    == evolveshift ((gamma0factor + gamma2factor betaF) rG[a] -  
                             shiftdriver B[a]),
*)

  (* result *) 
  Cif == addlinear, 
    ng[a,b] == pg[a,b] + dt rg[a,b],
    nK[a,b] == pK[a,b] + dt rK[a,b],
    nZ[a] == pZ[a] + dt rZ[a],
    nTheta == pTheta + dt rTheta,
    nalpha == palpha + dt ralpha,
    nbeta[a] == pbeta[a] + dt rbeta[a],
    nB[a] == pB[a] + dt rB[a], 

  Cif == else,
    ng[a,b] == rg[a,b],
    nK[a,b] == rK[a,b],
    nZ[a] == rZ[a],
    nTheta == rTheta,
    nalpha == ralpha,
    nbeta[a] == rbeta[a],
    nB[a] == rB[a],
  Cif == end
}


(* symmetries *)
 g[a_,b_] :=  g[b,a] /; !OrderedQ[{a,b}]
pg[a_,b_] := pg[b,a] /; !OrderedQ[{a,b}]
ng[a_,b_] := ng[b,a] /; !OrderedQ[{a,b}]
rg[a_,b_] := rg[b,a] /; !OrderedQ[{a,b}]

 K[a_,b_] :=  K[b,a] /; !OrderedQ[{a,b}]
pK[a_,b_] := pK[b,a] /; !OrderedQ[{a,b}]
nK[a_,b_] := nK[b,a] /; !OrderedQ[{a,b}]
rK[a_,b_] := rK[b,a] /; !OrderedQ[{a,b}]

ginv[a_,b_] := ginv[b,a] /; !OrderedQ[{a,b}]
K[a_,b_] := K[b,a] /; !OrderedQ[{a,b}]
KK[a_,b_] := KK[b,a] /; !OrderedQ[{a,b}]
R[a_,b_] := R[b,a] /; !OrderedQ[{a,b}]
Rphi[a_,b_] := Rphi[b,a] /; !OrderedQ[{a,b}]
lieg[a_,b_] := lieg[b,a] /; !OrderedQ[{a,b}]
lieK[a_,b_] := lieK[b,a] /; !OrderedQ[{a,b}]

ddf[a_,b_] := ddf[b,a] /; !OrderedQ[{a,b}]
cddf[a_,b_] := cddf[b,a] /; !OrderedQ[{a,b}]
ddalp[a_,b_] := ddalp[b,a] /; !OrderedQ[{a,b}]
cddalp[a_,b_] := cddalp[b,a] /; !OrderedQ[{a,b}]
ddpop[a_,b_] := ddpop[b,a] /; !OrderedQ[{a,b}]

gammado[c_,a_,b_] := gammado[c,b,a] /; !OrderedQ[{a,b}]
gamma[c_,a_,b_] := gamma[c,b,a] /; !OrderedQ[{a,b}]
dginv[a_,b_,c_] := dginv[b,a,c] /; !OrderedQ[{a,b}]

(* additional derivs on grid *)
dg[a_,b_,c_]     := dg[b,a,c]    /; !OrderedQ[{a,b}]
ddg[a_,b_,c_,d_] := ddg[b,a,c,d] /; !OrderedQ[{a,b}]
ddg[a_,b_,c_,d_] := ddg[a,b,d,c] /; !OrderedQ[{c,d}]
dK[a_,b_,c_]     := dK[b,a,c]    /; !OrderedQ[{a,b}]
ddalp[a_,b_]      := ddalp[b,a]    /; !OrderedQ[{a,b}]
ddbeta[c_,a_,b_]  := ddbeta[c,b,a] /; !OrderedQ[{a,b}]


(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be in Cfunctionfile 
*)

CFunctionFile = "Z4secondO_rhs.c"

(* the head of the function *)
BeginCFunction[] := Module[{},

  pr["#include \"sgrid.h\"\n"];
  pr["#include \"Z4secondO.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Log(x)     log((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];

  pr["\n\n\n"];

  pr["void Z4secondO_rhs(tVarList *unew, tVarList *upre, double dt, tVarList *ucur)\n"];
  pr["{\n"];
  pr["tGrid *grid = ucur->grid;\n"];
  pr["int bi;\n"];
  pr["\n"];

  pr["int addlinear = (dt != 0.0l);\n"];
  pr["/* int usepsi = 1; */\n"];

  pr["int useDD               = Getv(\"Z4secondO_useDD\", \"yes\");\n"];
  pr["double RtoRminusHfactor = Getd(\"Z4secondO_RtoRminusHfactor\");\n"];

  pr["int nonconstantlapse    =!Getv(\"Z4secondO_lapse\", \"constant\");\n"];
  pr["int oploglapse          = Getv(\"Z4secondO_lapse\", \"1+log\");\n"];
  pr["int withshift           = Getv(\"Z4secondO_lapse\", \"withshift\");\n"];
  pr["int harmoniclapse       = Getv(\"Z4secondO_lapse\", \"harmonic\");\n"];
  pr["int subtractK0          = Getv(\"Z4secondO_subtractK0\", \"yes\");\n"];

  pr["double gamma0factor    = Getv(\"Z4secondO_shift\", \"gamma0\");\n"];
  pr["double gamma2factor    = Getv(\"Z4secondO_shift\", \"gamma2\");\n"];
  pr["double shift_st        = Getd(\"Z4secondO_shift_stop_time\");\n"];
  pr["double evolveshift     = Cal( (grid->time>=shift_st)*(shift_st>=0), 0.0, 1.0);\n"];

  pr["double lapsepsipower   = Getd(\"Z4secondO_lapsepsipower\");\n"];
  pr["double lapseharmonicf  = Getd(\"Z4secondO_lapseharmonicf\");\n"];
  pr["double lapseharmonicm  = Getd(\"Z4secondO_lapseharmonicm\");\n"];
  pr["double shiftpsipower   = Getd(\"Z4secondO_shiftpsipower\");\n"];
  pr["double shiftalphapower = Getd(\"Z4secondO_shiftalphapower\");\n"];
  pr["double shiftgammacoeff = Getd(\"Z4secondO_shiftgammacoeff\");\n"];
  pr["double shiftdriver     = Getd(\"Z4secondO_shiftdriver\");\n"];
  pr["\n"];

  pr["for(bi = 0; bi < grid->nboxes; bi++)\n"];
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

  prdecvl[{K0}, "K_initial"];
  prdecvl[{psi, dpop[a], ddpop[a,b]}, "psiandderivs"];
  prdecvl[{ g[a,b],  K[a,b],  Theta,  Z[a],  alpha,  beta[a],  B[a]}, "ucur"];
  prdecvl[{ng[a,b], nK[a,b], nTheta, nZ[a], nalpha, nbeta[a], nB[a]}, "unew"];
  prdecvl[{pg[a,b], pK[a,b], pTheta, pZ[a], palpha, pbeta[a], pB[a]}, "upre"];

  prdecvlindices[{g[a,b], K[a,b], Theta, Z[a], alpha, beta[a], B[a]}, "ucur"];

  prdecvarname[{dg[a,b,c]},     "ADMvars_dgxxx"];
  prdecvarname[{ddg[a,b,c,d]},  "ADMvars_ddgxxxx"];
  prdecvarname[{dK[a,b,c]},     "ADMvars_dKxxx"];
  prdecvarname[{dTheta[a]},     "Z4secondO_dThetax"];
  prdecvarname[{dZ[a,b]},       "Z4secondO_dZxx"];
  prdecvarname[{dalp[a]},       "Z4secondO_dalpx"];
  prdecvarname[{ddalp[a,b]},    "Z4secondO_ddalpxx"];
  prdecvarname[{dbeta[a,b]},    "Z4secondO_dbetaxx"];
  prdecvarname[{ddbeta[a,b,c]}, "Z4secondO_ddbetaxxx"];
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

<< "../../Math/MathToC/TensorEquationsToC.m"

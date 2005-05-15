(* BSSN_rhs.m 
   Bernd Bruegmann 10/98, 10/02
   Wolfgang Tichy  4/2004       *)

(* compute right hand side of BSSN equations *)


(* variables *)
variables = {g[a,b],  A[a,b],  G[a],  K,  phi,   alpha,  beta[a],  B[a], 
             ng[a,b], nA[a,b], nG[a], nK, nphi, nalpha, nbeta[a], nB[a],
             pg[a,b], pA[a,b], pG[a], pK, pphi, palpha, pbeta[a], pB[a],
	      alphaDensity,
	     nalphaDensity,
	     palphaDensity,
	     psi, dpop[a], ddpop[a,b], K0,
             dgt[a,b,c], ddgt[a,b,c,d], dAt[a,b,c], dphi[a], ddphi[a,b], 
             dGt[a,b], dK[a], dalp[a], ddalp[a,b], dbeta[a,b], ddbeta[a,b,c] }

(* compute in this order *)
tocompute = {

  Cinstruction == "FirstAndSecondDerivsOf_Sab(box, index_g11, \
                    Ind(\"ADMvars_dgxxx\"), Ind(\"ADMvars_ddgxxxx\"));",
  Cinstruction == "FirstDerivsOf_Sab(box, index_A11, Ind(\"ADMvars_dKxxx\"));",

  Cinstruction == "FirstAndSecondDerivsOf_S(box, index_phi, \
                    Ind(\"BSSN_dphix\"), Ind(\"BSSN_ddphixx\"));",
  Cinstruction == "FirstDerivsOf_Sa(box, index_G1, Ind(\"BSSN_dGxx\"));",
  Cinstruction == "FirstDerivsOf_S(box, index_K, Ind(\"BSSN_dKx\"));",

  Cinstruction == "FirstAndSecondDerivsOf_S(box, index_alpha, \
                    Ind(\"BSSN_dalpx\"), Ind(\"BSSN_ddalpxx\"));",
  Cinstruction == "FirstAndSecondDerivsOf_Sa(box, index_beta1, \
                    Ind(\"BSSN_dbetaxx\"), Ind(\"BSSN_ddbetaxxx\"));",

  (* loop of all points *)
  Cinstruction == "forallpoints(box, ijk) {",

  (* transfer partial derivatives into Bernd's vars *)
  df[a] == dphi[a],
  ddf[a,b] == ddphi[a,b],
  da[a] == dalp[a],
  dda[a,b] == ddalp[a,b],
  db[a,b] == dbeta[b,a],
  ddb[a,b,c] == ddbeta[c,a,b],
  delg[c,a,b] == dgt[a,b,c],
  deldelg[a,b,c,d] == ddgt[c,d,a,b],
  delG[a,b] == dGt[b,a],
  (* dK[a] == dK[a], *)
  dA[a,b,c] == dAt[b,c,a],

  (* former shift adv terms *)
  betadf == beta[d] dphi[d],
  betadg[a,b] == beta[d] dgt[a,b,d],
  betadA[a,b] == beta[d] dAt[a,b,d],
  betadK == beta[d] dK[d],
  betadG[a] == beta[d] dGt[a,d],

  (* inverse conformal metric *)
  detginv == 1/matrixdet[g],
  ginv[a,b] == detginv matrixinvdet[g,a,b],

  (* derivatives of conformal metric *)
  dginv[a,b,c] == - ginv[b,d] ginv[c,e] delg[a,d,e],

  (* connection of conformal metric *) 
  gammado[c,a,b] == 1/2 (delg[a,b,c] + delg[b,a,c] - delg[c,a,b]),
  gamma[c,a,b] == ginv[c,d] gammado[d,a,b], 
  Gfromg[a] == ginv[b,c] gamma[a,b,c],

  (* curvature of conformal metric *)
  R[a,b] == ginv[c,d] ( -1/2 deldelg[c,d,a,b] +
                    gamma[e,c,a] gammado[b,e,d] +
                    gamma[e,c,b] gammado[a,e,d] +
                    gamma[e,a,d] gammado[e,c,b]) +
               1/2 (g[c,a] delG[b,c] + g[c,b] delG[a,c] +
               Gfromg[c] (gammado[a,b,c] + gammado[b,a,c])),

  (* conformal factor and its derivatives *)
  f == phi,
  logpsi == 0,
  Cif == usepsi,
    logpsi == log[psi],
    f == f + logpsi,
    df[a] == df[a] + dpop[a],
    ddf[a,b] == ddf[a,b] - dpop[a] dpop[b] + ddpop[a,b],
    betadf == betadf + beta[a] dpop[a],
  Cif == end,
  cddf[a,b] == ddf[a,b] - gamma[c,a,b] df[c],
  trcddf == ginv[a,b] cddf[a,b],
  psim4 == exp[-4 f],

  (* curvature contribution from conformal factor *)
  Rphi[a,b] == 4 df[a] df[b] - 2 cddf[a,b] - 2 g[a,b] trcddf -
               4 g[a,b] ginv[c,d] df[c] df[d],

  (* derivatives of the lapse *)
  cdda[a,b] == dda[a,b] - gamma[c,a,b] da[c] -
               2 (df[a] da[b] + df[b] da[a] - g[a,b] ginv[c,d] df[c] da[d]),
  trcdda == psim4 (ginv[a,b] cdda[a,b]),

  (* set K zero *)
  K == K forceKzerofactor,
  dK[a] == dK[a] forceKzerofactor,

  (* tf part of conformal extrinsic curvature *)
  AA[a,b] == ginv[c,d] A[a,c] A[d,b],
  AA == ginv[a,b] AA[a,b],
  Ainv[a,b] == ginv[a,c] ginv[b,d] A[c,d],
  divAinv[a] == 2/3 ginv[a,b] dK[b] - 6 Ainv[a,b] df[b] -
                  gamma[a,b,c] Ainv[b,c],

  (* Ricci scalar *)
  R == AA - 2 K K / 3,

  (* the shift terms *)
  divbeta == delta[a,b] db[a,b],
  totdivbeta == 2/3 divbeta,
  ootddivbeta[a] == 1/3 delta[b,c] ddb[a,b,c],
  
  lieg[a,b] == betadg[a,b] - g[a,b] totdivbeta +
                 db[a,c] g[b,c] + db[b,c] g[a,c],

  lieA[a,b] == betadA[a,b] - A[a,b] totdivbeta +
                 db[a,c] A[b,c] + db[b,c] A[a,c],

  lieK == betadK,

  (* liephi is NOT a real Lie deriv.! It is simply:
     $ liephi = \nabla_i^{phys. metric} \beta^i /6 
              = \partial_i \beta^i /6 + \beta^i \partial_i BSSNphi $ *)
  liephi == betadf + divbeta/6,

  pseudolieG[a] == ginv[b,c] ddb[b,c,a] + ginv[a,b] ootddivbeta[b] -
                   Gfromg[b] db[b,a] + Gfromg[a] totdivbeta + 
                   betadG[a],

  (* right hand sides *)
  rg[a,b] == - 2 alpha A[a,b] + lieg[a,b],

  rA[a,b] == psim4 ( - cdda[a,b] + alpha (R[a,b] + Rphi[a,b])) -
             1/3 g[a,b] (- trcdda + alpha R) +
             alpha (K A[a,b] - 2 AA[a,b]) + lieA[a,b],

  rG[a] == - 2 Ainv[a,b] da[b] - 2 alpha divAinv[a] + pseudolieG[a], 

  rK == - trcdda + alpha (AA + K K/3) + lieK, 

  rphi == - alpha K/6 + liephi,


  (* gauge conditions, avoid unnecessary if statements, use flags *)
  ralpha0 == (6 oplogwithshift liephi - 
              alpha (K - subtractK0 K0)) psi^lapsepsipower,
  ralpha  == ralpha0 * nonconstantlapse *
             (oploglapse lapseharmonicf +
              oploglapse2 8/3/(3-alpha) + harmoniclapse alpha),

  Cif == densitizedLapse,
    (* exponentials of total conformal factor used in densitized lapse *)
    E6alphaDensityWeightf == exp[6 alphaDensityWeight f],
    ooE6alphaDensityWeightf == 1.0/E6alphaDensityWeightf,

    ralphaDensity  == ooE6alphaDensityWeightf ralpha - 
                      6 alphaDensityWeight alphaDensity rphi,

    Cif == densitizedoplogWithoutShift,
       ralphaDensity  == alphaDensity * 
           ( (lapseharmonicf/alpha) psi^lapsepsipower - alphaDensityWeight ) *
           (-alpha (K - subtractK0 K0)),
    Cif == end,

    ralphaDensity  == ralphaDensity * nonconstantlapse,
  Cif == end,

  betaF == shiftgammacoeff alpha^shiftalphapower / psi^shiftpsipower,
  rbeta[a] == evolveshift (gamma2factor + gamma0factor betaF) B[a],
  rB[a]    == evolveshift ((gamma0factor + gamma2factor betaF) rG[a] -  
                            shiftdriver B[a]),

  (* result *) 
  Cif == addlinear, 
    ng[a,b] == pg[a,b] + dt rg[a,b],
    nA[a,b] == pA[a,b] + dt rA[a,b],
    nG[a] == pG[a] + dt rG[a],
    nK == forceKzerofactor (pK + dt rK),
    nphi == pphi + dt rphi,
    Cif == densitizedLapse,
      nalphaDensity == palphaDensity + dt ralphaDensity,
      (* compute new alpha from new alphaDensity *)
      nf == nphi + logpsi,
      nalpha == nalphaDensity exp[6 alphaDensityWeight nf],
    Cif == else,
      nalpha == palpha + dt ralpha,
      (* compute new alphaDensity from new alpha *)
      (* nf == nphi + logpsi,
         nalphaDensity == nalpha exp[-6 alphaDensityWeight nf], *)
    Cif == end,
    nbeta[a] == pbeta[a] + dt rbeta[a],
    nB[a] == pB[a] + dt rB[a], 

    (* normalize detg and subtract trace from As *)
    detnginv == 1/matrixdet[ng],
    Cif == subtractA,
      traceA == detnginv matrixinvdet[ng,a,b] nA[a,b],
      aux == -traceA/3,
      nA[a,b] == nA[a,b] + aux ng[a,b],
    Cif == end,
    Cif == normalizedetg,
      aux == detnginv^(1/3),
      ng[a,b] == aux ng[a,b], 
    Cif == end,

  Cif == else,
    ng[a,b] == rg[a,b],
    nA[a,b] == rA[a,b],
    nG[a] == rG[a],
    nK == forceKzerofactor rK,
    nphi == rphi,
    Cif == densitizedLapse,
      nalphaDensity == ralphaDensity,
      nalpha == E6alphaDensityWeightf (ralphaDensity + 
                6 alphaDensityWeight rphi alphaDensity),
      Cinstruction == "errorexit(\"error in nalphaDensity: we should use the new E6alphaDensityWeightf, but we are using the old one instead.\");",
    Cif == else,
      nalpha == ralpha,
      (* nalphaDensity == ooE6alphaDensityWeightf (ralpha -
                       6 alphaDensityWeight rphi alpha),     *)
      (* error in nalphaDensity: we should use the new ooE6alphaDensityWeightf, but we are using the old one instead. *)
    Cif == end,
    nbeta[a] == rbeta[a],
    nB[a] == rB[a],
  Cif == end

  (* constraints and other diagnostic functions *)
(*
  Cif == constraints,
    hamil == psim4 ginv[a,b] (R[a,b] + Rphi[a,b]) - AA + 2/3 K K,
    Cif == divAfromD,
      Adiag[a,b] == ginv[a,c] A[c,b],
      diffeo[a] == ginv[a,b] ginv[c,d] dA[c,d,b] +
                   Adiag[b,c] dginv[b,c,a] - Adiag[a,b] G[b] - divAinv[a],
    Cif == else,
      diffeo[a] == divAinv[a] - 2/3 ginv[a,b] dK[b] + 6 Ainv[a,b] df[b] +
                   gamma[a,b,c] Ainv[b,c],
    Cif == end,
  Cif == end
*)
}


(* symmetries *)
 g[a_,b_] :=  g[b,a] /; !OrderedQ[{a,b}]
pg[a_,b_] := pg[b,a] /; !OrderedQ[{a,b}]
ng[a_,b_] := ng[b,a] /; !OrderedQ[{a,b}]
rg[a_,b_] := rg[b,a] /; !OrderedQ[{a,b}]

 A[a_,b_] :=  A[b,a] /; !OrderedQ[{a,b}]
pA[a_,b_] := pA[b,a] /; !OrderedQ[{a,b}]
nA[a_,b_] := nA[b,a] /; !OrderedQ[{a,b}]
rA[a_,b_] := rA[b,a] /; !OrderedQ[{a,b}]

ginv[a_,b_] := ginv[b,a] /; !OrderedQ[{a,b}]
Ainv[a_,b_] := Ainv[b,a] /; !OrderedQ[{a,b}]
Kinv[a_,b_] := Kinv[b,a] /; !OrderedQ[{a,b}]
K[a_,b_] := K[b,a] /; !OrderedQ[{a,b}]
AA[a_,b_] := AA[b,a] /; !OrderedQ[{a,b}]
R[a_,b_] := R[b,a] /; !OrderedQ[{a,b}]
Rphi[a_,b_] := Rphi[b,a] /; !OrderedQ[{a,b}]
lieg[a_,b_] := lieg[b,a] /; !OrderedQ[{a,b}]
lieA[a_,b_] := lieA[b,a] /; !OrderedQ[{a,b}]
result[a_,b_] := result[b,a] /; !OrderedQ[{a,b}]

ddf[a_,b_] := ddf[b,a] /; !OrderedQ[{a,b}]
cddf[a_,b_] := cddf[b,a] /; !OrderedQ[{a,b}]
dda[a_,b_] := dda[b,a] /; !OrderedQ[{a,b}]
cdda[a_,b_] := cdda[b,a] /; !OrderedQ[{a,b}]
ddpop[a_,b_] := ddpop[b,a] /; !OrderedQ[{a,b}]
ddb[a_,b_,c_] := ddb[b,a,c] /; !OrderedQ[{a,b}]

deldel[a_,b_,c_] := deldel[b,a,c] /; !OrderedQ[{a,b}]

delg[c_,a_,b_] := delg[c,b,a] /; !OrderedQ[{a,b}]
delgb[c_,a_,b_] := delgb[c,b,a] /; !OrderedQ[{a,b}]
dginv[c_,a_,b_] := dginv[c,b,a] /; !OrderedQ[{a,b}]
dA[c_,a_,b_] := dA[c,b,a] /; !OrderedQ[{a,b}]

deldelg[a_,b_,c_,d_] := deldelg[b,a,c,d] /; !OrderedQ[{a,b}]
deldelg[a_,b_,c_,d_] := deldelg[a,b,d,c] /; !OrderedQ[{c,d}]

codel[a_, (x_ /; NumberQ[x])] = 0
codelK[c_,a_,b_] := codelK[c,b,a] /; !OrderedQ[{a,b}]

(* additional derivs on grid *)
dgt[a_,b_,c_]     := dgt[b,a,c]    /; !OrderedQ[{a,b}]
ddgt[a_,b_,c_,d_] := ddgt[b,a,c,d] /; !OrderedQ[{a,b}]
ddgt[a_,b_,c_,d_] := ddgt[a,b,d,c] /; !OrderedQ[{c,d}]
dAt[a_,b_,c_]     := dAt[b,a,c]    /; !OrderedQ[{a,b}]
ddphi[a_,b_]      := ddphi[b,a]    /; !OrderedQ[{a,b}]
ddalp[a_,b_]      := ddalp[b,a]    /; !OrderedQ[{a,b}]
ddbeta[c_,a_,b_]  := ddbeta[c,b,a] /; !OrderedQ[{a,b}]


(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be in Cfunctionfile 
*)

CFunctionFile = "BSSN_rhs.c"

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

  pr["void BSSN_rhs(tVarList *unew, tVarList *upre, double dt, tVarList *ucur)\n"];
  pr["{\n"];
  pr["tGrid *grid = ucur->grid;\n"];
  pr["int bi;\n"];
  pr["\n"];

  pr["int addlinear = (dt != 0.0l);\n"];
  pr["int usepsi = 1;\n"];

  pr["double forceKzerofactor = Getv(\"BSSN_forceKzero\", \"no\");\n"];
  pr["int subtractA           = Getv(\"BSSN_subtractA\", \"yes\");\n"];
  pr["int normalizedetg       = Getv(\"BSSN_normalizedetg\", \"yes\");\n"];
  pr["int nonconstantlapse    =!Getv(\"BSSN_lapse\", \"constant\");\n"];
  pr["int oploglapse          = Getv(\"BSSN_lapse\", \"1+log\");\n"];
  pr["int oploglapse2         = Getv(\"BSSN_lapse\", \"1+log2\");\n"];
  pr["int oplogwithshift      = Getv(\"BSSN_lapse\", \"withshift\");\n"];
  pr["int harmoniclapse       = Getv(\"BSSN_lapse\", \"harmonic\");\n"];
  pr["int subtractK0          = Getv(\"BSSN_subtractK0\", \"yes\");\n"];

  pr["int densitizedLapse = !Getv(\"BSSN_densitizedLapse\", \"no\");\n"];
  pr["int densitizedoplogWithoutShift = Getv(\"BSSN_densitizedLapse\", \"1+log_withoutShift\");\n"];
  pr["double alphaDensityWeight = Getd(\"BSSN_alphaDensityWeight\");\n"];

  pr["double gamma0factor    = Getv(\"BSSN_shift\", \"gamma0\");\n"];
  pr["double gamma2factor    = Getv(\"BSSN_shift\", \"gamma2\");\n"];
  pr["double shift_st        = Getd(\"BSSN_shift_stop_time\");\n"];
  pr["double evolveshift     = Cal( (grid->time>=shift_st)*(shift_st>=0), 0.0, 1.0);\n"];

  pr["double lapsepsipower   = Getd(\"BSSN_lapsepsipower\");\n"];
  pr["double lapseharmonicf  = Getd(\"BSSN_lapseharmonicf\");\n"];
  pr["double shiftpsipower   = Getd(\"BSSN_shiftpsipower\");\n"];
  pr["double shiftalphapower = Getd(\"BSSN_shiftalphapower\");\n"];
  pr["double shiftgammacoeff = Getd(\"BSSN_shiftgammacoeff\");\n"];
  pr["double shiftdriver     = Getd(\"BSSN_shiftdriver\");\n"];
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
  prdecvl[{ g[a,b],  phi,  A[a,b],  K,  G[a],  alpha,  beta[a], B[a],  alphaDensity}, "ucur"];
  prdecvl[{ng[a,b], nphi, nA[a,b], nK, nG[a], nalpha, nbeta[a],nB[a], nalphaDensity}, "unew"];
  prdecvl[{pg[a,b], pphi, pA[a,b], pK, pG[a], palpha, pbeta[a],pB[a], palphaDensity}, "upre"];

  prdecvlindices[{g[a,b], phi, A[a,b], K, G[a], alpha, beta[a], B[a], alphaDensity}, "ucur"];

  prdecvarname[{dgt[a,b,c]},    "ADMvars_dgxxx"];
  prdecvarname[{ddgt[a,b,c,d]}, "ADMvars_ddgxxxx"];
  prdecvarname[{dAt[a,b,c]},    "ADMvars_dKxxx"];
  prdecvarname[{dphi[a]},       "BSSN_dphix"];
  prdecvarname[{ddphi[a,b]},    "BSSN_ddphixx"];
  prdecvarname[{dGt[a,b]},      "BSSN_dGxx"];
  prdecvarname[{dK[a]},         "BSSN_dKx"];
  prdecvarname[{dalp[a]},       "BSSN_dalpx"];
  prdecvarname[{ddalp[a,b]},    "BSSN_ddalpxx"];
  prdecvarname[{dbeta[a,b]},    "BSSN_dbetaxx"];
  prdecvarname[{ddbeta[a,b,c]}, "BSSN_ddbetaxxx"];
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

(* ADMconstraints.m 
   Wolfgang Tichy 1/2004, Bernd Bruegmann 11/02 *)

(* compute ADM constraints from ADM variables *)


(* variables *)
variables = {g[a,b], K[a,b], alpha, beta[a],
             rho, j[a], S[a,b], psi, dpop[a], ddpop[a,b], 
             ham, mom[a], trK, dtrKdt, normham, normmom[a], 
             dg[a,b,c], ddg[a,b,c,d], dK[a,b,c], dalp[a], ddalp[a,b]}



(* compute in this order *)
tocompute = {

  (* partial derivatives *)
  delg[c,a,b] == OD[g[a,b], c],          (* del[c,g[a,b]], *)
  deldelg[a,b,c,d] == OD2[g[c,d], a,b],  (* deldel[a,b,g[c,d]], *)
  delK[a,b,c] == OD[K[b,c], a],          (* del[a, K[b,c]], *)

  (* make transition to physical metric *)
  f == psi^4,
  delf[a] == 4 f dpop[a],
  deldelf[a,b] == 4 f (ddpop[a,b] + 3 dpop[a] dpop[b]), 
  deldelg[a,b,c,d] == f deldelg[a,b,c,d] + deldelf[a,b] g[c,d] +
                      delf[a] delg[b,c,d] + delf[b] delg[a,c,d],  
  delg[c,a,b] == f delg[c,a,b] + g[a,b] delf[c],

  (* inverse physical metric *)
  detginvf == 1/(f matrixdet[g]),
  ginv[a,b] == detginvf matrixinvdet[g,a,b],

  (* connection of physical metric *) 
  gammado[c,a,b] == 1/2 (delg[a,b,c] + delg[b,a,c] - delg[c,a,b]),
  gamma[c,a,b] == ginv[c,d] gammado[d,a,b], 

  (* curvature of physical metric *)
  R[a,b] == ginv[c,d] ( 1/2 (
    - deldelg[c,d,a,b] - deldelg[a,b,c,d] + 
      deldelg[a,c,b,d] + deldelg[b,c,a,d]) +
      gamma[e,a,c] gammado[e,b,d] - 
      gamma[e,a,b] gammado[e,c,d]),
  R == ginv[a,b] R[a,b],

  (* extrinsic curvature terms *)
  Kud[a,b] == ginv[a,c] K[c,b],
  K == Kud[a,b] delta[a,b],
  KudKud == Kud[a,b] Kud[b,a],
  codelK[a,b,c] == delK[a,b,c] - gamma[d,a,b] K[d,c] - gamma[d,a,c] K[b,d],
  cdKudd[a,b,c] == ginv[a,d] codelK[d,b,c],

  (* use matter terms, but only if the vars are enabled *)
  Cinstruction == "if(rho==NULL) {",
    hamrhs == 0,
  Cinstruction == "} else {",
    hamrhs == 16 Pi rho,
  Cinstruction == "}",

  Cinstruction == "if(j1==NULL) {",
    momrhs[a] == 0,
  Cinstruction == "} else {",
    momrhs[a] == 8 Pi j[a],
  Cinstruction == "}",

  (* Hamiltonian constraint *)
  ham == R + K K - KudKud - hamrhs,

  (* momentum constraint *)
  momu[a] == ginv[a,b] cdKudd[c,b,d] delta[c,d] - ginv[b,c] cdKudd[a,b,c],
  mom[a] == f g[a,b] momu[b] - momrhs[a],

  (* trace of extrinsic curvature *)
  trK == K,

  (* time derivative of trace of extrinsic curvature *)
  (* derivatives of the lapse *)
  cdda[a,b] == ddalp[a,b] - gamma[c,a,b] dalp[c],
  trcdda == ginv[a,b] cdda[a,b],
  dtrK[a] == ginv[b,c] codelK[a,b,c],
  betadtrK == beta[a] dtrK[a],
  Cinstruction == "if(S11==NULL) {",
    Sminus3rho == 0,
  Cinstruction == "} else {",
    Sminus3rho == ginv[a,b] S[a,b] - 3 rho,
  Cinstruction == "}",
  dtrKdt == -trcdda + alpha (R + K*K) + betadtrK + alpha * 4 PI Sminus3rho,

  (* compute normalized constraints *)
  Cif == normConstr,

    (* normalized Hamiltonian constraint *)
    (* normham == ham/(fabs[R] + K K + fabs[Kud[a,b] Kud[b,a]] + fabs[hamrhs]), *)
    Cif == TermByTerm,
      (* compute fabs of some terms in R separately *)
      RA[a,b] == ginv[c,d] * (-deldelg[c,d,a,b]),
      RB[a,b] == ginv[c,d] * ( deldelg[a,c,b,d]),
      RC[a,b] == ginv[c,d] * ( gamma[e,a,c] gammado[e,b,d]),
      RD[a,b] == ginv[c,d] * (-gamma[e,a,b] gammado[e,c,d]),
      RA == ginv[a,b] RA[a,b],
      RB == ginv[a,b] RB[a,b],
      RC == ginv[a,b] RC[a,b],
      RD == ginv[a,b] RD[a,b],
      (* hamnum == RA + RB + RC + RD + K K - KudKud - hamrhs, *)
      denom  == fabs[RA]+fabs[RB]+fabs[RC]+fabs[RD] + K K + fabs[KudKud] +
                fabs[hamrhs],
    Cif == else,
      denom  == fabs[R] + K K + fabs[KudKud] + fabs[hamrhs],
    Cif == end,

    normham == Cal[denom <= 0.0, 0.0, ham/denom],

    (* normalized momentum constraint *)
    (* normmom[a] == mom[a]/( f*( fabs[ cdKudd[c,a,d] delta[c,d] ] 
                                 +fabs[ ginv[b,c] codelK[a,b,c] ] ) ), *)
    Cif == TermByTerm,
      (* compute fabs of some terms in codelK separately *)
      codelKA[a,b,c] ==  delK[a,b,c],
      codelKB[a,b,c] == -gamma[d,a,b] K[d,c],
      codelKC[a,b,c] == -gamma[d,a,c] K[b,d],
      cdKdA[b] == ginv[a,c] codelKA[a,b,c],
      cdKdB[b] == ginv[a,c] codelKB[a,b,c],
      cdKdC[b] == ginv[a,c] codelKC[a,b,c],
      codelTrKA[a] == ginv[b,c] codelKA[a,b,c],
      codelTrKB[a] == ginv[b,c] codelKB[a,b,c],
      codelTrKC[a] == ginv[b,c] codelKC[a,b,c],
      (* note: mom[la] =  CD[ K[la,ue], le ] - CD[ TrK, la ] *)
      (* momnum[a] == cdKdA[a] + cdKdB[a] + cdKdC[a] -
                      codelTrKA[a] - codelTrKB[a] - codelTrKC[a] -momrhs[a], *)
      denom[a] ==  fabs[cdKdA[a]] + fabs[cdKdB[a]] + fabs[cdKdC[a]] +
                   fabs[codelTrKA[a]] + fabs[codelTrKB[a]] +
                   fabs[codelTrKC[a]] + fabs[momrhs[a]],
    Cif == else,
      denom[a]  == fabs[cdKudd[c,a,d] delta[c,d]] + 
                   fabs[ginv[b,c] codelK[a,b,c]] + fabs[momrhs[a]],
    Cif == end,

    denom == denom1 + denom2 + denom3,

    normmom[a] == Cal[denom <= 0.0, 0.0, mom[a]/denom],
 
  Cif == end
}


(* symmetries *)
g[a_,b_]       := g[b,a] /; !OrderedQ[{a,b}]
K[a_,b_]       := K[b,a] /; !OrderedQ[{a,b}]
R[a_,b_]       := R[b,a] /; !OrderedQ[{a,b}]
ginv[a_,b_]    := ginv[b,a] /; !OrderedQ[{a,b}]
ddpop[a_,b_]   := ddpop[b,a] /; !OrderedQ[{a,b}]
deldelf[a_,b_] := deldelf[b,a] /; !OrderedQ[{a,b}]
S[a_,b_]       := S[b,a] /; !OrderedQ[{a,b}]

delg[c_,a_,b_] := delg[c,b,a] /; !OrderedQ[{a,b}]
delK[c_,a_,b_] := delK[c,b,a] /; !OrderedQ[{a,b}]
deldelg[a_,b_,c_,d_] := deldelg[b,a,c,d] /; !OrderedQ[{a,b}]
deldelg[a_,b_,c_,d_] := deldelg[a,b,d,c] /; !OrderedQ[{c,d}]
ddalp[a_,b_]      := ddalp[b,a]    /; !OrderedQ[{a,b}]

codelK[c_,a_,b_] := codelK[c,b,a] /; !OrderedQ[{a,b}]
cdKudd[c_,a_,b_] := cdKudd[c,b,a] /; !OrderedQ[{a,b}]

RA[a_,b_]         := RA[b,a] /; !OrderedQ[{a,b}]
RB[a_,b_]         := RB[b,a] /; !OrderedQ[{a,b}]
RC[a_,b_]         := RC[b,a] /; !OrderedQ[{a,b}]
RD[a_,b_]         := RD[b,a] /; !OrderedQ[{a,b}]
codelKA[c_,a_,b_] := codelKA[c,b,a] /; !OrderedQ[{a,b}]
codelKB[c_,a_,b_] := codelKB[c,b,a] /; !OrderedQ[{a,b}]
codelKC[c_,a_,b_] := codelKC[c,b,a] /; !OrderedQ[{a,b}]
cdKuddA[c_,a_,b_] := cdKuddA[c,b,a] /; !OrderedQ[{a,b}]
cdKuddB[c_,a_,b_] := cdKuddB[c,b,a] /; !OrderedQ[{a,b}]
cdKuddC[c_,a_,b_] := cdKuddC[c,b,a] /; !OrderedQ[{a,b}]

OD2[c_, a_,b_]   := OD2[c,b,a] /; !OrderedQ[{a,b}]

dg[a_,b_,c_]     := dg[b,a,c]  /; !OrderedQ[{a,b}]
dK[a_,b_,c_]     := dK[b,a,c]  /; !OrderedQ[{a,b}]
ddg[a_,b_,c_,d_] := ddg[b,a,c,d] /; !OrderedQ[{a,b}]
ddg[a_,b_,c_,d_] := ddg[a,b,d,c] /; !OrderedQ[{c,d}]

(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be in Cfunctionfile 
*)

CFunctionFile = "ADMconstraints.c"

(* the head of the function *)
BeginCFunction[] := Module[{},

  pr["#include \"sgrid.h\"\n"];
  pr["#include \"ADMvars.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Log(x)     log((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];

  pr["\n\n\n"];

  pr["void ADMconstraints(tVarList *u)\n"];
  pr["{\n"];
  pr["tGrid *grid = u->grid;\n"];
  pr["int bi;\n"];
  pr["\n"];

  pr["int normConstr = Getv(\"ADMvars_normalizedConstraints\", \"yes\");\n"];
  pr["int TermByTerm = GetvLax(\"ADMvars_ConstraintNorm\", \"TermByTerm\");\n"];
  pr["int useDD = GetvLax(\"ADMvars_useDD\", \"yes\");\n"];
  pr["\n"];

  pr["for(bi = 0; bi < grid->nboxes; bi++)\n"];
  pr["{\n"];
  pr["tBox *box = grid->box[bi];\n"];
  pr["int ijk;\n\n"];

  pr["if(useDD) {\n"];
  pr["  allDerivsOf_Sab(box, Ind(\"gxx\"), Ind(\"ADMvars_dgxxx\"),
                        Ind(\"ADMvars_ddgxxxx\"));\n"];
  pr["  allDerivsOf_S(box, Ind(\"alpha\"), Ind(\"ADMvars_dalpx\"),
                      Ind(\"ADMvars_ddalpxx\"));\n"];
  pr["}\n"];
  pr["else {\n"];
  pr["  FirstAndSecondDerivsOf_Sab(box, Ind(\"gxx\"),Ind(\"ADMvars_dgxxx\"),
                                   Ind(\"ADMvars_ddgxxxx\"));\n"];
  pr["  FirstAndSecondDerivsOf_S(box, Ind(\"alpha\"),Ind(\"ADMvars_dalpx\"),
                                 Ind(\"ADMvars_ddalpxx\"));\n"];
  pr["}\n"];
  pr["FirstDerivsOf_Sab(box, Ind(\"Kxx\"), Ind(\"ADMvars_dKxxx\"));\n"];
  pr["\n"];

  pr["forallpoints(box, ijk)\n"];
  pr["{\n"];
];

(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  prdecvl[{psi, dpop[a], ddpop[a,b]}, "psiandderivs"];
  prdecvl[{g[a,b], K[a,b], alpha, beta[a], rho, j[a], S[a,b], ham, mom[a], trK, dtrKdt, normham, normmom[a]}, "u"];
  prdecvarname[{dg[a,b,c]},    "ADMvars_dgxxx"];
  prdecvarname[{ddg[a,b,c,d]}, "ADMvars_ddgxxxx"];
  prdecvarname[{dK[a,b,c]},    "ADMvars_dKxxx"];
  prdecvarname[{dalp[a]},      "ADMvars_dalpx"];
  prdecvarname[{ddalp[a,b]},   "ADMvars_ddalpxx"];

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

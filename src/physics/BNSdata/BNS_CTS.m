(* BNS_CTS.m 
   Wolfgang Tichy  12/2007       *)

(* compute right hand side of BSSN equations *)


(* variables *)
variables = {Psi, beta[a], alphaP, Sigma, FPsi, Fbeta[a], FalphaP ,FSigma,
              dPsi[a],   dbeta[a,b],   dalphaP[a],    dSigma[a],
             ddPsi[a,b],ddbeta[a,b,c],ddalphaP[a,b], ddSigma[a,b],
	     LPsi,Lbeta[a],LalphaP,LSigma, FLPsi,FLbeta[a],FLalphaP,FLSigma,
              dLPsi[a],   dLbeta[a,b],   dLalphaP[a],    dLSigma[a],
             ddLPsi[a,b],ddLbeta[a,b,c],ddLalphaP[a,b], ddLSigma[a,b],
	     g[a,b],  K[a,b], q, vRS[a]}

(* compute in this order *)
tocompute = {

  Cif == nonlin,
    Cinstruction == "FirstAndSecondDerivsOf_S(box, index_Psi, \ 
			Ind(\"BNSdata_Psix\"), Ind(\"BNSdata_Psixx\"));",
    Cinstruction == "FirstAndSecondDerivsOf_Sa(box, index_beta1, \
			Ind(\"BNSdata_betaxx\"), Ind(\"BNSdata_betaxxx\"));",
    Cinstruction == "FirstAndSecondDerivsOf_S(box, index_alphaP, \
			Ind(\"BNSdata_alphaPx\"), Ind(\"BNSdata_alphaPxx\"));",
    Cinstruction == "FirstAndSecondDerivsOf_S(box, index_Sigma, \
			Ind(\"BNSdata_Sigmax\"), Ind(\"BNSdata_Sigmaxx\"));",
  Cif == else,
    Cinstruction == "FirstAndSecondDerivsOf_S(box, index_LPsi, \ 
					index_dLPsi1, index_ddLPsi11);",
    Cinstruction == "FirstAndSecondDerivsOf_Sa(box, index_Lbeta1, \
					index_dLbeta11, index_ddLbeta111);",
    Cinstruction == "FirstAndSecondDerivsOf_S(box, index_LalphaP, \
					index_dLalphaP1, index_ddLalphaP11);",
    Cinstruction == "FirstAndSecondDerivsOf_S(box, index_LSigma, \
					index_dLSigma1, index_ddLSigma11);",
  Cif == end,

  (* loop of all points *)
  Cinstruction == "forallpoints(box, ijk) {",

  Cinstruction == "
        double xmax1 = grid->box[0]->x_of_X[1](
                        (void *) grid->box[0], 0, 0.0,0.0,0.0);
        double xmin1 = grid->box[0]->x_of_X[1](
                        (void *) grid->box[0], 0, 0.0,1.0,0.0);
        double xmax2 = grid->box[3]->x_of_X[1](
                        (void *) grid->box[3], 0, 0.0,1.0,0.0);
        double xmin2 = grid->box[3]->x_of_X[1](
                        (void *) grid->box[3], 0, 0.0,0.0,0.0);
        double R1  = 0.5*(xmax1-xmin1);
        double R2  = 0.5*(xmax2-xmin2);
	double rh1 = 0.0;
        double rh2 = 0.0;

        if(bi==0 || bi==5)  rh1 = -3.0/(R1*R1*R1);
        if(bi==3 || bi==4)  rh2 = -6.0/(R2*R2*R2); ",

  (* inverse conformal metric *)
  (*  detginv == 1/matrixdet[g],
      ginv[a,b] == detginv matrixinvdet[g,a,b], *)


  Cif == nonlin,
    (* equations for Psi, beta[a], alphaP, Sigma *)
    FPsi     == delta[b,c] ddPsi[b,c] - rh1, 
    Fbeta[a] == delta[b,c] ddbeta[a,b,c] - 0,
    FalphaP  == delta[b,c] ddalphaP[b,c] - 0, 
    FSigma   == delta[b,c] ddSigma[b,c] - rh2, 
  Cif == else,
    (* linearized equations for Psi, beta[a], alphaP, Sigma *)
    FLPsi     == delta[b,c] ddLPsi[b,c] - 0, 
    FLbeta[a] == delta[b,c] ddLbeta[a,b,c] - 0,
    FLalphaP  == delta[b,c] ddLalphaP[b,c] - 0, 
    FLSigma   == delta[b,c] ddLSigma[b,c] - 0, 
  Cif == end,

  Cinstruction == "} /* end of points loop */\n"
}


(* symmetries *)
g[a_,b_] :=  g[b,a] /; !OrderedQ[{a,b}]
K[a_,b_] :=  K[b,a] /; !OrderedQ[{a,b}]

ginv[a_,b_] := ginv[b,a] /; !OrderedQ[{a,b}]
Kup[a_,b_]  := Kup[b,a]  /; !OrderedQ[{a,b}]

ddPsi[a_,b_]     := ddPsi[b,a]    /; !OrderedQ[{a,b}]
ddbeta[a_,b_,c_] := ddbeta[a,c,b] /; !OrderedQ[{b,c}]
ddalphaP[a_,b_]  := ddalphaP[b,a] /; !OrderedQ[{a,b}]
ddSigma[a_,b_]   := ddSigma[b,a]  /; !OrderedQ[{a,b}]

ddLPsi[a_,b_]     := ddLPsi[b,a]    /; !OrderedQ[{a,b}]
ddLbeta[a_,b_,c_] := ddLbeta[a,c,b] /; !OrderedQ[{b,c}]
ddLalphaP[a_,b_]  := ddLalphaP[b,a] /; !OrderedQ[{a,b}]
ddLSigma[a_,b_]   := ddLSigma[b,a]  /; !OrderedQ[{a,b}]

(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be in Cfunctionfile 
*)

CFunctionFile = "BNS_CTS.c"

(* the head of the function *)
BeginCFunction[] := Module[{},

  pr["#include \"sgrid.h\"\n"];
  pr["#include \"BNSdata.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Log(x)     log((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];

  pr["\n\n\n"];

  pr["void BNS_CTS(tVarList *vlFu, tVarList *vlu, \ 
		   tVarList *vlJdu, tVarList *vldu, tVarList *vlduDerivs, \
		   int nonlin)\n"];
  pr["{\n"];
  pr["tGrid *grid = vlu->grid;\n"];
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

  prdecvl[{FPsi, Fbeta[a], FalphaP ,FSigma}, "vlFu"];
  prdecvl[{ Psi,  beta[a],  alphaP,  Sigma}, "vlu"];

  prdecvl[{FLPsi, FLbeta[a], FLalphaP ,FLSigma}, "vlJdu"];
  prdecvl[{ LPsi, Lbeta[a], LalphaP, LSigma}, "vldu"];
  prdecvl[{dLPsi[a],ddLPsi[a,b], dLbeta[a,b],ddLbeta[a,b,c], dLalphaP[a],ddLalphaP[a,b], dLSigma[a],ddLSigma[a,b]}, "vlduDerivs"];

  prdecvlindices[{ Psi,  beta[a],  alphaP,  Sigma}, "vlu"];
  prdecvlindices[{LPsi, Lbeta[a], LalphaP, LSigma}, "vldu"];
  prdecvlindices[{dLPsi[a],ddLPsi[a,b], dLbeta[a,b],ddLbeta[a,b,c], dLalphaP[a],ddLalphaP[a,b], dLSigma[a],ddLSigma[a,b]}, "vlduDerivs"];

  prdecvarname[{dPsi[a]},       "BNSdata_Psix"];
  prdecvarname[{ddPsi[a,b]},    "BNSdata_Psixx"];
  prdecvarname[{dbeta[a,b]}, 	"BNSdata_betaxx"];
  prdecvarname[{ddbeta[a,b,c]}, "BNSdata_betaxxx"];
  prdecvarname[{dalphaP[a]},    "BNSdata_alphaPx"];
  prdecvarname[{ddalphaP[a,b]}, "BNSdata_alphaPxx"];
  prdecvarname[{dSigma[a]},     "BNSdata_Sigmax"];
  prdecvarname[{ddSigma[a,b]},  "BNSdata_Sigmaxx"];

  prdecvarname[{g[a,b]}, "gxx"];
  prdecvarname[{K[a,b]}, "Kxx"];
  prdecvarname[{q},      "BNSdata_q"];
  prdecvarname[{vRS[a]}, "BNSdata_vRSx"];

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


(************************************************************************)
(* now we are ready to go *)

<< "../../Math/MathToC/TensorEquations3dToC.m"

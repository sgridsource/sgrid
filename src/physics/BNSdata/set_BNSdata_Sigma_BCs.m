(* BNS_CTS.m 
   Wolfgang Tichy  12/2007       *)

(* compute residuals of BNS ham, mom, alphaP and Sigma eqns *)


(* variables *)
variables = {Psi, B[a], alphaP, Sigma, FPsi, FB[a], FalphaP ,FSigma,
             dSigma[a],
	     lPsi,lB[a],lalphaP,lSigma, FlPsi,FlB[a],FlalphaP,FlSigma,
              dlPsi[a],   dlB[a,b],   dlalphaP[a],    dlSigma[a],
             ddlPsi[a,b],ddlB[a,b,c],ddlalphaP[a,b], ddlSigma[a,b],
             q, wB[a], dq[a], x, y}

constvariables = {OmegaCrossR[a]}

(* compute in this order *)
tocompute = {

  (* do nothing in boxes 4 and 5 inside stars *)
  Cif == (bi>=4),
    Cinstruction == "continue;",
  Cif == end,

  (* Use Sigma=0 as BC if corot *)
  Cif == ( ((bi==0 || bi==1) && corot1) || ((bi==2 ||bi==3) && corot2) ),

    Cif == nonlin, (* non-linear case *)
      (* go over A=0 plane *)
      Cinstruction == "forplane1(i,j,k, n1,n2,n3, 0){ ijk=Index(i,j,k);",
        FSigma  == Sigma,  (* set Sigma=0 *)
      Cinstruction == "} /* end forplane1 */",
      (* go over A=1 plane *)
      Cinstruction == "forplane1(i,j,k, n1,n2,n3, n1-1){ ijk=Index(i,j,k);",
        FSigma  == Sigma,  (* set Sigma=0 *)
      Cinstruction == "} /* end forplane1 */",

    Cif == else,   (* linear case *)
      (* go over A=0 plane *)
      Cinstruction == "forplane1(i,j,k, n1,n2,n3, 0){ ijk=Index(i,j,k);",
        FlSigma  == lSigma,
      Cinstruction == "} /* end forplane1 */",
      (* go over A=1 plane *)
      Cinstruction == "forplane1(i,j,k, n1,n2,n3, n1-1){ ijk=Index(i,j,k);",
        FlSigma  == lSigma,  (* set Sigma=0 *)
      Cinstruction == "} /* end forplane1 */",

    Cif == end,

    Cinstruction == "continue; /* for corot we are done with this box */",
  Cif == end,

  (* Use Sigma=0 as BC on A=1 if (bi==1 || bi==2) for corot and general cases *)
  Cif == ( bi==1 || bi==2 ),
    Cif == nonlin, (* non-linear case *)
      Cinstruction == "forplane1(i,j,k, n1,n2,n3, n1-1){ ijk=Index(i,j,k);",
        FSigma  == Sigma,  (* set Sigma=0 *)
      Cinstruction == "} /* end forplane1 */",
    Cif == else,   (* linear case *)
      Cinstruction == "forplane1(i,j,k, n1,n2,n3, n1-1){ ijk=Index(i,j,k);",
        FlSigma  == lSigma,  (* set Sigma=0 *)
      Cinstruction == "} /* end forplane1 */",
    Cif == end,
  Cif == end,

  (* Put the condition that 
     Sigma(A=0) in box 1/2 equals Sigma(A=0) in box 0/3
     in FSigma and FlSigma at position with i=1.
     Note: FSigma and FlSigma at position with i=0 are already set by 
           set_BNSdata_BCs such that the normal derivs in adjacent boxes
           are equal. *)
  Cif == ( bi==1 || bi==2 ),
    Cinstruction == "int    ind, indin, biin, n1in,n2in,n3in;",
    Cinstruction == "double Sig, Sigin;",
    Cinstruction == "double *Sigmain;",
    Cinstruction == "double *lSigmain;",
    Cinstruction == "if(bi==1) biin=0; else biin=3;",
    Cinstruction == "n1in = grid->box[biin]->n1;
                     n2in = grid->box[biin]->n2;
                     n3in = grid->box[biin]->n3;\n
                     Sigmain  = grid->box[biin]->v[index_Sigma];
                     lSigmain = grid->box[biin]->v[index_lSigma];",
    Cinstruction == "if(n2in!=n2 || n3in!=n3) errorexit(\"we need n2in=n2 and n3in=n3\");",
    Cif == nonlin, (* non-linear case *)
      Cinstruction == "forplane1(i,j,k, n1,n2,n3, 1){ ijk=Index(i,j,k);",
        Cinstruction == "ind   = Ind_n1n2(0,j,k,n1,n2);
                         indin = Ind_n1n2(0,j,k,n1in,n2in);
                         Sig   = Sigma[ind];
                         Sigin = Sigmain[indin];",
        Cif == (SigmaZeroInOutBoxAtA0B0 && j==0),
          FSigma  == Sig,        (* set Sigma=0 at A=B=0 in outer box *)
        Cif == else,
          FSigma  == Sig-Sigin,  (* set Sigmas equal at A=0*)
        Cif == end,
      Cinstruction == "} /* end forplane1 */",
    Cif == else,   (* linear case *)
      Cinstruction == "forplane1(i,j,k, n1,n2,n3, 1){ ijk=Index(i,j,k);",
        Cinstruction == "ind   = Ind_n1n2(0,j,k,n1,n2);
                         indin = Ind_n1n2(0,j,k,n1in,n2in);
                         Sig   = lSigma[ind];
                         Sigin = lSigmain[indin];",
        Cif == (SigmaZeroInOutBoxAtA0B0 && j==0),
          FlSigma  == Sig,        (* set Sigma=0 at A=B=0 in outer box *)
        Cif == else,
          FlSigma  == Sig-Sigin,  (* set Sigmas equal at A=0*)
        Cif == end,
      Cinstruction == "} /* end forplane1 */",
    Cif == end,

  (* if not b=1 or 2, i.e. bi==0 || bi==3 *)
  Cif == else,
    (* if we get here bi=0 or 3 and there is no corot in this box *)
    Cinstruction == "FirstDerivsOf_S(box,  Ind(\"BNSdata_q\"), \
			                 Ind(\"BNSdata_qx\"));",
    (* non-linear case: *)
    Cif == nonlin,
      Cinstruction == "FirstDerivsOf_S(box, index_Sigma, \
                                       Ind(\"BNSdata_Sigmax\"));",
      (* go over A=0 plane *)
      Cinstruction == "forplane1(i,j,k, n1,n2,n3, 0){ ijk=Index(i,j,k);",

        (* Omega \times r term *)
        OmegaCrossR1 == - Omega y,
        OmegaCrossR2 == + Omega (x-xCM),
        OmegaCrossR3 == 0,

        (* shift in rotating frame (in the inertial frame beta^i = B^i) *)
        beta[a] == B[a] + OmegaCrossR[a],

        (* some abbreviations *)
        Psi2  == Psi*Psi,
        Psi4  == Psi2*Psi2,
        Psim2 == 1/Psi2,
        Psim4 == Psim2*Psim2,
        Psim6 == Psim4*Psim2,
        alpha == alphaP/Psi,
        alpha2 == alpha alpha,
        (* terms needed *)
        h == (n+1) q + 1,
        h2 == h h,
        DSigmaUp[a] == Psim4 dSigma[a],
        dSigmaUp[a] == dSigma[a],
        w[a] == Psim6 wB[a],
        wBDown[a] == wB[a],
        wDown[a] == Psim2 wBDown[a],
        L2 == h2 + (wDown[c] + dSigma[c]) (w[c] + DSigmaUp[c]),
        uzerosqr == L2/(alpha2 h2),
        uzero == sqrt[uzerosqr],

        FSigma == dSigmaUp[c] dq[c] - h uzero Psi4 beta[c] dq[c],
        (* add extra term with wB *)
        FSigma == FSigma + Psim2 wB[c] dq[c],
(*
FSigma == Sigma - 10,
FSigma == dSigmaUp[c] dq[c] - 100 dq2,
FSigma == dSigmaUp[c] dq[c] - h uzero Psi4 beta[c] dq[c],
*)

      Cinstruction == "} /* end forplane1 */",

      (* set Sigma and to zero at A=0, B=0 (one point at xout1/2) *)
      Cif == SigmaZeroAtA0B0,
        Cinstruction == "i=0;  j=0;",
        Cinstruction == "for(k=0; k<n3; k++){ ijk=Index(i,j,k);",
          FSigma == Sigma,
        Cinstruction == "} /* end for k  */",
      Cif == end,

    (* linear case: *)
    Cif == else,
      Cinstruction == "FirstDerivsOf_S(box, index_lSigma, index_dlSigma1);",

      (* go over A=0 plane *)
      Cinstruction == "forplane1(i,j,k, n1,n2,n3, 0){ ijk=Index(i,j,k);",

        (* Omega \times r term *)
        OmegaCrossR1 == - Omega y,
        OmegaCrossR2 == + Omega (x-xCM),
        OmegaCrossR3 == 0,

        (* shift in rotating frame (in the inertial frame beta^i = B^i) *)
        beta[a] == B[a] + OmegaCrossR[a],

        (* some abbreviations *)
        Psi2  == Psi*Psi,
        Psi3  == Psi2*Psi,
        Psi4  == Psi2*Psi2,
        Psim2 == 1/Psi2,
        Psim3 == 1/Psi3,
        Psim4 == Psim2*Psim2,
        Psim6 == Psim4*Psim2,
        Psim8 == Psim4*Psim4,
        Psim5 == Psim6*Psi,
        Psim7 == Psim5*Psim2,
        Psim9 == Psim7*Psim2,
        (* non-linear terms needed *)
        h == (n+1) q + 1,
        h2 == h h,
        DSigmaUp[a] == Psim4 dSigma[a],
        dSigmaUp[a] == dSigma[a],
        w[a] == Psim6 wB[a],
        wBDown[a] == wB[a],
        wDown[a] == Psim2 wBDown[a],
        L2 == h2 + (wDown[c] + dSigma[c]) (w[c] + DSigmaUp[c]),
        uzerosqr == L2/(alpha2 h2),
        uzero == sqrt[uzerosqr],
        alpha == alphaP/Psi,
        alpha2 == alpha alpha,

        (* linearized terms *)
        lq     == 0,
        dlq[a] == 0,
        lh   == 0,
        lLnh == 0,
        (* wB remains const under linearization *)
        lwB[a] == 0,
        dlwB[a,b] == 0,
        lalpha == lalphaP/Psi - alphaP lPsi/Psi2,   
        (* dSigmaUp[a] == dSigma[a], *)
        dlSigmaUp[a] == dlSigma[a],
        lL2 == 2*(Psim8 wBDown[c] lwB[c] +
                     Psim6 (lwB[c] dSigma[c] + wB[c] dlSigma[c]) +
                     Psim4 dSigmaUp[c] dlSigma[c]  - 
                  (8 Psim9 wBDown[c] wB[c] + 12 Psim7 wB[c] dSigma[c] +
                   4 Psim5 dSigma[c] dSigmaUp[c]) lPsi + 2 h2 lLnh),
        luzerosqr == (lL2 - 2 L2 (lalpha/alpha + lLnh))/(alpha2 h2),
        luzero == luzerosqr/(2 uzero),
        lhuzeroPsi4beta[a] == h (luzero Psi4 beta[a] +
                                 4 uzero Psi3 lPsi beta[a] +
                                 uzero Psi4 lB[a]) + lh uzero Psi4 beta[a],
        
        FlSigma  == dlSigmaUp[c] dq[c] - lhuzeroPsi4beta[c] dq[c] +
                    dSigmaUp[c] dlq[c] - h uzero Psi4 beta[c] dlq[c],
        (* add extra term with wB *)
        FlSigma == FlSigma + Psim2 lwB[c] dq[c] - 2 Psim3 lPsi wB[c] dq[c] +
                             Psim2 wB[c] dlq[c],
(*
lL2 == 2*(dlSigma[c] (w[c] + DSigmaUp[c])),
luzerosqr == (lL2)/(alpha2 h2),
luzero == luzerosqr/(2 uzero), 
FlSigma  == lSigma,
FlSigma == dlSigmaUp[c] dq[c],
FlSigma == dlSigmaUp[c] dq[c] - h luzero Psi4 beta[c] dq[c],
*)

      Cinstruction == "} /* end forplane1 */",

      (* set Sigma and to zero at A=0, B=0 (one point at xout1/2) *)
      Cif == SigmaZeroAtA0B0,
        Cinstruction == "i=0;  j=0;",
        Cinstruction == "for(k=0; k<n3; k++){ ijk=Index(i,j,k);",
          FlSigma == lSigma,
        Cinstruction == "} /* end for k  */",
      Cif == end,

    Cif == end, (* end of nonlin/linear case *)

  Cif == end,

  Cinstruction == "/* end all */\n"
}


(* symmetries *)
ddPsi[a_,b_]     := ddPsi[b,a]    /; !OrderedQ[{a,b}]
ddB[a_,b_,c_]    := ddB[a,c,b]    /; !OrderedQ[{b,c}]
ddalphaP[a_,b_]  := ddalphaP[b,a] /; !OrderedQ[{a,b}]
ddSigma[a_,b_]   := ddSigma[b,a]  /; !OrderedQ[{a,b}]

ddlPsi[a_,b_]     := ddlPsi[b,a]    /; !OrderedQ[{a,b}]
ddlB[a_,b_,c_]    := ddlB[a,c,b]    /; !OrderedQ[{b,c}]
ddlalphaP[a_,b_]  := ddlalphaP[b,a] /; !OrderedQ[{a,b}]
ddlSigma[a_,b_]   := ddlSigma[b,a]  /; !OrderedQ[{a,b}]

(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be in Cfunctionfile 
*)

CFunctionFile = "set_BNSdata_Sigma_BCs.c"

(* the head of the function *)
BeginCFunction[] := Module[{},

  pr["#include \"sgrid.h\"\n"];
  pr["#include \"BNSdata.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Log(x)     (log((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];

  pr["\n\n\n"];

  pr["void set_BNSdata_Sigma_BC(tVarList *vlFu, tVarList *vlu, \ 
		   tVarList *vlJdu, tVarList *vldu, tVarList *vlduDerivs, \
		   int nonlin)\n"];
  pr["{\n"];

  pr["int corot1 = Getv(\"BNSdata_rotationstate1\",\"corotation\");\n"];
  pr["int corot2 = Getv(\"BNSdata_rotationstate2\",\"corotation\");\n"];
  pr["int SigmaZeroAtA0B0 = Getv(\"BNSdata_Sigma_surface_BCs\",\"ZeroAt00\");\n"];
  pr["int SigmaZeroInOutBoxAtA0B0 = Getv(\"BNSdata_Sigma_surface_BCs\",\"ZeroInOuterBoxAt00\");\n"];
  pr["int noBCs = Getv(\"BNSdata_Sigma_surface_BCs\",\"none\");\n"];
  pr["double n = Getd(\"BNSdata_n\");\n"];
  pr["double kappa = Getd(\"BNSdata_kappa\");\n"];
  pr["double Omega = Getd(\"BNSdata_Omega\");\n"];
  pr["double xCM = Getd(\"BNSdata_x_CM\");\n"];
  pr["\n"];

  pr["tGrid *grid = vlu->grid;\n"];
  pr["int bi;\n"];
  pr["\n"];
  pr["\n"];

  pr["/* do nothing if noBCs, i.e. BNSdata_Sigma_surface_BCs = none */\n"];
  pr["if(noBCs) return;\n\n\n"];

  pr["forallboxes(grid,bi)\n"];
  pr["{\n"];
  pr["tBox *box = grid->box[bi];\n"];
  pr["int ijk;\n\n"];
  pr["int n1 = box->n1;\n"];
  pr["int n2 = box->n2;\n"];
  pr["int n3 = box->n3;\n"];
  pr["int i,j,k;\n\n"];

  pr["\n"];
  pr["\n"];
];

(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  prdecvl[{FPsi, FB[a], FalphaP ,FSigma}, "vlFu"];
  prdecvl[{ Psi,  B[a],  alphaP,  Sigma}, "vlu"];

  prdecvl[{FlPsi, FlB[a], FlalphaP ,FlSigma}, "vlJdu"];
  prdecvl[{ lPsi, lB[a], lalphaP, lSigma}, "vldu"];
  prdecvl[{dlPsi[a],ddlPsi[a,b], dlB[a,b],ddlB[a,b,c], dlalphaP[a],ddlalphaP[a,b], dlSigma[a],ddlSigma[a,b]}, "vlduDerivs"];

  prdecvlindices[{ Psi,  B[a],  alphaP,  Sigma}, "vlu"];
  prdecvlindices[{lPsi, lB[a], lalphaP, lSigma}, "vldu"];
  prdecvlindices[{dlPsi[a],ddlPsi[a,b], dlB[a,b],ddlB[a,b,c], dlalphaP[a],ddlalphaP[a,b], dlSigma[a],ddlSigma[a,b]}, "vlduDerivs"];

  prdecvarname[{dSigma[a]},     "BNSdata_Sigmax"];

  prdecvarname[{x},      "x"];
  prdecvarname[{y},      "y"];

  (* prdecvarname[{g[a,b]}, "gxx"]; prdecvarname[{K[a,b]}, "Kxx"]; *)
  prdecvarname[{q},       "BNSdata_q"];
  prdecvarname[{wB[a]},   "BNSdata_wBx"];
  prdecvarname[{dq[a]},   "BNSdata_qx"];

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

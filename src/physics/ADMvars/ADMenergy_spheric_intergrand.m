(* ADMconstraints.m 
   Bernd Bruegmann 11/02, Wolfgang Tichy 1/2004 *)

(* compute ADM constraints from ADM variables *)


(* variables *)
variables = {Eadm, gb[a,b], x[a], psi, dpop[a], dgb[a,b,c]}


(* compute in this order *)
tocompute = {

  (* partial derivatives *)
  delg[c,a,b] == OD[gb[a,b], c],          (* del[c,gb[a,b]], *)

  (* make transition to physical metric *)
  f == psi^4,
  delf[a] == 4 f dpop[a],
  delg[c,a,b] == f delg[c,a,b] + gb[a,b] delf[c],

  (* det of gb and physical g *)
  detgb == matrixdet[gb],
  detg  == f^3 detgb,
 
  (* inverse physical metric *)
  ginv[a,b] == matrixinvdet[gb,a,b]/(detgb*f),

  (* radius^2 *)
  r2 == x1*x1 + x2*x2 + x3*x3,

  Cif == r2,
    (* normal *)
    n[a] == x[a]/sqrt[r2],
    (* Integrand *)
    integrand == (ginv[a,b] ginv[c,d]) * (delg[d,b,c] - delg[b,c,d]) *
                 (detg) * n[a] / (16.0 * PI),
  Cif == else,
    integrand == 0.0,
  Cif == end,

  Eadm == integrand
}


(* symmetries *)
gb[a_,b_]      := gb[b,a] /; !OrderedQ[{a,b}]
ginv[a_,b_]    := ginv[b,a] /; !OrderedQ[{a,b}]

delg[c_,a_,b_] := delg[c,b,a] /; !OrderedQ[{a,b}]
dgb[a_,b_,c_]  := dgb[b,a,c]  /; !OrderedQ[{a,b}]

(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be in Cfunctionfile 
*)

CFunctionFile = "ADMenergy_spheric_intergrand.c"

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

  pr["void ADMenergy_spheric_intergrand(tVarList *u)\n"];
  pr["{\n"];
  pr["tGrid *grid = u->grid;\n"];
  pr["int bi;\n"];
  pr["\n"];

  pr["for(bi = 0; bi < grid->nboxes; bi++)\n"];
  pr["{\n"];
  pr["tBox *box = grid->box[bi];\n"];
  pr["int ijk;\n\n"];

  pr["FirstDerivsOf_Sab(box, Ind(\"gxx\"), Ind(\"ADMvars_dgxxx\"));\n"];
  pr["\n"];

  pr["forallpoints(box, ijk)\n"];
  pr["{\n"];
];

(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  prdecvl[{Eadm, gb[a,b], x[a], psi, dpop[a]}, "u"];
  prdecvarname[{dgb[a,b,c]}, "ADMvars_dgxxx"];

];    
(* auxillary variables are automatically inserted here *)

(* custom C-code which is placed directly after the variable declarations *)
InitializationCommands[] := Module[{},

  pr["/* Jetzt geht's los! */\n"];

];

(* the end or tail of the function (we need at least a }) *)
EndCFunction[] := Module[{},

  pr["} /* end of points */\n\n"];

(*
  pr["snprintf(str, 999, \"box%d_Coordinates\", box->b);\n"];
  pr["if( Getv(str, \"SphericalDF\") )\n"];
  pr["{ double *Eadm = vlldataptr(u, box, 0);\n"];
  pr["  spec_sphericalDF2dIntegral(box, Eadm, Eadm);\n"];
  pr["}\n"];
  pr["else\n"];
  pr["errorexits(\"ADMenergy_spheric_intergrand: I don't know how to do surface integrals \"\n"];
  pr["\"in %s coordinates\", str);\n\n"];
*)

  pr["} /* end of boxes */\n"];
  pr["\n\n"];
  pr["}  /* end of function */\n\n"];
];


(* to turn off optimization for development set optimizeflag = False *)
optimizeflag = True;


(************************************************************************)
(* now we are ready to go *)

<< "../../Math/MathToC/TensorEquationsToC.m"

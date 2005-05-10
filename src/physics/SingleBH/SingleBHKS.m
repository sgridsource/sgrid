(*  SingleBHKS.m 
   Wolfgang Tichy 12/2002 *)

(* Compute gb[a,b], K[a,b], psi from x,y,z and M=BHmass1 *)


(* variables *)
variables = { x[a], gb[a,b], K[a,b], psi, dpsiopsi[a], ddpsiopsi[a,b], 
	      alpha, beta[a] }

(* compute in this order *)
tocompute = {
 
 rho == sqrt[ delta[a,b] x[a] x[b] ] ,
 aS[a] == If[a==1,Sx,If[a==2,Sy,Sz]],
 adotx == delta[a,b] aS[a] x[b],
 a2 == delta[a,b] aS[a] aS[b],
 Omega == sqrt[(rho^2 - a2)^2 + 4 adotx^2],
 r == sqrt[(rho^2 - a2 + Omega)/2],
 dr[a] == 1/(2 r) ( 1 + (rho^2 - a2)/Omega) x[a] + adotx/(Omega r) aS[a],
 acrossx[a] == epsilon[a,b,c] aS[b] x[c],
 ahat[a,b] == epsilon[a,c,b] aS[c],

 H == M r^3/(r^4 + adotx^2),
 dH[a] == H/r (( 3 adotx^2 - r^4) dr[a] - 2 r adotx aS[a])/(r^4 + adotx^2),

 l[a] == (r^2 x[a] - r acrossx[a] + adotx aS[a])/(r (r^2 + a2)),

 dl[a,b] == (r^2 delta[a,b] - r ahat[a,b] + aS[a] aS[b])/(r^3 + a2 r) + (((a2 r^2 - r^4) x[a] + 2 r^3 acrossx[a] - (3 r^2 + a2) adotx aS[a]) dr[b])/(r^3 + a2 r)^2,

 (* now we can construct the initial data *)

 (* the metric *)
 g[a,b] == delta[a,b] + 2 H l[a] l[b],

 (* inverse metric *)
 detgup == 1/matrixdet[g], 
 gup[a,b] == detgup matrixinvdet[g,a,b],

 (* Lapse and shift *)
 alpha == 1/sqrt[1 + 2 H],
 beta[a] == 2 H/(1 + 2 H) l[a],

 K[a,b] == alpha ( l[a] dH[b] + l[b] dH[a] + H dl[a,b] + H dl[b,a] + 2 H^2 (l[a] l[c] dl[b,c] + l[b] l[c] dl[a,c]) + 2 H l[a] l[b] l[c] dH[c]),

 (* Conformal Factor *)
 Cif == ConformalFactor,

   (* set cconformal factor psi=1+2M/r   <-adhoc choice *)
   psi == 1+2M/r,
   dpsiopsi[a] == ( -(2M/r^2) x[a]/r )/psi ,  (* = ( -(2M/r^3) x[a] )/psi *)
   ddpsiopsi[a,b] == ( (6M/r^3) (x[a]/r) (x[b]/r) -(2M/r^3)delta[a,b] )/psi,

 Cif == else,

   (* set conformal factor to 1 *)
   psi ==1,
   dpsiopsi[a] == 0,
   ddpsiopsi[a,b] == 0,

 Cif == end,

 (* background metric *)
 gb[a,b] == g[a,b]/psi^4
 
}


(* symmetries *)

epsilon[1,2,3] := 1
epsilon[a_,b_,c_] := -epsilon[a,c,b] /; !OrderedQ[{b,c}]
epsilon[a_,b_,b_] := 0
epsilon[a_,b_,c_] := -epsilon[c,b,a] /; !OrderedQ[{a,c}]
epsilon[a_,b_,a_] := 0
epsilon[a_,b_,c_] := -epsilon[b,a,c] /; !OrderedQ[{a,b}]
epsilon[a_,a_,c_] := 0

g[a_,b_]   := g[b,a]   /; !OrderedQ[{a,b}]
gb[a_,b_]  := gb[b,a]  /; !OrderedQ[{a,b}]
gup[a_,b_] := gup[b,a] /; !OrderedQ[{a,b}]
ahat[a_,b_] := -ahat[b,a] /; !OrderedQ[{a,b}]
ahat[a_,a_] := 0

K[a_,b_] := K[b,a] /; !OrderedQ[{a,b}]

ddpsiopsi[a_,b_] := ddpsiopsi[b,a] /; !OrderedQ[{a,b}]


(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be in Cfunctionfile 
*)

CFunctionFile = "SingleBHKS.c"

(* the head of the function *)
BeginCFunction[] := Module[{},

  pr["#include \"sgrid.h\"\n"];
  pr["#include \"SingleBH.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Log(x)     log((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];

  pr["\n\n\n"];

  pr["void SingleBHKS(tGrid *grid, int i_x, int i_gb, int i_K, 
                           int i_psi, int i_dpsiopsi, int i_ddpsiopsi,
                           int i_alpha, int i_beta)\n"];
  pr["{\n"];
  pr["int bi;\n"];
  pr["\n"];

  pr["double M = Getd(\"BHmass1\");\n"];
  pr["double Sx = Getd(\"BHsx1\");\n"];
  pr["double Sy = Getd(\"BHsy1\");\n"];
  pr["double Sz = Getd(\"BHsz1\");\n"];   
  pr["int ConformalFactor = Getv(\"SingleBH_ConformalFactor\",\"yes\");\n"];
  pr["\n"];

  pr["for(bi = 0; bi < grid->nboxes; bi++)\n"];
  pr["{\n"];
  pr["tBox *box = grid->box[bi];\n"];
  pr["int ijk;\n\n"];

  pr["forallpoints(box, ijk)\n"];
  pr["{\n"];
];

(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  prdecvar[{x[a]},		"i_x"];
  prdecvar[{gb[a,b]},    	"i_gb"];
  prdecvar[{K[a,b]},		"i_K"];
  prdecvar[{psi},		"i_psi"];
  prdecvar[{dpsiopsi[a]},	"i_dpsiopsi"];
  prdecvar[{ddpsiopsi[a,b]},	"i_ddpsiopsi"];
  prdecvar[{alpha},		"i_alpha"];
  prdecvar[{beta[a]},		"i_beta"];
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

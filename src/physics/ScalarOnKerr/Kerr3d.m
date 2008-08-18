(* Kerr3d.m 
   Wolfgang Tichy 9/2007 *)

(* Compute Kerr alpha, beta[a], g[a,b] from x,y,z and M=BHmass1 *)


(* variables *)
variables = { x[a], alpha, beta[a], g[a,b], gup[a,b], K[a,b], TrK,
	      Gam[a,b,c], dalpha[a], dbeta[a,b] }

(* compute in this order *)
tocompute = {
 
 rho == sqrt[ delta[a,b] x[a] x[b] ] ,
 aS[a] == If[a==1,Sx,If[a==2,Sy,Sz]],
 Cif == (M!=0.0),
   aS[a] == aS[a]/M,
 Cif == end,

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

 (* dl[a,b] = d_b l_a *)
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

 (* deriv of lapse *)
 dalpha[a] == -dH[a] / ((1 + 2 H)*sqrt[1 + 2 H]),

 (* deriv of shift *)
 dbeta[a,b] == ( 2/(1 + 2 H) - 4 H/(1 + 2 H)^2 ) dH[b] l[a] +
               2 H/(1 + 2 H) dl[a,b],

 (* extrinsic curv *) 
 K[a,b] == alpha ( l[a] dH[b] + l[b] dH[a] + H dl[a,b] + H dl[b,a] + 2 H^2 (l[a] l[c] dl[b,c] + l[b] l[c] dl[a,c]) + 2 H l[a] l[b] l[c] dH[c]),

 (* Tr K *)
 TrK == K[a,b] gup[a,b],

 (* dg[a,b,c] = d_c g_ab *)
 dg[a,b,c] == 2 dH[c] l[a] l[b] + 2 H dl[a,c] l[b] + 2 H l[a] dl[b,c],

 (* lower Christoffel *)
 Gamdo[a,b,c] == ( dg[a,b,c] + dg[a,c,b] - dg[b,c,a] )/2,

 (* standard upper Christoffel *)
 Gam[a,b,c] == gup[a,d] Gamdo[d,b,c]
 
}


(* symmetries *)

epsilon[1,2,3] := 1
epsilon[a_,b_,c_] := -epsilon[a,c,b] /; !OrderedQ[{b,c}]
epsilon[a_,b_,b_] := 0
epsilon[a_,b_,c_] := -epsilon[c,b,a] /; !OrderedQ[{a,c}]
epsilon[a_,b_,a_] := 0
epsilon[a_,b_,c_] := -epsilon[b,a,c] /; !OrderedQ[{a,b}]
epsilon[a_,a_,c_] := 0

g[a_,b_]     := g[b,a]   /; !OrderedQ[{a,b}]
gup[a_,b_]   := gup[b,a] /; !OrderedQ[{a,b}]
dg[a_,b_,c_] := dg[b,a,c]  /; !OrderedQ[{a,b}]
ahat[a_,b_]  := -ahat[b,a] /; !OrderedQ[{a,b}]
ahat[a_,a_]  := 0
Gamdo[a_,b_,c_] := Gamdo[a,c,b]  /; !OrderedQ[{b,c}]
Gam[a_,b_,c_]   := Gam[a,c,b]    /; !OrderedQ[{b,c}]



(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be in Cfunctionfile 
*)

CFunctionFile = "Kerr3d.c"

(* the head of the function *)
BeginCFunction[] := Module[{},

  pr["#include \"sgrid.h\"\n"];
  pr["#include \"ScalarOnKerr.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Log(x)     log((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];

  pr["\n\n\n"];

  pr["void Kerr3d(tGrid *grid, int i_x, int i_alpha, int i_beta, int i_g,
		  int i_K, int i_TrK, int i_gup,
                  int i_Gam, int i_dalpha, int i_dbeta)\n"];
  pr["{\n"];
  pr["int bi;\n"];
  pr["\n"];

  pr["double M = Getd(\"BHmass\");\n"];
  pr["double Sx = Getd(\"BHsx\");\n"];
  pr["double Sy = Getd(\"BHsy\");\n"];
  pr["double Sz = Getd(\"BHsz\");\n"];   
  pr["\n"];

  pr["forallboxes(grid,bi)\n"];
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
  prdecvar[{alpha},		"i_alpha"];
  prdecvar[{beta[a]},		"i_beta"];
  prdecvar[{g[a,b]},    	"i_g"];
  prdecvar[{K[a,b]},    	"i_K"];
  prdecvar[{TrK},               "i_TrK"];
  prdecvar[{gup[a,b]},		"i_gup"];
  prdecvar[{Gam[a,b,c]},	"i_Gam"];
  prdecvar[{dalpha[a]},         "i_dalpha"];
  prdecvar[{dbeta[a,b]},        "i_dbeta"];
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

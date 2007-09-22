(* Kerr.m 
   Wolfgang Tichy 9/12/2007 *)

(* Compute metric and Gamma from x,y,z and M=BHmass1 *)

(* variables *)
variables = { xp[a], g[a,b], gup[a,b], Gam[a,b,c], G[a] }

(* compute in this order *)
tocompute = {

 x[a] == If[a==0, 0, xp[a]],

 rho == sqrt[ x1*x1 + x2*x2 + x3*x3 ],
 aS[a] == If[a==1,Sx,If[a==2,Sy,If[a==3,Sz,0]]],

 Cif == (M!=0.0),
   aS[a] == aS[a]/M,
 Cif == end,

 adotx == aS1*x1 + aS2*x2 + aS3*x3,
 a2    == aS1*aS1 + aS2*aS2 + aS3*aS3,
 Omega == sqrt[(rho^2 - a2)^2 + 4 adotx^2],
 r == sqrt[(rho^2 - a2 + Omega)/2],

 dr[a] == 1/(2 r) ( 1 + (rho^2 - a2)/Omega) x[a] + adotx/(Omega r) aS[a],
 dr0 == 0,

 acrossx[a] == epsilon3D[a,b,c] aS[b] x[c],
 ahat[a,b] == epsilon3D[a,c,b] aS[c],

 H == M r^3/(r^4 + adotx^2),
 dH[a] == H/r (( 3 adotx^2 - r^4) dr[a] - 2 r adotx aS[a])/(r^4 + adotx^2),
 dH0 == 0,

 (* l_a = -k_a^{in}  in MTW p. 903 *)
 l[a] == (r^2 x[a] - r acrossx[a] + adotx aS[a])/(r (r^2 + a2)),
 l0   == +1,
 lup[a] == delta[a,b]*If[a==0,-1,1] l[b],

 (* dl[a,b] = d_b l_a *)
 dl[a,b] == (r^2 delta[a,b] - r ahat[a,b] + aS[a] aS[b])/(r^3 + a2 r) + (((a2 r^2 - r^4) x[a] + 2 r^3 acrossx[a] - (3 r^2 + a2) adotx aS[a]) dr[b])/(r^3 + a2 r)^2,
 dl0[b] == 0,
 dl10 == 0,
 dl20 == 0,
 dl30 == 0,

 (* now we can construct the Kerr metric *)

 (* the metric *)
 g[a,b] == delta[a,b]*If[a==0,-1,1] + 2 H l[a] l[b],

 (* the inverse metric *)
 gup[a,b] == delta[a,b]*If[a==0,-1,1] - 2 H lup[a] lup[b],

 (* dg[a,b,c] = d_c g_ab *)
 dg[a,b,c] == 2 dH[c] l[a] l[b] + 2 H dl[a,c] l[b] + 2 H l[a] dl[b,c],

 (* lower Christoffel *)
 Gamdo[a,b,c] == ( dg[a,b,c] + dg[a,c,b] - dg[b,c,a] )/2,

 (* standard upper Christoffel, and its contraction *)
 Gam[a,b,c] == gup[a,d] Gamdo[d,b,c],
 G[a] == gup[b,c] Gam[a,b,c]
}


(* symmetries *)
epsilon3D[0,b_,c_] := 0 (* dummy def. *)
epsilon3D[1,2,3] := 1
epsilon3D[a_,b_,c_] := -epsilon3D[a,c,b] /; !OrderedQ[{b,c}]
epsilon3D[a_,b_,b_] := 0
epsilon3D[a_,b_,c_] := -epsilon3D[c,b,a] /; !OrderedQ[{a,c}]
epsilon3D[a_,b_,a_] := 0
epsilon3D[a_,b_,c_] := -epsilon3D[b,a,c] /; !OrderedQ[{a,b}]
epsilon3D[a_,a_,c_] := 0

g[a_,b_]     := g[b,a]     /; !OrderedQ[{a,b}]
gup[a_,b_]   := gup[b,a]   /; !OrderedQ[{a,b}]
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

CFunctionFile = "Kerr.c"

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

  pr["void Kerr(tGrid *grid, int i_x, int i_g, int i_gup, int i_Gam, int i_G)\n"];
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

  pr["double *xp1 = box->v[i_x+0];\n"];
  pr["double *xp2 = box->v[i_x+1];\n"];
  pr["double *xp3 = box->v[i_x+2];\n"];
  prdecvar[{g[a,b]},    	"i_g"];
  prdecvar[{gup[a,b]},		"i_gup"];
  prdecvar[{Gam[a,b,c]},	"i_Gam"];
  prdecvar[{G[a]},		"i_G"];
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

(* use 4d tensors the default is 3 *)
TensorEquationsDim = 4;

(************************************************************************)
(* now we are ready to go *)

<< "../../Math/MathToC/TensorEquationsToC.m"

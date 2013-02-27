(* tensorrules.m *)
(* Bernd Bruegmann, 2/96, 10/02 & Wolfgang Tichy 2.2.2005 *)

(* Mathematica script, to be read by TensorEquationsToC.m:

   Define rules for tensor formulas like  
   - covariant derivatives
   - determinants
   - factoring out common factors like the metric

   Several entries here are obsolete, make current.
*)  


(* keep track of time *)
timerevents = {}
timerold = 0
timer[a_] := Module[{},
  timernew = TimeUsed[];
  timerevents = Append[timerevents, {a, timernew, timernew-timerold}];
  timerold = timernew;
]
prtimer[] := 
  Print["Timer:  ", PaddedForm[TableForm[timerevents], {8, 2}]]

timer["starting tensorrules.m"]



(* set number of dimensions
   the default is 3, resulting in 1,2,3
   if dimension variable is set to 4, start at 0, resulting in 0,1,2,3
   if dimension variable is set to 2, start at 1 end at 2, resulting in 1,2
*)
If [TensorEquationsDim === 4, imin = 0, imin = 1]
If [TensorEquationsDim === 2, imax = 2, imax = 3]


(**************************************************************************)
(* miscellaneous *)

(* this save 11 multiplies or so but is cumbersome *)
(*
  gg = Array[g, {3,3}];
  gginvdet = Inverse[gg] Det[gg];
  ginvdet[a_Integer,b_Integer] := gginvdet[[a,b]]

  gginv[a,b] == ginvdet[a,b],
  detg == g[1,c] gginv[1,c],
  ginv[a,b] == gginv[a,b] / detg,
*)

(* compute determinant and inverse *)
(* the number of dimensions can be 3d, 4d or 2d *)
matrixarray[g_] := Array[g, {imax-imin+1,imax-imin+1}, {imin,imin}]
matrixdet[g_] := Det[matrixarray[g]]
matrixinv[g_, a_Integer, b_Integer] :=
  Inverse[matrixarray[g]] [[a-imin+1, b-imin+1]]
matrixinvdet[g_, a_Integer, b_Integer] :=
  (Det[matrixarray[g]] Inverse[matrixarray[g]])[[a-imin+1, b-imin+1]]


(* Kronecker delta *)
delta[a_,b_] := 0 /; a != b
delta[a_,b_] := 1 /; a == b


(* The totally antisymmetric matrix in 3d *)
epsmatrix3d[a_Integer,b_Integer,c_Integer] := (b-a)(c-b)(c-a)/2
epsmatrix3d[a_,b_,c_] := -epsmatrix3d[b,a,c] /; !OrderedQ[{a,b}]
epsmatrix3d[a_,b_,c_] := -epsmatrix3d[a,c,b] /; !OrderedQ[{b,c}]


timer["after tensorrules.m"]

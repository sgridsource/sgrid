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

timer["starting tensorrules3d.m"]





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

matrixdet[g_] := Det[Array[g, {3, 3}]]
matrixinv[g_, a_Integer, b_Integer] := 
  Inverse[Array[g, {3, 3}]] [[a, b]]
matrixinvdet[g_, a_Integer, b_Integer] := 
  (Det[Array[g, {3, 3}]] Inverse[Array[g, {3, 3}]])[[a, b]]

delta[a_,b_] := 0 /; a != b
delta[a_,b_] := 1 /; a == b

(* The totally antisymmetric matrix *)
epsmatrix[a_Integer,b_Integer,c_Integer] := (b-a)(c-b)(c-a)/2
epsmatrix[a_,b_,c_] := -epsmatrix[b,a,c] /; !OrderedQ[{a,b}]
epsmatrix[a_,b_,c_] := -epsmatrix[a,c,b] /; !OrderedQ[{b,c}]


timer["after tensorrules3d.m"]

(* ontogrid.m *)
(* Wolfgang Tichy 2.2.2005 *)

(* Mathematica script, to be read by TensorEquations3dToC.m:

   Given equations with explicit indices, discretize for numerical grid:
   - put variables onto grid, g11 becomes g11[ijk] 
*)  

(***************************************************************************)
(* partial derivatives del and deldel *)
PartialDerivatives = {
  del[m_,0] :> 0,
  deldel[m_,n_,0] :> 0,
  del[m_, x_[ijk]]  :> ToExpression["d"<>ToString[Round[m]]<>ToString[x]<>"[ijk]"],
  deldel[m_, n_, x_[ijk]] :> ToExpression["d"<>ToString[Round[m]]<>"d"<>ToString[Round[n]]<>ToString[x]<>"[ijk]"]
} 



(***************************************************************************)
(* put variables on grid *)
gridvars = Join[vars];
putvarsgrid = Map[(# -> (#)[ijk])&, gridvars]
dif = all /. putvarsgrid;
If[False, prlist["dif = all /. putvarsgrid", dif]]

(* Derivatives del and deldel *)
dif = dif //. PartialDerivatives;
dif = DeleteCases[dif, True]

timer["after ontogrid.m"]

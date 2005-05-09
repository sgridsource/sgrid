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
  deldel[m_, n_, x_[ijk]] :> ToExpression["d"<>ToString[Round[m]]<>"d"<>ToString[Round[n]]<>ToString[x]<>"[ijk]"],
  OD[0,m_] :> 0,
  OD2[0,m_,n_] :> 0,
  OD[x_[ijk],m_]  :> ToExpression["d"<>ToString[x]<>ToString[Round[m]]<>"[ijk]"],
  OD2[x_[ijk],m_, n_] :> ToExpression["dd"<>ToString[x]<>ToString[Round[m]]<>ToString[Round[n]]<>"[ijk]"]
} 



(***************************************************************************)
(* put variables on grid *)
gridvars = Join[vars];
putvarsgrid = Map[(# -> (#)[ijk])&, gridvars]
dif = all /. putvarsgrid;
If[False, prlist["dif = all /. putvarsgrid", dif]]

(* Derivatives del, deldel and OD, OD2 *)
dif = dif //. PartialDerivatives;
dif = DeleteCases[dif, True]

timer["after ontogrid.m"]

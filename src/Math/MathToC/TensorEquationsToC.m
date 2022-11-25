(* TensorEquations3dToC.m *)
(* Bernd Bruegmann, 2/96, 10/02 & Wolfgang Tichy 2.2.2005 *)

(* Mathematica script:
   Given equations in tensor notation, write corresponding C code.
*)

(* Mathics is not looking relative to input files by default. *)
(* Add MathematicaToC directory to the Path, if we are running on Mathics *)
If[ValueQ[MathicsVersion],
  PrependTo[$Path, DirectoryName[$InputFileName]]
]

(*************************)
(* apply tensor formulas *)

<< tensorrules.m


(************************************************)
(* replace abstract indices by explicit indices *)

<< expandindices.m


(*************************)
(* put variables on grid *)

<< ontogrid.m

(************************************************)
(* optimize for best floating point performance *)

<< optimize.m


(****************)
(* write C code *)

<< writeC.m

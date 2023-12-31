(* expandindices3d.m *)
(* Wolfgang Tichy 2005, Bernd Bruegmann, 2/96, 10/02 *)

(* Mathematica script, to be read by TensorEquationsToC.m:

   Given equations in abstract tensor notation, expand indices:
   - expand contractions to sums with explicit indices
   - expand free indices by writing several equations with explicit indices
   - glue explicit indices to variable name, g[1,2] becomes g12
*)  




(*************************************************************************)
(* expand abstract indices *)

expcontraction = {
  s_[a___,b_Symbol,c___,d_Symbol,e___] \
    t_[f___,b_Symbol,g___] u_[h___,d_Symbol,i___] ->
    Sum[s[a,b,c,d,e] t[f,b,g] u[h,d,i], {b,1,3}, {d,1,3}],

  s_[a___,b_Symbol,c___] t_[d___,b_Symbol,e___] ->
    Sum[s[a,b,c] t[d,b,e], {b,1,3}]
}


(* careful: Union sorts, so apply to elements for which this is ok or wanted *)
expfreesort[x_,a_Symbol] := Union[Flatten[Table[x, {a,1,3}]]]
expfreesort[x_,{a_Symbol}] := expfreesort[x,a]
expfreesort[x_,{a_Symbol,b__}] := 
  Union[Flatten[Table[expfreesort[x,{b}], {a,1,3}]]]

expexpsign = -x_ :> x
expfreeexp[x_Symbol] := x
expfreeexp[x_[a__]] := 
  DeleteCases[Union[Flatten[expfreesort[x[a], {a}]/.expexpsign]], 0]

(* careful: when called on equations, LHS and RHS must have same syms *)
expequsign = -x_ == y_ :> x == -y
expfreeequ[x_Symbol==y_] := x==y
expfreeequ[delt[x_Symbol]==y_] := delt[x] == y 
expfreeequ[x_[a__]==y_] := 
  DeleteCases[Union[Flatten[expfreesort[x[a]==y,{a}]/.expequsign]], True]
expfreeequ[delt[x_[a__]]==y_] := 
  DeleteCases[Union[Flatten[expfreesort[delt[x[a]]==y,{a}]/.expequsign]], True]

(* expand free indices:
   for equations, expansion is based on lhs
   for expressions, only top level list is handled
*)
expfreeindices = {
  x_ :> Flatten[Map[expfreeequ,x]] /; !FreeQ[x,y_==z_],
  x_ :> Flatten[Map[expfreeexp,x]] /; FreeQ[x,y_==z_]
}



glue[a_,b_] := (* SequenceForm[a,b]; *) 
 ToExpression[ToString[a] <> ToString[b]]; 
glueall[x_] := x //. y_Symbol[a_,b___] :> glue[y,a][b] /. y_[] :> y
gluevariables[variables_] := Module[{varexpanded, vartoglue, varglued},
  varexpanded = Flatten[Map[expfreeexp, variables]];
  vartoglue = Cases[varexpanded, x_[a__]];
  varglued = Map[glueall, vartoglue]; 
  Table[vartoglue[[i]] -> varglued[[i]], {i, 1, Length[varglued]}]
]


(* utilities *)

prlist[a_,b_] := Block[{},
  Print[a, "\n"];
  Map[(Print[#,"\n"])&, b];
  Print["\n"];
];

colpow = {
  Power[x_,2] -> pow2[x], 
  Power[x_,-2] -> pow2inv[x]
}



(***************************************************************************)
(* transform equations into form that is to be differenced *)


(* substitutions *)
all = Join[tocompute]
prlist["all = Join[tocompute]", all]

(* do as often as there are contractions *)
all = ExpandAll[ExpandAll[all] //. expcontraction]
all = ExpandAll[ExpandAll[all] //. expcontraction]
all = ExpandAll[ExpandAll[all] //. expcontraction]
If[False, prlist["all = ExpandAll[all] /. expcontraction ...", all]]

(* expand free indices *)
all = all /. expfreeindices
If[False, prlist["all = all /. expfreeindices", all]]


(* various sets of variables *)

(* determine auxiliary variables that require a local definition:
   - these are all variables that appear on the lhs of an assignment
     but are not in the list of variables
   - treat the indices as abstract by using patterns
   - note that the level specification {2} also ensures that scalars do 
     not become patterns!
*)
(* old: auxvariables = Complement[Map[#[[1]]&, tocompute], variables]; *)

lhsvariables = Map[#[[1]]&, tocompute];
lhsvariables = DeleteCases[lhsvariables, Cif | Cinstruction];
patternofvariables = Map[(Pattern[#, Blank[]]) &, variables, {2}];
auxvariables = DeleteCases[lhsvariables, 
                           patternofvariables /. List -> Alternatives];
auxvariables = Union[auxvariables];   (* sort and remove duplicates *)

constvariables = {r, lr[a]} (* totally obsolete *)
gluevar = gluevariables[Join[variables, auxvariables, constvariables]];
vars = variables /. expfreeindices /. gluevar;
auxvars = auxvariables /. expfreeindices /. gluevar;
constvars = constvariables /. expfreeindices /. gluevar;

(* remove all vars from auxvars and constvars *)
auxvars = DeleteCases[auxvars, vars /. List -> Alternatives]
constvars = DeleteCases[constvars, vars /. List -> Alternatives]

Print["vars"];
Print[vars];
Print["auxvars"];
Print[auxvars];
Print["constvars"];
Print[constvars];

all = all /. gluevar
If[False, prlist["all /. gluevar", all]]


timer["after expandindices3d.m"]

(* expandindices.m *)
(*  Wolfgang Tichy 2005, Bernd Bruegmann, 2/96, 10/02 *)

(* Mathematica script, to be read by TensorEquationsToC.m:

   Given equations in abstract tensor notation, expand indices:
   - expand contractions to sums with explicit indices
   - expand free indices by writing several equations with explicit indices
   - glue explicit indices to variable name, g[1,2] becomes g12
*)  


(* set number of dimensions
   the default is 3, resulting in 1,2,3
   if dimension variable is set to 4, start at 0, resulting in 0,1,2,3
   if dimension variable is set to 2, start at 1 end at 2, resulting in 1,2
*)
If [TensorEquationsDim === 4, imin = 0, imin = 1]
If [TensorEquationsDim === 2, imax = 2, imax = 3]


(*************************************************************************)
(* expand abstract indices *)

expcontraction = {
  s_[a___,b_Symbol,c___,d_Symbol,e___] \
    t_[f___,b_Symbol,g___] u_[h___,d_Symbol,i___] ->
    Sum[s[a,b,c,d,e] t[f,b,g] u[h,d,i], {b,imin,imax}, {d,imin,imax}],

  s_[a___,b_Symbol,c___] t_[d___,b_Symbol,e___] ->
    Sum[s[a,b,c] t[d,b,e], {b,imin,imax}]
}


(* careful: Union sorts, so apply to elements for which this is ok or wanted *)
expfreesort[x_,a_Symbol] := Union[Flatten[Table[x, {a,imin,imax}]]]
expfreesort[x_,{a_Symbol}] := expfreesort[x,a]
expfreesort[x_,{a_Symbol,b__}] := 
  Union[Flatten[Table[expfreesort[x,{b}], {a,imin,imax}]]]

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

(* add localvariables to auxvariables *)
auxvariables = Flatten[Join[{localvariables}, auxvariables]];
auxvariables = Union[auxvariables];   (* sort and remove duplicates *)
locvariables = localvariables;
Clear[localvariables];
auxvariables = DeleteCases[auxvariables, localvariables];
localvariables = locvariables;
Clear[locvariables];

(* check whether constvariables is set, if not set it to {dummy} *)
convariables = constvariables;
Clear[constvariables];
If[convariables === constvariables, 
  constvariables = {dummy}, constvariables = convariables];
Clear[convariables];

(* proceed *)
gluevar = gluevariables[Join[variables, auxvariables, constvariables]];
vars = variables /. expfreeindices /. gluevar;
auxvars = auxvariables /. expfreeindices /. gluevar;
constvars = constvariables /. expfreeindices /. gluevar;

(* remove all vars from auxvars and constvars *)
auxvars = DeleteCases[auxvars, vars /. List -> Alternatives]
constvars = DeleteCases[constvars, vars /. List -> Alternatives]

(* sort and remove duplicates *)
vars = Union[vars];
auxvars = Union[auxvars];
constvars = Union[constvars];

Print["vars"];
Print[vars];
Print["constvars"];
Print[constvars];
Print["auxvars"];
Print[auxvars];

all = all /. gluevar
If[False, prlist["all /. gluevar", all]]


timer["after expandindices.m"]

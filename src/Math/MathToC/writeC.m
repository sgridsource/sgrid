(* writeC.m *)
(* Bernd Bruegmann, 2/96, 10/02 & Wolfgang Tichy 2.2.2005 *)

(* Mathematica script, to be read by TensorEquationsToC.m:

   Given equations, write C code:
   - write to file CFunctionFile
   - write header of function
   - write definitions of C variables
   - write equations as C assignments
   - write end of function
*)


(***************************************************************************)
(* write C code *)

Off[Part::partd]
Off[Part::pspec]

Print["Writing to ", CFunctionFile, "\n"];
DeleteFile[CFunctionFile];
filepointer = OpenAppend[CFunctionFile];
pr[x_] := Module[{s},
  WriteString[filepointer, x];
]

(* write comments *)
pr["/* "<>CFunctionFile<>" */\n"];
pr["/* Copyright (C) 2005-2022 Wolfgang Tichy, "]
If[$VersionNumber < 5.0,
  pr[ToString[Date[][[3]]]<>"."<>
     ToString[Date[][[2]]]<>"."<>
     ToString[Date[][[1]]]];,
  pr[ToString[DateList[][[3]]]<>"."<>
     ToString[DateList[][[2]]]<>"."<>
     ToString[DateList[][[1]]]];
]
pr[" */\n"];
pr["/* Produced with Mathematica "<>$Version<>"*/\n"];
pr["\n"];

(* write the function head *)
BeginCFunction[];


(***************************************************************************)
(* write definition of C variables *)

prdecvar[vars_, name_] := Module[{cvars, s, nf},
  cvars = vars /. expfreeindices /. gluevar;
  For[nf = 0, nf < Length[cvars], nf++,
    s = StringForm["double *`` = box->v[``+``];\n",
                    cvars[[nf+1]], name, nf];
    pr[s];
  ];
]

prdecvarname[vars_, name_] := Module[{cvars, s, nf},
  cvars = vars /. expfreeindices /. gluevar;
  s = StringForm["int index_`` = Ind(\"``\");\n", name, name];
  pr[s];
  For[nf = 0, nf < Length[cvars], nf++,
      (* box should not be hardwired *)
      s = StringForm["double *`` = box->v[index_`` + ``];\n",
                    cvars[[nf+1]], name, nf];
    pr[s];
  ];
]

prdecvl[vars_, name_] := Module[{cvars, s, nf},
  cvars = vars /. expfreeindices /. gluevar;
  For[nf = 0, nf < Length[cvars], nf++,
    s = StringForm["double *`` = vlldataptr(``, box, ``);\n",
                    cvars[[nf+1]], name, nf];
    pr[s];
  ];
]

prdecvlindices[vars_, name_] := Module[{cvars, s, nf},
  cvars = vars /. expfreeindices /. gluevar;
  For[nf = 0, nf < Length[cvars], nf++,
      s = StringForm["int index_`` = (``)->index[``];\n",
                    cvars[[nf+1]], name, nf];
    pr[s];
  ];
]


(* module to write var declarations for auxvars *)
writeCdef[vars_, auxvars_] := Module[{cvars, s, nf, nfields, nl},

  nfields = Length[auxvars];
  For[nf = 1, nf <= nfields, nf++,
    s = StringForm["double ``;\n", auxvars[[nf]]];
    pr[s];
  ];
  pr["\n\n\n"]
]

(* write var declarations for global vars in variabledeclarations[] *)
variabledeclarations[];
pr["\n"];

(* insert custom InitializationCommands[] *)
(* in InitializationCommands[] we can put an OpenMP loop and all
   auxvars are then thread local. *)
InitializationCommands[];

(* write var declarations for auxvars *)
writeCdef[vars, auxvars]


(***************************************************************************)
(* write equations *)

writeCassignCForm[e_] := Module[{x, nfields, nf, acc, one},

  If[e[[1]] =!= Cif && e[[1]] =!= Cinstruction,

    x = ToString[InputForm[e]];
    x = StringReplace[x, "[ijk]" -> "[[ijk]]"];
    x = ToExpression[x];
    x = x /. colpow;
    (* use N to get numerical expressions and get rid of 1. X and such *)
    acc = 20;
    one = N[1,acc];
    x = N[x,acc];
    x = x /. -1. xxx_ -> -xxx;
    x = x /. -one xxx_ -> -xxx;
    (* x = x /.  1. xxx_ -> xxx;
       x = x /. -one xxx_ -> -xxx;
       x = x /. -one xxx_ -> -xxx;
       x = x /. -one xxx_ -> -xxx;
       x = x /.  one xxx_ -> xxx;
       x = x /.  one xxx_ -> xxx; *)
    (* x = x /.  one xxx_ -> xxx; *)
    Print[ CForm[ x[[1]] ] ];

    PutAppend[CForm[ x[[1]] ], CFunctionFile];
    pr["=\n"];        (* to avoid assignment in CForm *)
    PutAppend[CForm[ x[[2]] ], CFunctionFile];
    pr[";\n\n"];
  ];

  If[e[[1]] === Cif,
    If[e[[2]] === end,
      pr["}\n"];
      pr[StringForm["/* if (``) */\n\n\n", lastif]]
    ];
    If[e[[2]] === else,
      pr["\n"];
      pr[StringForm["} else { /* if (!``) */\n\n", lastif]]
    ];
    If[e[[2]] =!= end && e[[2]] =!= else,
      pr["\n\n"];
      pr["/* conditional */\n"];
      pr[StringForm["if (``) {\n\n", e[[2]]]];
      lastif = e[[2]]
    ]
  ];

  (* Cinstruction: a piece of C-code *)
  If[e[[1]] === Cinstruction,
      pr["\n"];
      pr[StringForm["``\n\n", e[[2]]]];
  ];

]

(*************************************************************************)
(* HR: Mathics 4.0.0. does not implement CForm, so we use the code from
       BAMPS instead *)

(* BAMPS version of writeCassign *)

(* rules for a self-made MatlabForm, which does not exist in Mathematica *)
(* misnomer, adopted from writeMatlab *)

tomatlab[e_ == f_] := ToString[e] <> " = " <> tomatlab[f] <> ";"

tomatlab[r_Rational] :=
"(" <> ToString[N[Numerator[r]]] <> "/" <> ToString[N[Denominator[r]]] <> ")"

(* was "(" <> ToString[n] <> "/" <> ToString[m] <> ")"
   but for C, need e.g. -1/2 -> -1./2.
*)
(* works but strange: +-0.5000
   tovector[Rational[n_, m_]] := "("<>ToString[N[n/m,17]]<>")" /; n/m < 0 *)

tomatlab[Plus[e_, f__]] := "(" <> tomatlabplus[e, f] <> ")"
tomatlabplus[e_, f__] := tomatlab[e] <> " + " <> tomatlabplus[f]
tomatlabplus[e_] := tomatlab[e]

tomatlab[Times[-1, e__]] := "(-" <> tomatlab[Times[e]] <> ")"
tomatlab[Times[n_, e_]] :=
 tomatlab[n] <> "*" <> tomatlab[e] /; NumberQ[n]
tomatlab[Times[e_, f__]] := tomatlab[e] <> "*" <> tomatlab[Times[f]]

tomatlab[Power[e_, -1]] := "1/("<>tomatlab[e]<>")"
tomatlab[Power[e_, 2]] := tomatlab[e]<>"*"<>tomatlab[e]
tomatlab[Power[e_, 3]] := tomatlab[e]<>"*"<>tomatlab[e]<>"*"<>tomatlab[e]
tomatlab[Power[e_, p_]] := "Power(" <> tomatlab[e] <> "," <> tomatlab[p] <> ")"
tomatlabpower[e_, 0.5] := "sqrt("<>tomatlab[e]<>")"
tomatlabpower[e_,-0.5] := "(1/sqrt("<>tomatlab[e]<>"))"
tomatlabpower[e_,p_] := "Power(" <> tomatlab[e] <> "," <> tomatlab[p] <> ")"
tomatlab[Power[e_, p_Rational]] := tomatlabpower[e,N[p]]
(* "(" <> tomatlab[e] <> ".^" <> "(" <> tomatlab[p] <> "))" *)

(* unfinished, handle mathematica True *)
(*
tomatlab[Sign[e_]] := "((" <> tomatlab[e] <> " > 0) - (" <> tomatlab[e] <> " < 0))"
*)

(*fmf 02/19 implemented abs for rgrmhd module*)
tomatlab[Abs[e_]] := "fabs("<>tomatlab[e]<>")"

tomatlab[sin[e_]] := "sin(" <> tomatlab[e] <> ")"
tomatlab[cos[e_]] := "cos(" <> tomatlab[e] <> ")"
tomatlab[acos[e_]]:= "acos(" <> tomatlab[e] <> ")"
tomatlab[atan2[e_, p_]] := "atan2(" <> tomatlab[e] <> "," <> tomatlab[p] <> ")"
tomatlab[fabs[e_]]:= "fabs(" <> tomatlab[e] <> ")"
tomatlab[exp[e_]] := "exp(" <> tomatlab[e] <> ")"
tomatlab[tanh[e_]] := "tanh(" <> tomatlab[e] <> ")"
tomatlab[cosh[e_]] := "cosh(" <> tomatlab[e] <> ")"
tomatlab[log[e_]] := "log(" <> tomatlab[e] <> ")"
tomatlab[erf[e_]] := "erf(" <> tomatlab[e] <> ")"
tomatlab[sqrt[e_]] := "sqrt(" <> tomatlab[e] <> ")"

tomatlab[e_] := ToString[e]


(* write assignment *)
writeCassignBAMPS[e_] := Module[{},

  If[e[[1]] =!= Cif && e[[1]] =!= Cinstruction,
    (* info *)
    Print[CForm[e[[1]]]];

    (* write lhs *)
    pr[tomatlab[e]];
    pr["\n\n"];
  ];

  If[e[[1]] === Cif,
    If[e[[2]] === end,
      pr["}\n"];
      pr[StringForm["/* if (``) */\n\n\n", lastif]]
    ];
    If[e[[2]] === else,
      pr["\n"];
      pr[StringForm["} else { /* if (!``) */\n\n", lastif]]
    ];
    If[e[[2]] =!= end && e[[2]] =!= else,
      pr["\n\n"];
      pr["/* conditional */\n"];
      pr[StringForm["if (``) {\n\n", e[[2]]]];
      lastif = e[[2]]
    ]
  ];

  (* Cinstruction: a piece of C code *)
  If[e[[1]] === Cinstruction,
      pr["\n"];
      pr[StringForm["``\n\n", e[[2]]]];
  ];

];
(*************************************************************************)

(* choose writeCassign version *)
If[ValueQ[MathicsVersion],
   writeCassign = writeCassignBAMPS,
   writeCassign = writeCassignCForm
]

(* use writeCassign *)
Map[writeCassign, dif];

(***************************************************************************)
(* write the tail or end of the function *)

EndCFunction[];

(* count operations (Note: Mathics 5 has Count) *)

char = ReadList[CFunctionFile, Character];
nmul = Count[char, "*"];
ndiv = Count[char, "/"];
nsum = Count[char, "+"] + Count[char, "-"];
count = StringForm[
  "/* nvars = ``, n* = ``,  n/ = ``,  n+ = ``, n = ``, O = `` */\n", 
  Length[vars], nmul, ndiv, nsum, nmul+ndiv+nsum, 
  If [optimizeflag, 1, 0, 0]
];
Print[""];
Print[CFunctionFile];
Print[count];
pr["/* "<>CFunctionFile<>" */\n"];
pr[count];

Close[filepointer];

timer["after writeC.m"]
prtimer[]

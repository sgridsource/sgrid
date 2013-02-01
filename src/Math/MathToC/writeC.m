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

pr["/* "<>CFunctionFile<>" */\n"];
pr["/* Copyright (C) 2005-2008 Wolfgang Tichy, "<>
       ToString[Date[][[3]]]<>"."<>
       ToString[Date[][[2]]]<>"."<>
       ToString[Date[][[1]]]<>" */\n"];
pr["/* Produced with Mathematica */\n"];
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

writeCassign[e_] := Module[{x, nfields, nf, acc, one},

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
      pr[StringForm["`` \n\n", e[[2]]]];
  ];

]


Map[writeCassign, dif];


(***************************************************************************)
(* write the tail or end of the function *)

EndCFunction[];

(* count operations *)

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


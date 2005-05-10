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
pr["/* Copyright (C) 2005 Wolfgang Tichy & Bernd Bruegmann, "<>
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


(* write var declarations *)
writeCdef[vars_, auxvars_] := Module[{cvars, s, nf, nfields, nl},

  variabledeclarations[];
  pr["\n"];
  
  nfields = Length[auxvars];
  For[nf = 1, nf <= nfields, nf++,
    s = StringForm["double ``;\n", auxvars[[nf]]];
    pr[s];
  ];
  pr["\n\n\n"]
]

writeCdef[vars, auxvars]

(* insert custom InitializationCommands *)
InitializationCommands[];

(***************************************************************************)
(* write equations *)

writeCassign[e_] := Module[{x, nfields, nf},

  If[e[[1]] =!= Cif && e[[1]] =!= Cinstruction,

    x = ToString[InputForm[e]];
    x = StringReplace[x, "[ijk]" -> "[[ijk]]"];
    x = ToExpression[x];
    x = x /. colpow;
    Print[ CForm[ N[x[[1]],20] /. -1. xxx_ -> -xxx ] ];

    PutAppend[CForm[ N[x[[1]],20] /. -1. xxx_ -> -xxx ], CFunctionFile];
    pr["=\n"];        (* to avoid assignment in CForm *)
    PutAppend[CForm[ N[x[[2]],20] /. -1. xxx_ -> -xxx ], CFunctionFile];
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


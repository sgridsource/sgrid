(* Output function func[x] into file name in region={x,xmin,xmax} 
   using Npoints points. Ex: PlotData["lll", f[x,0] , {x,0,3}, 10] *)
PlotData[name_, func_, region_, Npoints_] := Module[{fp, i, x, y},
  fp = OpenWrite[name];
  For[i=0, i < Npoints, i++, 
    x = region[[2]] + i*(region[[3]] - region[[2]])/(Npoints-1);
    y = func /. region[[1]]->x;
    Print[CForm[N[x,17]], " \t\t\t", CForm[N[y,17]]];
    WriteString[fp, CForm[N[x,17]], " \t\t\t", CForm[N[y,17]]];
    WriteString[fp, "\n"];
  ];
  Close[fp];
];

(* Output two funcs into file name in region={x,xmin,xmax} 
   using Npoints points. Ex: PlotData["lll", f[x,0], g[x], {x,0,3}, 10] *)
PlotData2[name_, func1_, func2_, region_, Npoints_] := 
Module[{fp, i, x, y1, y2},
  fp = OpenWrite[name];
  For[i=0, i < Npoints, i++, 
    x = region[[2]] + i*(region[[3]] - region[[2]])/(Npoints-1);
    y1 = func1 /. region[[1]]->x;
    y2 = func2 /. region[[1]]->x;
    Print[CForm[N[x,17]], " \t", CForm[N[y1,17]], " \t", CForm[N[y2,17]]];
    WriteString[fp, CForm[N[x,17]], " \t", CForm[N[y1,17]], " \t", CForm[N[y2,17]]];
    WriteString[fp, "\n"];
  ];
  Close[fp];
];

(* out put a function on the points of SphericalDF if n2 and n3 are even *)
(* Example: 
   SphericalDFData["lll", Re[ YY[-2,2,2,t,p] ], t,p, 94,46]
*)
SphericalDFData[name_, func_, theta_, phi_, n2_, n3_] := Module[
  {fp, j,k, th,ph, dth,dph, y},
  dth = 2.0*Pi/n2;
  dph = 2.0*Pi/n3;

  fp = OpenWrite[name];
  WriteString[fp, "# SphericalDF:   ", "n2 = ",n2, "  n3 = ",n3, "\n"];
  WriteString[fp, "# gridtype = full", "\n"];
  WriteString[fp, "# ntheta = ", n2/2, "\n"];
  WriteString[fp, "# nphi   = ", n3, "\n"];
  WriteString[fp, "# iteration = 0", "\n"];
  WriteString[fp, "# time = 0", "\n"];

  For[k=0, k < n3, k++,
    For[j=0, j < n2/2, j++, 
      th = j*dth + 0.5*dth;
      ph = k*dph;
      y = func /. theta->th /. phi->ph;
      WriteString[fp, CForm[N[y,12]]," "];
    ];
    WriteString[fp, "\n"];
  ];
  Close[fp];
];

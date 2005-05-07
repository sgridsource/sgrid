(* optimize.m *)
(* Bernd Bruegmann, 2/96, 10/02, Wolfgang Tichy 2.2.2005 *)

(* Mathematica script, to be read by TensorEquations3dToC.m:

   Given equations, try to optimize for best floating point performance:
   - reduce number of multiplies and adds by factoring common terms
     this is quite effective because our typical equations often contain
     products like g_ij (v^i + w^j) that have been expanded by the script
   - much more could be done here
*)  

(***************************************************************************)
(* trivial optimization: factor as much as you can linearly *)


If[optimizeflag =!= False,

  (* to make formulas clearer, get out inverse metric first *)
  ginvvars = {ginv[a,b]} /. expfreeindices /. gluevar;
  collginv = x_ == y_ :> x == Collect[y, ginvvars];

  (* one of the deepest Mathematica insights I ever had,
     just try it without *)
  rule = a_ x_ + a_ y_ -> a (x + y);
  sule = a_ /; a < 0 -> minus (- a);
  tule = minus -> -1;

  dif = dif /. collginv;
  dif = dif //. rule;      (* sufficient for numbers with different signs *) 
  dif = dif /. sule //. rule /. tule;

]

timer["after optimize.m"]

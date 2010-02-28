(* from:
   http://demonstrations.wolfram.com/SpinWeightedSphericalHarmonics/ 
*)

(* Spin-Weighted Spherical Harmonics *)
\!\(Y[s_, l_, m_, \[Theta]_, \[Phi]_] := \(\((\(-1\))\)\^m\) 
      Simplify[\(\@\(\(\(\((l + m)\)!\) \(\((l - m)\)!\) \((2  l + 
                        1)\)\)\/\(\(\((l + s)\)!\) \(\((l - s)\)!\) 
                    4  \[Pi]\)\)\) \(\((Sin[\[Theta]\/2])\)\^\(2  
                l\)\) \(\[Sum]\_\(r = 0\)\%\(l - s\)\((Binomial[l - s, r] 
                Binomial[l + s, 
                  r + s - 
                    m] \(\((\(-1\))\)\^\(l - r - 
                      s\)\) \(\[ExponentialE]\^\(\[ImaginaryI]\ m\ \[Phi]\)\) \
\((Cot[\[Theta]\/2])\)\^\(2  r + s - m\))\)\), 
        Assumptions \[Rule] {\[Phi] \[Element] Reals, \[Theta] \[Element] 
              Reals}]; \)

(* Spin-Weighted Spherical Harmonics,
   without assuming that theta,phi are reals*)
\!\(YY[s_, l_, m_, \[Theta]_, \[Phi]_] := \(\((\(-1\))\)\^m\) 
      Simplify[\(\@\(\(\(\((l + m)\)!\) \(\((l - m)\)!\) \((2  l + 
                        1)\)\)\/\(\(\((l + s)\)!\) \(\((l - s)\)!\) 
                    4  \[Pi]\)\)\) \(\((Sin[\[Theta]\/2])\)\^\(2  
                l\)\) \(\[Sum]\_\(r = 0\)\%\(l - s\)\((Binomial[l - s, r] 
                Binomial[l + s, 
                  r + s - 
                    m] \(\((\(-1\))\)\^\(l - r - 
                      s\)\) \(\[ExponentialE]\^\(\[ImaginaryI]\ m\ \[Phi]\)\) \
\((Cot[\[Theta]\/2])\)\^\(2  r + s - m\))\)\) ]; \)

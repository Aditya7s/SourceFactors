(* ::Package:: *)

BeginPackage["KerrWKBFrequenciesNu`"]

rpeak::usage="rpeak[\[Mu],a] gives the position of the extrema of the WKB potential"
\[CapitalOmega]Rwkb::usage="\[CapitalOmega]Rwkb[\[Mu],a] gives the scaled real part of the WKB frequency"
\[CapitalOmega]Iwkb::usage="\[CapitalOmega]Iwkb[\[Mu],a] gives the scaled imaginary part of the WKB frequency"
\[Omega]WKB::usage="\[Omega]WKB[a,l,m,n] gives the WKB frequency in units M=1"


Begin["Private`"]

(*This polynomial determines the position of the peak of the potential*)

Poly[\[Mu]_,a_]:=2x^4*(x-3)^2 + 4 x^2 ((1-\[Mu]^2)*x^2 - 2x -3 (1-\[Mu]^2))*a^2 +(1-\[Mu]^2)*((2-\[Mu]^2)*x^2+ 2(2+\[Mu]^2)x + (2-\[Mu]^2))*a^4
	
(*These give the appropriate solutions for prograde and retrograde orbits*)

xp[\[Mu]_,a_]:=x/.(NSolve[Poly[\[Mu],a]==0,x][[5]]//First)

xr[\[Mu]_,a_]:=x/.(NSolve[Poly[\[Mu],a]==0,x][[6]]//First)

rpeak[\[Mu]_,a_]:=If[\[Mu]<0,xr[\[Mu],a],xp[\[Mu],a]]

\[CapitalOmega]Rwkb[\[Mu]_,a_]:=Module[{x},x=rpeak[\[Mu],a];
(Sqrt[2] (-1+x) )/Sqrt[4 x^2 (-3+x^2)-a^2 (-3+\[Mu]^2+x^2 (-3+\[Mu]^2)-2 x (1+\[Mu]^2))] 
]

\[CapitalOmega]Iwkb[\[Mu]_,a_]:=Module[{x,omR},
x=rpeak[\[Mu],a]; 
omR = (Sqrt[2](-1+x))/Sqrt[4 x^2 (-3+x^2)-a^2 (-3+\[Mu]^2+x^2 (-3+\[Mu]^2)-2 x (1+\[Mu]^2))];
(x^2-2x+a^2)*Sqrt[4*(6omR^2*x^2-1)+2a^2*omR^2*(3-\[Mu]^2)]/(2*x^4*omR-4*a*x*\[Mu]+a^2*x*omR(x*(3-\[Mu]^2)+2*(1+\[Mu]^2))+a^4*omR(1-\[Mu]^2))
]

\[Omega]WKB[a_,l_,m_,n_]:=(l+1/2)*\[CapitalOmega]Rwkb[m/(l+1/2),a]- I*(n+1/2)*\[CapitalOmega]Iwkb[m/(l+1/2),a]

End[]

EndPackage[]

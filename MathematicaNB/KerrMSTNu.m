(* ::Package:: *)

BeginPackage["KerrMSTNu`", {"KerrLeaverNu`"}]

nuappx::usage="nuappx[s,l,omega] gives a second order in omega expansion for \[Nu] which is pretty bad"

FindNu::usage="FindNu[a,s,m,omega,lambda,method:{{0,1},{0,1}},Nfrac:25,n0:-1] finds the generalized angular momentum. The method variable
	should either be a complex number or should be an interval in the complex plane (default).
	In the first case, it will use method as a seed for a root search, and in the second it will seek the root in the chosen
	interval."

FindNuInterval::usage="FindNuInverval[a,s,m,omega,lambda,interval:{{0,1},{0,1}},Nfrac:25,n0:-1] finds the generalized angular momentum 
   by starting in the center of the given interval in the complex plane. 
   There will always be a root in this interval or on the edges.
   Edges will occur when omega is real and small; this version may fail for these cases.
   Output is nu and two error measures. Nfrac gives the order of the continued fraction, and n0 the number of inversions."
   
FindNuSeed::usage="FindNuSeed[a,s,m,omega,lambda,nuseed:0.5+I*0.5,Nfrac:25,n0:-1] finds the generalized angular momentum with FindRoot 
	starting from a guess nuseed, with n0 inversions and Nfrac terms in the continued fraction."
   
fn::usage="fn[nu,a,s,m,omega,lambda,nf,Nfrac:20] produces the coefficients of the series solution, with 2nf+1 coefficients evenly around zero.
   The normalization is f0 = 1 and the coefficients are built up and down using CF methods to ensure convergence, provided nu is a sln to the CF.
   The order of the CF is Nfrac."
   
Knu::usage="Knu[\[Nu],a,s,m,\[Omega],\[Lambda],r:-1,nf:15,Nfrac:20] is a function used to build the solutions and scattering amplitudes.
    It is useful as a diagnostic output. It generates fn internally."
    
Btrans::usage="Btrans[nu,a,s,m,omega,lambda,method:{{0,1},{0,1}},n0:-1,nf:15,Nfrac:20] is the transmission amplitude. It generates nu and fn internally."

Binc::usage="Binc[nu,a,s,m,omega,lambda,method:{{0,1},{0,1}},n0:-1,nf:15,Nfrac:20,r:-1] is the incident amplitude. It generates nu and fn internally,
    and relies on Knu (hence the integer r, which can be anything)."
    
Bref::usage="Bref[nu,a,s,m,omega,lambda,method:{{0,1},{0,1}},n0:-1,nf:15,Nfrac:20] is the reflection amplitude. It generates nu and fn internally."


Begin["Private`"]

nuappx[s_,l_,omega_]:=l+1/(2l+1)*(-2-s^2/(l(l+1)) +( (l+1)^2-s^2)^2/((2l+1)(2l+2)(2l+3)) - (l^2 - s^2)^2/((2l-1)(2l)(2l+1)) )(2*omega)^2

FindNu[a_,s_,m_,omega_,lambda_,method_:{{0,1},{0,1}},Nfrac_:15,n0_:0]:=If[
	Length[method]==0,
	FindNuSeed[a,s,m,omega,lambda,method,Nfrac,n0],
	If[Dimensions[method]=={2,2},
		FindNuInterval[a,s,m,omega,lambda,method,Nfrac,n0],
		Null]
]

(*FindNu[a_,s_,m_,omega_,lambda_,l_:2,Nmax_:15,n0_:0,nuseed_:0.]:=Module[{\[Epsilon],\[Tau],\[Kappa],guess,alpha,beta,gamma,j,k,R,L,nu,nux,x,crosscheck},
(*Quantities used in def of recursion coeffs*)
\[Epsilon] = 2*omega;
\[Kappa]=Sqrt[1-a^2];
\[Tau] = (\[Epsilon] - m*a)/\[Kappa];
(*These are the coeffs of the recursion relation*)
alpha[\[Nu]_,n_]:=I*\[Epsilon]*\[Kappa]*(n+\[Nu]+1+s+I*\[Epsilon])*(n+\[Nu]+1+s-I*\[Epsilon])*(n+\[Nu]+1+I*\[Tau])/((n+\[Nu]+1)*(2n+2\[Nu]+3));
beta[\[Nu]_,n_]:=-lambda -s(s+1) + (n+\[Nu])*(n + \[Nu]+1)+ \[Epsilon]^2 +\[Epsilon](\[Epsilon]- m*a)+ \[Epsilon](\[Epsilon]- m*a)*(s^2+\[Epsilon]^2)/((n+\[Nu])*(n+\[Nu]+1));
gamma[\[Nu]_,n_]:= - I*\[Epsilon]*\[Kappa](n+\[Nu] - s + I*\[Epsilon])(n+\[Nu] - s - I*\[Epsilon])*(n+\[Nu]-I*\[Tau])/((n+\[Nu])*(2n+2\[Nu]-1));
(*Build the continued fraction Rn*)
R [nu_?NumericQ,n_]:=-gamma[nu,n]/(beta[nu,n]+ ContinuedFractionK[-alpha[nu,j]*gamma[nu,j+1],beta[nu,j+1],{j,n,Nmax-n}]);
(*Build the CF called Ln. Need a reversed CF here, since Ln is contructed by descending the index starting from n and going to -Infinity. To accomplish this, I use the regular CF in with index argument reversed.*)
L [nu_?NumericQ,n_]:= - alpha[nu,n]/(beta[nu,n]+ContinuedFractionK[-alpha[nu,-(k+1)]*gamma[nu,-k],beta[nu,-(k+1)],{k,-n,Nmax-n}]);
(*Now construct a guess using the approximate expression for nu*)
guess = If[nuseed\[Equal]0,nuappx[l,s,omega],nuseed];
nux= x/.FindRoot[beta[x,n0]+alpha[x,n0]*R[x,n0+1]+gamma[x,n0]*L[x,n0-1] \[Equal]0,{x,guess}];
crosscheck = R[nux,n0]*L[nux,n0-1]-1;
(*nux*)
{nux,beta[nux,n0]+alpha[nux,n0]*R[nux,n0+1]+gamma[nux,n0]*L[nux,n0-1]//Abs,crosscheck//Abs}
]
*)
FindNuInterval[a_,s_,m_,omega_,lambda_,interval_:{{0,1},{0,1}},Nfrac_:25,n0_:-1]:=Module[
{\[Epsilon],\[Tau],\[Kappa],guess,alpha,beta,gamma,j,k,R,L,nu,x,y,sln,nux,crosscheck},
(*Quantities used in def of recursion coeffs*)
\[Epsilon] = 2*omega;
\[Kappa]=Sqrt[1-a^2];
\[Tau] = (\[Epsilon] - m*a)/\[Kappa];
(*These are the coeffs of the recursion relation*)
alpha[\[Nu]_,n_]:=I*\[Epsilon]*\[Kappa]*(n+\[Nu]+1+s+I*\[Epsilon])*(n+\[Nu]+1+s-I*\[Epsilon])*(n+\[Nu]+1+I*\[Tau])/((n+\[Nu]+1)*(2n+2\[Nu]+3));
beta[\[Nu]_,n_]:=-lambda -s(s+1) + (n+\[Nu])*(n + \[Nu]+1)+ \[Epsilon]^2 +\[Epsilon](\[Epsilon]- m*a)+ \[Epsilon](\[Epsilon]- m*a)*(s^2+\[Epsilon]^2)/((n+\[Nu])*(n+\[Nu]+1));
gamma[\[Nu]_,n_]:= - I*\[Epsilon]*\[Kappa](n+\[Nu]-s+I*\[Epsilon])(n+\[Nu]-s-I*\[Epsilon])*(n+\[Nu]-I*\[Tau])/((n+\[Nu])*(2n+2\[Nu]-1));
(*Build the continued fraction Rn*)
R[nu_?NumericQ,n_]:=-gamma[nu,n]/(beta[nu,n]+ ContinuedFractionK[-alpha[nu,j]*gamma[nu,j+1],beta[nu,j+1],{j,n,Nfrac-n}]);
(*Build the CF called Ln. Need a reversed CF here, since Ln is contructed by descending the index starting from n and going to -Infinity. To accomplish this, I use the regular CF in with index argument reversed.*)
L[nu_?NumericQ,n_]:= - alpha[nu,n]/(beta[nu,n]+ContinuedFractionK[-alpha[nu,-(k+1)]*gamma[nu,-k],beta[nu,-(k+1)],{k,-n,Nfrac-n}]);
(*Now search the patch of the complex plane for nu*)
sln= {x,y}/.FindRoot[{Re[beta[x+I*y,n0]+alpha[x+I*y,n0]*R[x+I*y,n0+1]+gamma[x+I*y,n0]*L[x+I*y,n0-1]] ==0,Im[beta[x+I*y,n0]+alpha[x+I*y,n0]*R[x+I*y,n0+1]+gamma[x+I*y,n0]*L[x+I*y,n0-1]] ==0},{{x,(interval[[1,1]]+interval[[1,2]])/2,interval[[1,1]],interval[[1,2]]},{y,(interval[[2,1]]+interval[[2,2]])/2,interval[[2,1]],interval[[2,2]]}},AccuracyGoal->16,PrecisionGoal->16,Method->"Newton"];
nux=sln[[1]]+I*sln[[2]];
(*Make a crosscheck that should vanish*)
crosscheck = R[nux,n0]*L[nux,n0-1]-1;
(*Output two checks with nu*)
{nux,beta[nux,n0]+alpha[nux,n0]*R[nux,n0+1]+gamma[nux,n0]*L[nux,n0-1]//Abs,crosscheck//Abs}
]

FindNuSeed[a_,s_,m_,omega_,lambda_,nuseed_:0.5+I*0.5,Nfrac_:25,n0_:-1]:=Module[{\[Epsilon],\[Tau],\[Kappa],guess,alpha,beta,gamma,j,k,R,L,nu,nux,x,crosscheck},
(*Quantities used in def of recursion coeffs*)
\[Epsilon] = 2*omega;
\[Kappa]=Sqrt[1-a^2];
\[Tau] = (\[Epsilon] - m*a)/\[Kappa];
(*These are the coeffs of the recursion relation*)
alpha[\[Nu]_,n_]:=I*\[Epsilon]*\[Kappa]*(n+\[Nu]+1+s+I*\[Epsilon])*(n+\[Nu]+1+s-I*\[Epsilon])*(n+\[Nu]+1+I*\[Tau])/((n+\[Nu]+1)*(2n+2\[Nu]+3));
beta[\[Nu]_,n_]:=-lambda -s(s+1) + (n+\[Nu])*(n + \[Nu]+1)+ \[Epsilon]^2 +\[Epsilon](\[Epsilon]- m*a)+ \[Epsilon](\[Epsilon]- m*a)*(s^2+\[Epsilon]^2)/((n+\[Nu])*(n+\[Nu]+1));
gamma[\[Nu]_,n_]:= - I*\[Epsilon]*\[Kappa](n+\[Nu]-s+I*\[Epsilon])(n+\[Nu]-s-I*\[Epsilon])*(n+\[Nu]-I*\[Tau])/((n+\[Nu])*(2n+2\[Nu]-1));
(*Build the continued fraction Rn*)
R[nu_?NumericQ,n_]:=-gamma[nu,n]/(beta[nu,n]+ ContinuedFractionK[-alpha[nu,j]*gamma[nu,j+1],beta[nu,j+1],{j,n,Nfrac-n}]);
(*Build the CF called Ln. Need a reversed CF here, since Ln is contructed by descending the index starting from n and going to -Infinity. To accomplish this, I use the regular CF in with index argument reversed.*)
L[nu_?NumericQ,n_]:= - alpha[nu,n]/(beta[nu,n]+ContinuedFractionK[-alpha[nu,-(k+1)]*gamma[nu,-k],beta[nu,-(k+1)],{k,-n,Nfrac-n}]);
(*Now construct a guess using the approximate expression for nu*)
nux= x/.FindRoot[beta[x,n0]+alpha[x,n0]*R[x,n0+1]+gamma[x,n0]*L[x,n0-1] ==0,{x,nuseed}];
crosscheck = R[nux,n0]*L[nux,n0-1]-1;
(*nux*)
{nux,beta[nux,n0]+alpha[nux,n0]*R[nux,n0+1]+gamma[nux,n0]*L[nux,n0-1]//Abs,crosscheck//Abs}
]

(*Coefficients of series solution*)

fnBuildUp[nu_,a_,s_,m_,omega_,lambda_,nf_:15,Nfrac_:25]:=Module[{\[Epsilon],\[Tau],\[Kappa],guess,alpha,beta,gamma,R,f,j,k,p},
\[Epsilon] = 2*omega;
\[Kappa]=Sqrt[1-a^2];
\[Tau] = (\[Epsilon] - m*a)/\[Kappa];
(*These are the coeffs of the recursion relation*)
alpha[\[Nu]_,n_]:=I*\[Epsilon]*\[Kappa]*(n+\[Nu]+1+s+I*\[Epsilon])*(n+\[Nu]+1+s-I*\[Epsilon])*(n+\[Nu]+1+I*\[Tau])/((n+\[Nu]+1)*(2n+2\[Nu]+3));
beta[\[Nu]_,n_]:=-lambda -s(s+1) + (n+\[Nu])*(n + \[Nu]+1)+ \[Epsilon]^2 +\[Epsilon](\[Epsilon]- m*a)+ \[Epsilon](\[Epsilon]- m*a)*(s^2+\[Epsilon]^2)/((n+\[Nu])*(n+\[Nu]+1));
gamma[\[Nu]_,n_]:= - I*\[Epsilon]*\[Kappa](n+\[Nu] - s + I*\[Epsilon])(n+\[Nu] - s - I*\[Epsilon])*(n+\[Nu]-I*\[Tau])/((n+\[Nu])*(2n+2\[Nu]-1));
R [\[Nu]_?NumericQ,k_]:=-gamma[\[Nu],k]/(beta[\[Nu],k]+ ContinuedFractionK[-alpha[\[Nu],j]*gamma[\[Nu],j+1],beta[\[Nu],j+1],{j,k,Nfrac-k}]);
f = Range[nf+1]*0.0;
f[[1+0]]={0,1};
For[p=1,p<=nf,p++,f[[1+p]] = {p, f[[1+p-1,2]]*R[nu,p]};];
f
]

fnBuildDown[nu_,a_,s_,m_,omega_,lambda_,nf_:-15,Nfrac_:25]:=Module[{\[Epsilon],\[Tau],\[Kappa],guess,alpha,beta,gamma,L,f,j,k,p},
If[nf>0,Abort[],Null];
\[Epsilon] = 2*omega;
\[Kappa]=Sqrt[1-a^2];
\[Tau] = (\[Epsilon] - m*a)/\[Kappa];
(*These are the coeffs of the recursion relation*)
alpha[\[Nu]_,n_]:=I*\[Epsilon]*\[Kappa]*(n+\[Nu]+1+s+I*\[Epsilon])*(n+\[Nu]+1+s-I*\[Epsilon])*(n+\[Nu]+1+I*\[Tau])/((n+\[Nu]+1)*(2n+2\[Nu]+3));
beta[\[Nu]_,n_]:=-lambda -s(s+1) + (n+\[Nu])*(n + \[Nu]+1)+ \[Epsilon]^2 +\[Epsilon](\[Epsilon]- m*a)+ \[Epsilon](\[Epsilon]- m*a)*(s^2+\[Epsilon]^2)/((n+\[Nu])*(n+\[Nu]+1));
gamma[\[Nu]_,n_]:= - I*\[Epsilon]*\[Kappa](n+\[Nu] - s + I*\[Epsilon])(n+\[Nu] - s - I*\[Epsilon])*(n+\[Nu]-I*\[Tau])/((n+\[Nu])*(2n+2\[Nu]-1));
L [\[Nu]_?NumericQ,k_]:= -alpha[\[Nu],k]/(beta[\[Nu],k]+ContinuedFractionK[-alpha[\[Nu],-(j+1)]*gamma[\[Nu],-j],beta[\[Nu],-(j+1)],{j,-k,Nfrac-k}]);
f = Range[Abs[nf]+1]*0.0;
f[[1+0]]={0,1};
(*The indexing here looks odd, but the idea is the index of the array goes up, and the argument of L for f_p is p, and f_p = L(nu,p)f_(p+1) where f_(p+1) will be the previous entry in the array*)
For[p=1,p<= Abs[nf],p++,f[[1+p]] = {-p, f[[1+p-1,2]]*L[nu,-p]};];
(*So as to not repeat a coefficient I'll output from the second index on*)
f[[2;;]]//Reverse
]

fn[nu_,a_,s_,m_,omega_,lambda_,nf_:15,Nfrac_:25]:=Module[{fup,flow},
fup = fnBuildUp[nu,a,s,m,omega,lambda,nf,Nfrac];
flow = fnBuildDown[nu,a,s,m,omega,lambda,-nf,Nfrac];
Join[flow,fup]
]

(*Functions used to construct amplutudes*)

Knu[\[Nu]_,a_,s_,m_,\[Omega]_,\[Lambda]_,r_:-1,nf_:15,Nfrac_:25]:=Module[{\[Epsilon],\[Kappa],\[Tau],\[Epsilon]p,Nmax,ftab,n,Sum1,Sum2},
\[Epsilon] = 2*\[Omega];
\[Kappa]=Sqrt[1-a^2];
\[Tau] = (\[Epsilon] - m*a)/\[Kappa];
\[Epsilon]p = (\[Epsilon]+\[Tau])/2;
ftab = fn[\[Nu],a,s,m,\[Omega],\[Lambda],nf,Nfrac];
(*This arises from the D sum and starts at the start of ftab*)
Sum1 = Sum[(-1)^n*Pochhammer[\[Nu]+1+s - I*\[Epsilon],n]/((r-n)!*Pochhammer[r+2\[Nu]+2,n]*Pochhammer[\[Nu]+1-s + I*\[Epsilon],n])*ftab[[n+nf+1,2]],{n,-nf,r}];
(*This arises from the C sum and starts in the middle of ftab, running to the end*)
Sum2 = Sum[(-1)^n*(Gamma[n+r+2*\[Nu]+1]*Gamma[n+\[Nu]+1+s+I*\[Epsilon]]*Gamma[n+\[Nu]+1+I*\[Tau]])/((n-r)!Gamma[n+\[Nu]+1-s-I*\[Epsilon]]Gamma[n+\[Nu]+1-I*\[Tau]])*ftab[[n+nf+1,2]],{n,r,nf}];
Exp[I*\[Epsilon]*\[Kappa]]*(2*\[Epsilon]*\[Kappa])^(s-\[Nu]-r)*2^(-s)*I^r*Gamma[1-s-2*I*\[Epsilon]p]*Gamma[r+2\[Nu]+2]/(Gamma[r+\[Nu]+1-s+I*\[Epsilon]]Gamma[r+\[Nu]+1+I*\[Tau]]Gamma[r+\[Nu]+1+s+I*\[Epsilon]])*Sum2/Sum1
]

AnuP[ftab_,\[Nu]_,a_,s_,m_,\[Omega]_]:= Module[{\[Epsilon],\[Kappa],\[Tau]},
\[Epsilon] = 2*\[Omega];
\[Kappa]=Sqrt[1-a^2];
\[Tau] = (\[Epsilon] - m*a)/\[Kappa];
Exp[-Pi*\[Epsilon]/2]Exp[I*Pi*(\[Nu]+1-s)/2]*2^(-1+s-I*\[Epsilon])*Gamma[\[Nu]+1-s+I*\[Epsilon]]/Gamma[\[Nu]+1+s - I*\[Epsilon]]*Tr[ftab]
]

AnuM[ftab_,\[Nu]_,a_,s_,m_,\[Omega]_]:= Module[{\[Epsilon],\[Kappa],\[Tau],Nmax,n},
\[Epsilon] = 2*\[Omega];
\[Kappa]=Sqrt[1-a^2];
\[Tau] = (\[Epsilon] - m*a)/\[Kappa];
Nmax = (Length[ftab]-1)/2;
2^(-1-s+I*\[Epsilon])*Exp[-I*Pi*(\[Nu]+1+s)/2]Exp[-Pi*\[Epsilon]/2]*Sum[(-1)^n*Pochhammer[\[Nu]+1+s-I*\[Epsilon],n]/Pochhammer[\[Nu]+1-s+I*\[Epsilon],n]*ftab[[n+Nmax+1]],{n,-Nmax,Nmax}]
]

(*Amplitudes*)

Btrans[nu_,a_,s_,m_,omega_,lambda_,method_:{{0,1},{0,1}},n0_:-1,nf_:15,Nfrac_:25] :=Module[{\[Epsilon],\[Kappa],\[Tau],\[Epsilon]p,(*nu*),ftab},
\[Epsilon] = 2*omega;
\[Kappa]=Sqrt[1-a^2];
\[Tau] = (\[Epsilon] - m*a)/\[Kappa];
\[Epsilon]p = (\[Epsilon]+\[Tau])/2;
(*nu = FindNu[a,s,m,omega,lambda,method,Nfrac,n0][[1]];*)
ftab = fn[nu,a,s,m,omega,lambda,nf,Nfrac][[;;,2]];
(\[Epsilon]*\[Kappa]/omega)^(2s) Exp[I*\[Kappa]*\[Epsilon]p(1+2 Log[\[Kappa]]/(1+\[Kappa]))]*Tr[ftab]
]

Binc[nu_,a_,s_,m_,omega_,lambda_,method_:{{0,1},{0,1}},n0_:-1,nf_:15,Nfrac_:25,r_:-1]:= Module[{\[Epsilon],\[Kappa],\[Tau],(*nu*),ftab},
\[Epsilon] = 2*omega;
\[Kappa]=Sqrt[1-a^2];
\[Tau] = (\[Epsilon] - m*a)/\[Kappa];
(*nu = FindNu[a,s,m,omega,lambda,method,Nfrac,n0][[1]];*)
ftab = fn[nu,a,s,m,omega,lambda,nf,Nfrac][[;;,2]];
omega^(-1)*(Knu[nu,a,s,m,omega,lambda,r,nf,Nfrac]-I*Exp[-I*Pi*nu]*Sin[Pi(nu -s +I*\[Epsilon])]/Sin[Pi(nu+s-I*\[Epsilon])]*Knu[-nu-1,a,s,m,omega,lambda,r,nf,Nfrac])*AnuP[ftab,nu,a,s,m,omega]*Exp[-I*(\[Epsilon]*Log[\[Epsilon]]-(1 - \[Kappa])*\[Epsilon]/2)]
]

Bref[nu_,a_,s_,m_,omega_,lambda_,method_:{{0,1},{0,1}},n0_:-1,nf_:15,Nfrac_:25,r_:-1]:= Module[{\[Epsilon],\[Kappa],\[Tau],(*nu*),ftab},
\[Epsilon] = 2*omega;
\[Kappa]=Sqrt[1-a^2];
\[Tau] = (\[Epsilon] - m*a)/\[Kappa];
(*nu = FindNu[a,s,m,omega,lambda,method,Nfrac,n0][[1]];*)
ftab = fn[nu,a,s,m,omega,lambda,nf,Nfrac][[;;,2]];
omega^(-1-2*s)*(Knu[nu,a,s,m,omega,lambda,r,nf,Nfrac]+I*Exp[I*Pi*nu]*Knu[-nu-1,a,s,m,omega,lambda,r,nf,Nfrac])*AnuM[ftab,nu,a,s,m,omega]*Exp[I*(\[Epsilon]*Log[\[Epsilon]]-(1 - \[Kappa])*\[Epsilon]/2)]
]

End[]

EndPackage[]



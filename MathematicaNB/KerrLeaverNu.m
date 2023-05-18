(* ::Package:: *)

BeginPackage["KerrLeaverNu`", {"KerrWKBFrequenciesNu`"}]

Alm::usage="Alm[a,l,m,omega] gives the s=0 Spheroidal eigenvalue"
AlmS::usage="AlmS[s,l] = l(l+1) - s(s+1) is the Schwarzschild value of A"
sAlm::usage="sAlm[a,s,l,m,omega,Nmax:50] gives the spin-weighted spheroidal eigenvalue using Leaver's method without inversion"
sAlmExpansion::usage="sAlmExpansion[a,s,l,m,omega] gives the series expansion for the spin-weighted spheroidal harmonics up to sixth order in a*omega"
\[Lambda]Teuk::usage="\[Lambda]Teuk[a,m,omega,A] = A-2*a*m*\[Omega]+a^2*\[Omega]^2 is the radial eigenvalue of the Teuk eqn"

\[Delta]lmext::usage="\[Delta]lmext[l,m] gives the scalar-case nearly extremal quantity which controls the nature of the modes"
s\[Delta]lmext::usage="s\[Delta]lmext[s,l,m] gives the nearly-extremal quantity which controls the nature of the modes"

LeaverSearchInvS::usage="LeaverSearchInvS[s,l,m,n,\[Omega]seed:0.,Nmax:100] finds the QNM frequency for Schw using Leaver's method"
LeaverSearchInv::usage="LeaverSearchInv[spin,s,l,m,n,\[Omega]seed:0.,Nmax:400] finds the QNM frequency of Kerr for spin parameter spin(=a) using Leaver's method. 
	WKB is used to seed by default."

\[Omega]ZDM::usage="\[Omega]ZDM[a,s,l,m,n] gives the ZDM prediction for the QNM frequencies"
LeaverInvScaled::usage="LeaverInvScaled[a,s,l,m,nOvertone,nInv,guess:0.,Nmax:400] uses Leaver's method to search for \[Omega]t = (\[Omega] - m/2)/Sqrt[\[Epsilon]]. For guess = 0 it is seeded with the ZDM prediction."
LeaverSearchHS::usage="LeaverSearchHS[a,s,l,m,nOvertone,guess:0.,Nmax:400] uses a modified Leaver's method with nOvertone inversions to search for \[Omega]QNM and is best at high spin. 
	If guess=0, the ZDM prediction is used and may be bad"


Begin["Private`"]

(*Begin with the angular stuff*)

(*The s=0 case:*)

Alm[a_,l_,m_,omega_]:=N[SpheroidalEigenvalue[l,m,I*a*omega] - (a*omega)^2]

AlmS[s_,l_]:=l(l+1)-s(s+1)

(*I need a few functions for seeding things*)

Aseed[a_,s_,l_,m_,omega_]:=Alm[a,l,m,omega]-s(s+1);

sAlm[a_,s_,l_,m_,omega_,Nmax_:50]:=Module[{A,\[Gamma],\[Alpha],\[Beta],\[Alpha]t,\[Beta]t,\[Gamma]t,j,test,thisA,x},
\[Gamma] = a*omega;
\[Alpha]=Abs[m+s];
\[Beta]=Abs[m-s];
(*Important: I differ from Yanbei and Fackerell in that my A is the actual Alm eigenvalue*)
(*These three are the coefficients in the three term recursion relation for the Jacobi Polynomial expansion*)
\[Alpha]t[j_]:=-2*\[Gamma](j+\[Alpha]+1)(j+\[Beta]+1)(2j + \[Alpha]+\[Beta]+2s +2)/((2j + \[Alpha]+\[Beta]+2)(2j+\[Alpha]+\[Beta]+3));
\[Beta]t[j_,A_]:=A + s(s+1)+\[Gamma]^2 + 1/4 -(j +(\[Alpha]+\[Beta]+1)/2)^2 + 8s ( m s)\[Gamma]/((2j+\[Alpha]+\[Beta])(2j+\[Alpha]+\[Beta]+2));
\[Gamma]t[j_]:=2\[Gamma] j(j+\[Alpha]+\[Beta])(2j+\[Alpha]+\[Beta]-2s)/((2j+\[Alpha]+\[Beta]-1)(2j+\[Alpha]+\[Beta]));
(*This evaluates the continued fraction, which must be set to zero for the correct eigenvalue*)
test[A_?NumericQ]:=\[Beta]t[0,A]+ContinuedFractionK[-\[Alpha]t[j]*\[Gamma]t[j+1],\[Beta]t[j+1,A],{j,0,Nmax}];
thisA=x/.FindRoot[test[x]==0,{x,Aseed[a,s,l,m,omega]},MaxIterations->200];
{thisA,test[thisA]}
]

(* Alternative method using the expansion in Berti, Cardoso and Casals, this is used as a check. *)
(*Fan Zhang originally coded*)
(*Aaron implemented a simple fix for l=2, which otherwise gives singular results. *)

sAlmExpansion[a_,s_,l_,m_,omega_]:=Module[{\[Gamma],AlphaPlusBeta,AlphaMinusBeta,h,f0,f1,f2,f3,f4,f5,f6,ll},
\[Gamma] = a*omega;
AlphaPlusBeta=2Max[Abs[m],Abs[s]];
AlphaMinusBeta=2m s /Max[Abs[m],Abs[s]];
h[L_]:=(L^2-(1/4)AlphaPlusBeta^2)(L^2-(1/4)AlphaMinusBeta^2)(L^2-s^2)/2/(L-1/2)/L^3/(L+1/2);
f0=l(l+1)-s(s+1);
f1=-2 m s^2/l/(l+1);
f2=h[l+1]-h[l]-1;
f3=2h[l]m s^2/(l-1)/l^2/(l+1)-2h[l+1]m s^2/l/(l+1)^2/(l+2);
f4=m^2 s^4 (4h[l+1]/l^2/(l+1)^4/(l+2)^2-4h[l]/(l-1)^2/l^4/(l+1)^2)-(l+2)h[l+1]h[l+2]/2/(l+1)/(2l+3)+h[l+1]^2/2/(l+1)+h[l]h[l+1]/2/(l^2+l)-h[l]^2/2/l+(l-1)h[l-1]h[l]/(4l^2-2l);
f5= If[l==2,Limit[m^3 s^6 (8h[ll]/ll^6/(ll+3)^3/(ll-1)^3-8h[ll+1]/ll^3/(ll+1)^6/(ll+2)^3)+m s^2 h[ll](-h[ll+1](7ll^2+7ll+4)/ll^3/(ll+2)/(ll+1)^3/(ll-1)-h[ll-1](3ll-4)/ll^3/(ll+1)/(2ll-1)/(ll-2))+m s^2 ((3ll+7)h[ll+1]h[ll+2]/ll/(ll+1)^3/(ll+3)/(2ll+3)-3h[ll+1]^2/ll/(ll+1)^3/(ll+2)+3h[ll]^2/ll^3/(ll-1)/(ll+1)),ll->2],m^3 s^6 (8h[l]/l^6/(l+3)^3/(l-1)^3-8h[l+1]/l^3/(l+1)^6/(l+2)^3)+m s^2 h[l](-h[l+1](7l^2+7l+4)/l^3/(l+2)/(l+1)^3/(l-1)-h[l-1](3l-4)/l^3/(l+1)/(2l-1)/(l-2))+m s^2 ((3l+7)h[l+1]h[l+2]/l/(l+1)^3/(l+3)/(2l+3)-3h[l+1]^2/l/(l+1)^3/(l+2)+3h[l]^2/l^3/(l-1)/(l+1))];
f6=If[l==2,Limit[(16 m^4 s^8/ll^4/(ll+1)^4)(h[ll+1]/(ll+1)^4/(ll+2)^4-h[ll]/ll^4/(ll-1)^4)+(4m^2 s^4/ll^2/(ll+1)^2)(-(3ll^2+14ll+17)h[ll+1]h[ll+2]/(ll+1)^3/(ll+2)/(ll+3)^2/(2ll+3)+ 3h[ll+1]^2/(ll+1)^3/(ll+2)^2-3h[ll]^2/ll^3/(ll-1)^2)+(4m^2 s^4/ll^2/(ll+1)^2)((11ll^4+22ll^3+31ll^2+20ll+6)h[ll]h[ll+1]/ll^3/(ll-1)^2/(ll+1)^3/(ll+2)^2+(3ll^2-8ll+6)h[ll-1]h[ll]/ll^3/(ll-2)^2/(ll-1)/(2ll-1))+(h[ll+1]h[ll+2]/4/(ll+1)/(2ll+3)^2)((ll+3)h[ll+3]/3+((ll+2)/(ll+1))((ll+2)h[ll+2]-(7ll+10)h[ll+1]+(3ll^2+2ll-3)h[ll]/ll))+h[ll+1]^3/2/(ll+1)^2-h[ll]^3/2/ll^2+h[ll]h[ll+1]/4/ll^2/(ll+1)^2((2ll^2+4ll+3)h[ll]-(2ll^2+1)h[ll+1]-(ll^2-1)(3ll^2+4ll-2)h[ll-1]/(2ll-1)^2)+(h[ll-1]h[ll]/4/ll^2/(2ll-1)^2)((ll-1)(7ll-3)h[ll]-(ll-1)^2 h[ll-1]-ll (ll-2)h[ll-2]/3),ll->2],(16 m^4 s^8/l^4/(l+1)^4)(h[l+1]/(l+1)^4/(l+2)^4-h[l]/l^4/(l-1)^4)+(4m^2 s^4/l^2/(l+1)^2)(-(3l^2+14l+17)h[l+1]h[l+2]/(l+1)^3/(l+2)/(l+3)^2/(2l+3)+ 3h[l+1]^2/(l+1)^3/(l+2)^2-3h[l]^2/l^3/(l-1)^2)+(4m^2 s^4/l^2/(l+1)^2)((11l^4+22l^3+31l^2+20l+6)h[l]h[l+1]/l^3/(l-1)^2/(l+1)^3/(l+2)^2+(3l^2-8l+6)h[l-1]h[l]/l^3/(l-2)^2/(l-1)/(2l-1))+(h[l+1]h[l+2]/4/(l+1)/(2l+3)^2)((l+3)h[l+3]/3+((l+2)/(l+1))((l+2)h[l+2]-(7l+10)h[l+1]+(3l^2+2l-3)h[l]/l))+h[l+1]^3/2/(l+1)^2-h[l]^3/2/l^2+h[l]h[l+1]/4/l^2/(l+1)^2((2l^2+4l+3)h[l]-(2l^2+1)h[l+1]-(l^2-1)(3l^2+4l-2)h[l-1]/(2l-1)^2)+(h[l-1]h[l]/4/l^2/(2l-1)^2)((l-1)(7l-3)h[l]-(l-1)^2 h[l-1]-l (l-2)h[l-2]/3)];
f0 +f1 \[Gamma]^1+f2 \[Gamma]^2+f3 \[Gamma]^3+f4 \[Gamma]^4+f5 \[Gamma]^5+f6 \[Gamma]^6
]

\[Lambda]Teuk[a_,m_,omega_,A_]:=A-2*a*m*omega+a^2*omega^2

(*For nearly extremal holes*)

\[Delta]lmext[l_,m_]:=Sqrt[7/4*m^2 - Alm[1,l,m,m/2]-1/4]

s\[Delta]lmext[s_,l_,m_]:=Sqrt[7/4 m^2 -First[sAlm[1,s,l,m,m/2]]-(s+1/2)^2]

(*On to frequencies*)

LeaverSearchInvS[s_,l_,m_,n_,\[Omega]seed_:0.,Nmax_:100]:=Module[{\[Omega]s,\[Omega],\[Alpha]r,\[Beta]r,\[Gamma]r,Almx,j,lhs,rhs,test,\[Omega]x,x},
(*This function takes l,m,s,a and a seed guess for the QNM frequency and gives the result of a Leaver's method root find for \[Omega], using the guess. It also outputs the error in the CF equation, by evaluating the CF*)
(*The continued fraction is inverted at position n. If n=0, it will be the normal CF expansion*)
(*If a seed of 0 is provided, it will use the WKB formula to generate the seed*)
(*Gives size of CF expansion*)
\[Omega]s=If[\[Omega]seed==0,\[Omega]WKB[0.0,l,m,n],\[Omega]seed];
(*Set up CF*)
\[Alpha]r[j_,\[Omega]_]:=(1+j) (1+j-s-4 I \[Omega]);
\[Beta]r[j_,\[Omega]_]:=-1-Almx-2 j^2-s+6 I  \[Omega]+(24) \[Omega]^2+(\[Omega]) (2 I+8 \[Omega])+j (-2+16 I  \[Omega])/.Almx->l(l+1)-s(s+1);
\[Gamma]r[j_,\[Omega]_]:=(j-4 I \[Omega]) (j +s -4 I \[Omega]);
lhs[\[Omega]_]:=\[Beta]r[n,\[Omega]]+ContinuedFractionK[-\[Alpha]r[j- 1,\[Omega]]*\[Gamma]r[j,\[Omega]],\[Beta]r[j,\[Omega]],{j,n+1,Nmax}];
rhs[\[Omega]_]:=-ContinuedFractionK[-\[Alpha]r[n-j,\[Omega]]*\[Gamma]r[n-j+1,\[Omega]],\[Beta]r[n-j,\[Omega]],{j,1,n}];
(*This CF is equal to zero for QNMs*)
test[\[Omega]_?NumericQ]:=lhs[\[Omega]]-rhs[\[Omega]];
(*Seed the root finder with a decent guess*)
(*Use Find Root*)
\[Omega]x=x/.FindRoot[test[x]==0,{x,\[Omega]s}];
{\[Omega]x,test[\[Omega]x]}
]

(*For moderate and low spin Kerr*)

LeaverSearchInv[spin_,s_,l_,m_,n_,\[Omega]seed_:0.,Nmax_:400]:=Module[{a,\[Omega]s,b,\[Omega],\[Alpha]r,\[Beta]r,\[Gamma]r,Almx,j,lhs,rhs,test,\[Omega]x,x},
(*This function takes l,m,s,a and a seed guess for the QNM frequency and gives the result of a Leaver's method root find for \[Omega], using the guess. It also outputs the error in the CF equation, by evaluating the CF*)
(*The continued fraction is inverted at position n. If n=0, it will be the normal CF expansion*)
(*If a seed of 0 is provided, it will use the WKB formula to generate the seed*)
(*Gives size of CF expansion*)
(*To correct the units:*)
a = spin/2;
\[Omega]s=If[\[Omega]seed==0,\[Omega]WKB[l,m,n,a],\[Omega]seed];
(*Set up CF*)
b=(1 - 4 a^2)^(1/2);
\[Alpha]r[j_,\[Omega]_]:=((1+j) (I (2 a m-\[Omega])+b (1+j-s-I \[Omega])))/b;
\[Beta]r[j_,\[Omega]_]:=1/b (-(2 a m-\[Omega]) (I+2 I j+2 \[Omega])+b^2 \[Omega] (I+2 I j+2 \[Omega])-b (1+Almx+2 j+2 j^2+s-2 I \[Omega]-4 I j \[Omega]+2 a m \[Omega]-4 \[Omega]^2+a^2 \[Omega]^2))/.Almx->First[sAlm[a,s,l,m,\[Omega]]];
\[Gamma]r[j_,\[Omega]_]:=-(((I j+2 \[Omega]) (-2 a m+\[Omega]+b (I j+I s+\[Omega])))/b);
lhs[\[Omega]_]:=\[Beta]r[n,\[Omega]]+ContinuedFractionK[-\[Alpha]r[j- 1,\[Omega]]*\[Gamma]r[j,\[Omega]],\[Beta]r[j,\[Omega]],{j,n+1,Nmax}];
rhs[\[Omega]_]:=-ContinuedFractionK[-\[Alpha]r[n-j,\[Omega]]*\[Gamma]r[n-j+1,\[Omega]],\[Beta]r[n-j,\[Omega]],{j,1,n}];
(*This CF is equal to zero for QNMs*)
test[\[Omega]_?NumericQ]:=lhs[\[Omega]]-rhs[\[Omega]];
(*Seed the root finder with a decent guess*)
(*Use Find Root*)
\[Omega]x=x/.FindRoot[test[x]==0,{x,\[Omega]s}];
{\[Omega]x/2 ,test[\[Omega]x]}
]

(*For high spins this is best. First we get the value of the rescaled corotating frequency, then produce the final answer*)

\[Omega]ZDM[a_,s_,l_,m_,n_]:=m/2+Sqrt[1-a]/(Sqrt[2])*(-s\[Delta]lmext[s,l,m])-I*Sqrt[1-a]/(Sqrt[2])*(n+1/2)

LeaverInvScaled[a_,s_,l_,m_,nOvertone_,guess_:0.,Nmax_:400]:=Module[{\[Epsilon]\[Epsilon],\[Omega],\[Alpha],\[Alpha]r,\[Beta]r,\[Gamma]r,f,g,test,\[Epsilon],almx,bx,n,\[Omega]txpred,\[Omega]tx,\[Omega]t,InvertedTest},
\[Epsilon]\[Epsilon]=Sqrt[1-a];
\[Epsilon]=1-a;
\[Alpha]r[n_,\[CapitalOmega]t_]:=((1+n) (Sqrt[2-\[Epsilon]\[Epsilon]^2]-I m Sqrt[2-\[Epsilon]\[Epsilon]^2]+n Sqrt[2-\[Epsilon]\[Epsilon]^2]-s Sqrt[2-\[Epsilon]\[Epsilon]^2]-2 I \[CapitalOmega]t))/Sqrt[2-\[Epsilon]\[Epsilon]^2]+((1+n) \[Epsilon]\[Epsilon] (-I m-2 I Sqrt[2-\[Epsilon]\[Epsilon]^2] \[CapitalOmega]t))/Sqrt[2-\[Epsilon]\[Epsilon]^2];
If[s!=0,
\[Beta]r[n_,\[CapitalOmega]t_]:=-1-2 n^2-s-almx+m (-1+\[Epsilon]\[Epsilon]^2) (m+2 \[Epsilon]\[Epsilon] \[CapitalOmega]t)+I (2+\[Epsilon]\[Epsilon] Sqrt[2-\[Epsilon]\[Epsilon]^2]) (m+2 \[Epsilon]\[Epsilon] \[CapitalOmega]t)+1/4 (15+2 \[Epsilon]\[Epsilon]^2-\[Epsilon]\[Epsilon]^4+8 \[Epsilon]\[Epsilon] Sqrt[2-\[Epsilon]\[Epsilon]^2]) (m+2 \[Epsilon]\[Epsilon] \[CapitalOmega]t)^2+((m \[Epsilon]\[Epsilon]+2 \[CapitalOmega]t) (I+2 m+4 \[Epsilon]\[Epsilon] \[CapitalOmega]t))/Sqrt[2-\[Epsilon]\[Epsilon]^2]+n (-2+(2 I (m \[Epsilon]\[Epsilon]+2 \[CapitalOmega]t))/Sqrt[2-\[Epsilon]\[Epsilon]^2]+2 I (2+\[Epsilon]\[Epsilon] Sqrt[2-\[Epsilon]\[Epsilon]^2]) (m+2 \[Epsilon]\[Epsilon] \[CapitalOmega]t))/.almx->First[sAlm[a,s,l,m,m/2+\[CapitalOmega]t*Sqrt[\[Epsilon]]]],
\[Beta]r[n_,\[CapitalOmega]t_]:=-1-2 n^2-s-almx+m (-1+\[Epsilon]\[Epsilon]^2) (m+2 \[Epsilon]\[Epsilon] \[CapitalOmega]t)+I (2+\[Epsilon]\[Epsilon] Sqrt[2-\[Epsilon]\[Epsilon]^2]) (m+2 \[Epsilon]\[Epsilon] \[CapitalOmega]t)+1/4 (15+2 \[Epsilon]\[Epsilon]^2-\[Epsilon]\[Epsilon]^4+8 \[Epsilon]\[Epsilon] Sqrt[2-\[Epsilon]\[Epsilon]^2]) (m+2 \[Epsilon]\[Epsilon] \[CapitalOmega]t)^2+((m \[Epsilon]\[Epsilon]+2 \[CapitalOmega]t) (I+2 m+4 \[Epsilon]\[Epsilon] \[CapitalOmega]t))/Sqrt[2-\[Epsilon]\[Epsilon]^2]+n (-2+(2 I (m \[Epsilon]\[Epsilon]+2 \[CapitalOmega]t))/Sqrt[2-\[Epsilon]\[Epsilon]^2]+2 I (2+\[Epsilon]\[Epsilon] Sqrt[2-\[Epsilon]\[Epsilon]^2]) (m+2 \[Epsilon]\[Epsilon] \[CapitalOmega]t))/.almx->Alm[a,l,m,m/2+\[CapitalOmega]t*Sqrt[\[Epsilon]]]];
\[Gamma]r[n_,\[CapitalOmega]t_]:=-(1/Sqrt[2-\[Epsilon]\[Epsilon]^2])(2 m^2 Sqrt[2-\[Epsilon]\[Epsilon]^2]+3 I m n Sqrt[2-\[Epsilon]\[Epsilon]^2]-n^2 Sqrt[2-\[Epsilon]\[Epsilon]^2]+2 I m s Sqrt[2-\[Epsilon]\[Epsilon]^2]-n s Sqrt[2-\[Epsilon]\[Epsilon]^2]+4 m \[CapitalOmega]t+2 I n \[CapitalOmega]t)-1/Sqrt[2-\[Epsilon]\[Epsilon]^2] \[Epsilon]\[Epsilon] (2 m^2+I m n+8 m Sqrt[2-\[Epsilon]\[Epsilon]^2] \[CapitalOmega]t+6 I n Sqrt[2-\[Epsilon]\[Epsilon]^2] \[CapitalOmega]t+4 I s Sqrt[2-\[Epsilon]\[Epsilon]^2] \[CapitalOmega]t+8 \[CapitalOmega]t^2)-(\[Epsilon]\[Epsilon]^2 (4 m \[CapitalOmega]t+8 Sqrt[2-\[Epsilon]\[Epsilon]^2] \[CapitalOmega]t^2))/Sqrt[2-\[Epsilon]\[Epsilon]^2];
f[n_,\[Omega]t_]:=-\[Alpha]r[n-1,\[Omega]t]\[Gamma]r[n,\[Omega]t];
g[n_,\[Omega]t_]:=\[Beta]r[n,\[Omega]t];
(*This gives the analytically found, near extremal \[Omega]-m/2, when found using Yanbei's formula/Detweiler's expansion, and after scaling by \[Epsilon] = 1-a. It is used to seed the find root algorithm if there is no guess, guess=0. It depends on the overtone number n so that I can search out higher overtone modes.*)
\[Omega]txpred=If[s!=0,(-s\[Delta]lmext[l,m,s])/Sqrt[2]-I (nOvertone+1/2)/Sqrt[2.],(-\[Delta]lmext[l,m])/Sqrt[2]-I (nOvertone+1/2)/Sqrt[2.]];
InvertedTest[\[Omega]t_?NumericQ,kk_]:=g[kk,\[Omega]t]+ContinuedFractionK[f[n,\[Omega]t],g[n,\[Omega]t],{n,kk+1,Nmax}]  + Module[{ii,tempfrac},tempfrac=0;For[ii=1,ii<=kk,ii++,
tempfrac=f[ii,\[Omega]t]/(g[ii-1,\[Omega]t]+tempfrac)];
tempfrac];
(*If guess is given as 0, then the analytically calculated scaling is used, otherwise the guess is used*)
If[guess==0,
\[Omega]tx=\[Omega]t/.FindRoot[InvertedTest[\[Omega]t,nOvertone]==0,{\[Omega]t,\[Omega]txpred}],
\[Omega]tx=\[Omega]t/.FindRoot[InvertedTest[\[Omega]t,nOvertone]==0,{\[Omega]t,(guess-m/2)/Sqrt[\[Epsilon]]}]
];
(*This gives back \[Omega]t and an error, based on the evaluating the continued fractions expansion*)
{\[Omega]tx,InvertedTest[\[Omega]tx,nOvertone]}
]

(*This function gives the final answer for the high spin search*)
LeaverSearchHS[a_,s_,l_,m_,nOvertone_,guess_:0.,Nmax_:400]:=Module[{omegatout},
omegatout=LeaverInvScaled[a,s,l,m,nOvertone,guess,Nmax];
{m/2+ omegatout[[1]]*Sqrt[1-a],omegatout[[2]]}
]

(*Remember stick to the convention a, s, l, m, omega, extras, fixed quantities*)
End[]

EndPackage[]




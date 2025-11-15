(* ::Package:: *)

(* ::Title:: *)
(*External Data: Puncture and Effective Source*)


Clear["Global`*"]
NotebookDirectory[]
SetDirectory[NotebookDirectory[]]
Print["Start Notebook"]
LaunchKernels[];
$ProcessorCount


(* ::Section:: *)
(*Set Parameters*)


\[Lambda]=2rh; (*Hyperboloidal Length Scale*)
CoordFlagValue=2; (*Coordinate Map Type*)


q=1;  (*Scalar Charge*)
M=1;  (*BH mass*)

r0=10 M; (*Particle orbital radius*)

m=2; (*m-mode*)
nbarmax=2; (*Puncture Order*)
NN=50; (*Spectral Resolution for Puncture/Effective Source*)
Prec=4 MachinePrecision; (*Internal digits for Puncture/Effective Source*)


(* ::Chapter::Closed:: *)
(*Preamble*)


(* ::Section:: *)
(*DumpSave load*)


(*Print["Load Patrick"]*)
With[{dumpsaveFilePath = "philmnSing_up_to_n14.mx"},
	If[FileExistsQ[dumpsaveFilePath],Get[dumpsaveFilePath],Print["No DumpSave to load"]]]


(* ::Section::Closed:: *)
(*Spectral Routines*)


(*Print["Spec Routines"]*)


(* ::Subsection:: *)
(*Numerical Grids*)


Grid$Gauss[i_,Nz_]:=Cos[\[Pi] (i+1/2) /(Nz+1)]
Grid$Lobatto[i_,Nz_]:=Cos[\[Pi]*i /Nz]
Grid$RadauRHS[i_,Nz_]:=Cos[2 \[Pi]*i /(2 Nz+1)]
Grid$RadauLHS[i_,Nz_]:=Cos[\[Pi]-2 \[Pi]*i /(2Nz+1)]


(* ::Subsection:: *)
(*Chebyshev Coefficients*)


ChebGaussCoefficients[f_, Prec_]:=Module[
{nf, Nf,c,x},
nf=Length@f;
Nf=nf-1;

x[j_]:=N[ Grid$Gauss[j,Nf] ,Prec];

c=Table[  2 Sum[ f[[j+1]] * ChebyshevT[ m,x[j] ],{j,0,Nf}]  /nf        ,{m,0,Nf}]
];

ChebLobattoCoefficients[f_, Prec_]:=Module[
{nf, Nf,c,x},
nf=Length@f;
Nf=nf-1;

x[j_]:=N[ Grid$Lobatto[j,Nf] ,Prec];

c=Table[  ( 2-KroneckerDelta[m,Nf]) /(2 Nf) (  f[[1]]+(-1)^m f[[nf]] + 2 Sum[f[[j+1]]*ChebyshevT[m,x[j]],{j,1,Nf-1}]  )        ,{m,0,Nf}]
];

ChebRadauRHSCoefficients[f_, Prec_]:=Module[
{nf, Nf,c,x},
nf=Length@f;
Nf=nf-1;

x[j_]:=N[ Grid$RadauRHS[j,Nf] ,Prec];

c=Table[  4/(2*Nf+1)*( f[[1]]/2  +  Sum[f[[j+1]]*ChebyshevT[m,x[j]],{j,1,Nf}]  )        ,{m,0,Nf}]
];

ChebRadauLHSCoefficients[f_, Prec_]:=Module[
{nf, Nf,c,x},
nf=Length@f;
Nf=nf-1;

x[j_]:=N[ Grid$RadauLHS[j,Nf] ,Prec];

c=Table[  4/(2*Nf+1)*( (-1)^m f[[nf]]/2  +  Sum[f[[j+1]]*ChebyshevT[m,x[j]],{j,1,Nf}]  )        ,{m,0,Nf}]
];


(* ::Subsection:: *)
(*Spectral Differentiation matrices*)


DiffMatrixLobatto[a_,b_,Nz_,Prec_]:=Module[
{k,i,j,dx, \[CapitalDelta],Dy},

k[i_]:=If[i*(i-Nz)==0, 2,1];
dx[i_,j_]:=If[i!=j,
k[i]*(-1)^(i-j)/(k[j]*(Grid$Lobatto[i,Nz]-Grid$Lobatto[j,Nz])),
0
];
\[CapitalDelta]=b-a;
Dy=N[2/\[CapitalDelta] * Table[dx[i,j], {i,0,Nz}, {j,0,Nz}],Prec ];

For[i=0, i<=Nz, i++,
Dy[[i+1,i+1]]= -Sum[ Dy[[i+1,j+1]], {j,0,Nz} ];
];
Dy
]


(* ::Subsection:: *)
(*Chebyshev Interpolation*)


ChebInterpolation[c_,a_,b_,x_]:=Module[
{nc, Nc},
nc=Length@c;
Nc=nc-1;

y=(2 x-(a+b))/(b-a);

c[[1]]/2+Sum[c[[i+1]]*ChebyshevT[i,y],{i,1,Nc}]
];

ChebyLobattoInterpolationVector[xl_,Nl_,x\[Alpha]_,Prec]:=Module[
{\[Alpha],l},
Table[
If[l==0,
1/(2Nl)+Sum[(2-KroneckerDelta[j,Nl])ChebyshevT[j,x\[Alpha]]/(2Nl),{j,1,Nl}],
If[l==Nl,
1/(2Nl)+Sum[(2-KroneckerDelta[j,Nl])(-1)^j ChebyshevT[j,x\[Alpha]]/(2Nl),{j,1,Nl}],
1/(Nl)+Sum[(2-KroneckerDelta[j,Nl]) ChebyshevT[j,xl[[l+1]]]ChebyshevT[j,x\[Alpha]]/(Nl),{j,1,Nl}]
]
],
{l,0,Nl}
]

];

ChebyLobattoInterpolationMatrix[xl_,Nl_,x\[Alpha]_,N\[Alpha]_,Prec]:=Module[
{\[Alpha],l},
Table[
If[l==0,
1/(2Nl)+Sum[(2-KroneckerDelta[j,Nl])ChebyshevT[j,x\[Alpha][[\[Alpha]+1]]]/(2Nl),{j,1,Nl}],
If[l==Nl,
1/(2Nl)+Sum[(2-KroneckerDelta[j,Nl])(-1)^j ChebyshevT[j,x\[Alpha][[\[Alpha]+1]]]/(2Nl),{j,1,Nl}],
1/(Nl)+Sum[(2-KroneckerDelta[j,Nl]) ChebyshevT[j,xl[[l+1]]]ChebyshevT[j,x\[Alpha][[\[Alpha]+1]]]/(Nl),{j,1,Nl}]
]
],
{\[Alpha],0,N\[Alpha]},{l,0,Nl}
]

]


(* ::Subsection::Closed:: *)
(*Chebyshev Integration*)


DefiniteIntegral[c_,a_,b_]:=Module[
{kf,nc,Nc, Int},
nc=Length@c;
Nc=nc-1;
kf=Floor[Nc/2];
Int=c[[1]]- 2Sum[c[[2 k +1]]/(4 k ^2 -1),{k,1,kf}];
(b-a)*Int/2
];

ChebyLobattoHilbertIntegration[z_,Nz_,Prec_]:=Module[
{i,j},
Table[
If[i!=j,
0,
If[i==0 || i==Nz,
(1/Nz - Sum[(2-KroneckerDelta[2k,Nz])/(Nz(4k^2-1)),{k,1,Floor[Nz/2]}]),
(2/Nz - 2Sum[(2-KroneckerDelta[2k,Nz])*ChebyshevT[2k,z[[i+1]]]/(Nz(4k^2-1)),{k,1,Floor[Nz/2]}])
]
],
 {i,0,Nz}, {j,0,Nz}
]
]


(* ::Chapter::Closed:: *)
(*\[CapitalPhi]_m in BH frame*)


(*Print["punc Routines"]*)


(* ::Section:: *)
(*\[ScriptCapitalI](\!\(\*OverscriptBox[\(n\), \(_\)]\),\!\(\*OverscriptBox[\(l\), \(_\)]\),\!\(\*OverscriptBox[\(m\), \(_\)]\),m)*)


(* ::Input::Initialization:: *)
Clear[my\[ScriptCapitalI]Num]
my\[ScriptCapitalI]Num[\[Rho]_?NumericQ][zc_][nbar_?NumericQ,0,0,m_?NumericQ] := my\[ScriptCapitalI]Num[\[Rho]][zc][nbar,0,0,m] = 
Block[{cBlock = \[Rho] Sqrt[\[Rho]^2+zc^2], zBlock = (\[Rho]^2+zc^2/2)/(\[Rho] Sqrt[\[Rho]^2+zc^2]), mRegBlock = If[1+nbar/2 <=0,-Abs[m],Abs[m]]},  (2 \[Pi] cBlock^(nbar/2) LegendreP[nbar/2,mRegBlock,3,zBlock]/(-1)^mRegBlock)/Pochhammer[1+nbar/2,mRegBlock]]
(* my\[ScriptCapitalI]Num[\[Rho]_?NumericQ][zc_?NumericQ][n_?NumericQ,1,0,m_?NumericQ] := my\[ScriptCapitalI]Num[\[Rho]][zc][n,1,0,m] = zc (my\[ScriptCapitalI]Num[\[Rho]][zc][n-1,0,0,m-1/2] - my\[ScriptCapitalI]Num[\[Rho]][zc][n-1,0,0,m+1/2])/(2 I); *)
(*my\[ScriptCapitalI]Num[\[Rho]_?NumericQ][zc_?NumericQ][n_?NumericQ,l_?NumericQ,lm1_?NumericQ,m_?NumericQ]  := (my\[ScriptCapitalI]Num[\[Rho]][zc][n,l,l-1,m] = (-1)^(l-1) Factorial2[2 l-3] (2 l-1) \[Rho]^(l-1) my\[ScriptCapitalI]Num[\[Rho]][zc][n-l+1,1,0,m]) /; lm1 \[Equal] l-1*)
my\[ScriptCapitalI]Num[\[Rho]_?NumericQ][zc_?NumericQ][nbar_?NumericQ,lbar_?NumericQ,lm1_?NumericQ,mbar_?NumericQ]  := 0 /; lm1 == lbar-1
my\[ScriptCapitalI]Num[\[Rho]_?NumericQ][zc_?NumericQ][nbar_?NumericQ,lbar_?NumericQ,lbar_?NumericQ,mbar_?NumericQ] := my\[ScriptCapitalI]Num[\[Rho]][zc][nbar,lbar,lbar,mbar]= (-1)^lbar Factorial2[2 lbar-1] \[Rho]^lbar my\[ScriptCapitalI]Num[\[Rho]][zc][nbar-lbar,0,0,mbar];
my\[ScriptCapitalI]Num[\[Rho]_?NumericQ][zc_?NumericQ][n_?NumericQ,lbar_?NumericQ,mbar_?NumericQ,m_?NumericQ] := (my\[ScriptCapitalI]Num[\[Rho]][zc][n,lbar,mbar,m] = (If[mbar+2 > lbar,0,-my\[ScriptCapitalI]Num[\[Rho]][zc][n,lbar,mbar+2,m]] + If[mbar+2 > lbar-2,0,my\[ScriptCapitalI]Num[\[Rho]][zc][n,lbar-2,mbar+2,m]] + 
If[mbar > lbar-2,0,(lbar+mbar) (lbar+mbar-1) my\[ScriptCapitalI]Num[\[Rho]][zc][n,lbar-2,mbar,m]])/((lbar-mbar-1) (lbar-mbar))) /; mbar < lbar-1


(* ::Section:: *)
(*Subscript[\[CapitalPhi], m]*)


(* ::Text:: *)
(*NOTE: clm is using floating precision. This makes computations a bit faster.*)
(*I checked, with nbar=6, that the scalar, its derivatives and the effective source agree to machine precision to each other, when comparing floating-point precision with arbitrary precision, taking 100 digits into account (specifically, the maximum difference between floating and arbitrary is (~10)^-17 for scalar & derivatives, and ~(10^-15) for effective source).*)


(* ::Input::Initialization:: *)
Clear[clm]
clm[l_,m_] := Sqrt[(2 l+1)/(4 \[Pi]) Factorial[l-m]/Factorial[l+m]] (* 100 digit arbitrary precision! *)


(* ::Input::Initialization:: *)
Clear[\[CapitalPhi]mNumVec]
\[CapitalPhi]mNumVec[\[CapitalDelta]rVec_?ListQ,\[Rho]Vec_?ListQ][r0para_?NumericQ,Mpara_?NumericQ,qpara_:1][nbar_Integer,m_Integer] := Block[{q,F0,r0,i,mbar,lbar,\[Phi]lmnMatrixTilde,F0Block = Sqrt[1-2 Mpara/r0para],yVec = (\[Rho]Vec^2 - \[CapitalDelta]rVec^2/(1-2 Mpara/r0para))/r0para^2, zcBlock = -2 r0para Sqrt[(1-2 Mpara/r0para)/(1-3 Mpara/r0para)]},
\[Phi]lmnMatrixTilde =Table[2 \[Phi]lmn[lbar,mbar,nbar]/(1+KroneckerDelta[0,mbar]) /. {q -> qpara,r0 -> r0para,F0 -> F0Block},{lbar,0,3 nbar+3},{mbar,0,lbar}];

ParallelTable[Sum[If[\[Phi]lmnMatrixTilde[[lbar+1,mbar+1]] == 0,0,\[Phi]lmnMatrixTilde[[lbar+1,mbar+1]] clm[lbar,mbar]
Cos[mbar ArcTan[\[CapitalDelta]rVec[[i]]/F0Block,- r0para Sqrt[yVec[[i]]]]] my\[ScriptCapitalI]Num[\[Rho]Vec[[i]]][zcBlock][nbar,lbar,mbar,m]],
{lbar,0,3 nbar+3},{mbar,lbar,0,-1}]/(2 \[Pi])
,{i,1,Length[\[CapitalDelta]rVec]},Method->"CoarsestGrained"]]


(* ::Chapter::Closed:: *)
(*Seff_m in BH frame*)


(*Print["Seff Routines"]*)


(* ::Text:: *)
(*The effective source is just applying some operator \!\(-*)
(*\*SubscriptBox[\(\[Square]\), \(m\)]\) on Subscript[\[CapitalPhi], m]: \!\( *)
(*\*SubsuperscriptBox[\(S\), \(eff\), \(m\)]\  = \ \(\(-\( *)
(*\*SubscriptBox[\(\[Sum]\), \( *)
(*\*OverscriptBox[\(l\), \(_\)] *)
(*\*OverscriptBox[\(m\), \(_\)] *)
(*\*OverscriptBox[\(n\), \(_\)]\)]\ *)
(*\*SubscriptBox[*)
(*OverscriptBox[\(\[CapitalPhi]\), \(_\)], \( *)
(*\*OverscriptBox[\(l\), \(_\)] *)
(*\*OverscriptBox[\(m\), \(_\)] *)
(*\*OverscriptBox[\(n\), \(_\)]\)]\ *)
(*\*SubscriptBox[\(c\), \( *)
(*\*OverscriptBox[\(l\), \(_\)] *)
(*\*OverscriptBox[\(m\), \(_\)]\)]\ *)
(*\*SubscriptBox[\(\[Square]\), \(m\)]\((\(cos( *)
(*\*OverscriptBox[\(m\), \(_\)] *)
(*\*OverscriptBox[\(\[Phi]\), \(^\)])\)\ \(\[ScriptCapitalI]( *)
(*\*OverscriptBox[\(n\), \(_\)], *)
(*\*OverscriptBox[\(l\), \(_\)], *)
(*\*OverscriptBox[\(m\), \(_\)], m)\))\)\)\)\  = \ \(-\( *)
(*\*SubscriptBox[\(\[Sum]\), \( *)
(*\*OverscriptBox[\(l\), \(_\)] *)
(*\*OverscriptBox[\(m\), \(_\)] *)
(*\*OverscriptBox[\(n\), \(_\)]\)]\ *)
(*\*SubscriptBox[\(\[CapitalPhi]\), \( *)
(*\*OverscriptBox[\(l\), \(_\)] *)
(*\*OverscriptBox[\(m\), \(_\)] *)
(*\*OverscriptBox[\(n\), \(_\)]\)]\ *)
(*\*SubscriptBox[\(c\), \( *)
(*\*OverscriptBox[\(l\), \(_\)] *)
(*\*OverscriptBox[\(m\), \(_\)]\)]\ {*)
(*\*SubscriptBox[\(\[Square]\), \(m\)]\((cos( *)
(*\*OverscriptBox[\(m\), \(_\)] *)
(*\*OverscriptBox[\(\[Phi]\), \(^\)]))\)\ \[ScriptCapitalI]\  + \ \(cos( *)
(*\*OverscriptBox[\(m\), \(_\)] *)
(*\*OverscriptBox[\(\[Phi]\), \(^\)])\)\ *)
(*\*SubscriptBox[\(\[Square]\), \(m\)]\((\[ScriptCapitalI])\)\)\)\)\))}.*)


(* ::Input::Initialization:: *)
Clear[mModeBoxOperator0,mModeBoxOperator2,mModeBoxOperator2]
my\[CapitalOmega] = Sqrt[M/r0^3];
fTerm = 1-2 M/r;
mModeBoxOperator0[func_][m_] := m^2 (my\[CapitalOmega]^2/fTerm  -1/(r^2(1-y))) # & [func]
mModeBoxOperator1[func_] := (1+fTerm )/r D[#,\[CapitalDelta]r] + 2 (1-3 y) D[#,y]/r^2 & [func]
mModeBoxOperator2[func_] := fTerm  D[#,{\[CapitalDelta]r,2}]  +4 y (1-y) D[#,{y,2}]/r^2 & [func]


(* ::Section::Closed:: *)
(*\!\( *)
(*\*SubsuperscriptBox[\(\[Square]\), \(m\), \(0\)]\((cos( *)
(*\*OverscriptBox[\(m\), \(_\)] *)
(*\*OverscriptBox[\(\[Phi]\), \(^\)]))\)\)*)


(* ::Input::Initialization:: *)
Clear[box0Cosmbar\[Phi]bar]
box0Cosmbar\[Phi]bar[\[CapitalDelta]r_,y_][r0_,f0_][m_,mbar_] :=m^2 (1/((-1+y) (r0+\[CapitalDelta]r)^2)-((-1+f0) (r0+\[CapitalDelta]r))/(2 r0^2 (f0 r0+\[CapitalDelta]r))) Cos[mbar ArcTan[\[CapitalDelta]r/Sqrt[f0],-r0 Sqrt[y]]]


(* ::Section::Closed:: *)
(*\[Square]^1(cos(\!\(\*OverscriptBox[\(m\), \(_\)]\)\!\(\*OverscriptBox[\(\[Phi]\), \(^\)]\)))*)


(* ::Text:: *)
(*Here, you can change the value for ymin. Needed because d/dr cos(...)/y is formally 0/0 at y=0 and need to Taylor expand around y=0 to get correct value. For y<ymin, the leading order term*)
(*is used instead of evaluating d/dr cos(...)/y at y.*)


(* ::Input::Initialization:: *)
Clear[d\[CapitalDelta]rCosmbar\[Phi]bar,box1Cosmbar\[Phi]bar]
ymin = 10^(-15);
d\[CapitalDelta]rCosmbar\[Phi]bar[\[CapitalDelta]r_,y_][r0_,f0_][mbar_] := d\[CapitalDelta]rCosmbar\[Phi]bar[\[CapitalDelta]r,y][r0,f0][mbar] =   -((mbar r0 Sqrt[f0 y] Sin[mbar ArcTan[\[CapitalDelta]r,-r0 Sqrt[f0 y]]])/(f0 r0^2 y+\[CapitalDelta]r^2))
box1Cosmbar\[Phi]bar[\[CapitalDelta]r_,y_][r0_,f0_][mbar_] := 1/(r0+\[CapitalDelta]r)^2 ((r0+f0 r0+5 \[CapitalDelta]r) d\[CapitalDelta]rCosmbar\[Phi]bar[\[CapitalDelta]r,y][r0,f0][mbar]  - 
\[CapitalDelta]r If[y < ymin,(f0 mbar^2 r0^2)/\[CapitalDelta]r^3 If[\[CapitalDelta]r < 0,Cos[mbar \[Pi]],1],1/y d\[CapitalDelta]rCosmbar\[Phi]bar[\[CapitalDelta]r,y][r0,f0][mbar]])


(* ::Section::Closed:: *)
(*\[Square]^2(cos(\!\(\*OverscriptBox[\(m\), \(_\)]\)\!\(\*OverscriptBox[\(\[Phi]\), \(^\)]\)))*)


(* ::Input::Initialization:: *)
Clear[box2Cosmbar\[Phi]bar]
box2Cosmbar\[Phi]bar[\[CapitalDelta]r_,y_][r0_,f0_][mbar_] :=  If[y < ymin,0,-((mbar r0)/((r0+\[CapitalDelta]r)^2 (f0 r0^2 y+\[CapitalDelta]r^2)^2)) Sqrt[f0/y] (mbar r0 (r0 (f0 y)^(3/2) (r0+\[CapitalDelta]r)+\[CapitalDelta]r (r0 Sqrt[f0 y^3]+Sqrt[f0 y] \[CapitalDelta]r)) Cos[mbar ArcTan[\[CapitalDelta]r,-r0 Sqrt[f0 y]]]+\[CapitalDelta]r (f0 r0 y (r0-3 r0 y-2 \[CapitalDelta]r)+\[CapitalDelta]r (-2 r0 y+\[CapitalDelta]r-3 y \[CapitalDelta]r)) Sin[mbar ArcTan[\[CapitalDelta]r,-r0 Sqrt[f0 y]]])]


(* ::Section::Closed:: *)
(*\[Square]^1\[ScriptCapitalI](\!\(\*OverscriptBox[\(n\), \(_\)]\),\!\(\*OverscriptBox[\(l\), \(_\)]\),\!\(\*OverscriptBox[\(m\), \(_\)]\),m)*)


(* ::Input::Initialization:: *)
Clear[myd\[ScriptCapitalI]Num]
myd\[ScriptCapitalI]Num[\[Rho]_?NumericQ][zc_?NumericQ][nbar_?NumericQ,lbar_?NumericQ,lm1_?NumericQ,m_?NumericQ]  := 0 /; lm1 == lbar-1
myd\[ScriptCapitalI]Num[\[Rho]_?NumericQ][zc_?NumericQ][nbar_?NumericQ,lbar_?NumericQ,lbar_?NumericQ,m_?NumericQ] := myd\[ScriptCapitalI]Num[\[Rho]][zc][nbar,lbar,lbar,m]= 
(lbar my\[ScriptCapitalI]Num[\[Rho]][zc][nbar,lbar,lbar,m] + (nbar-lbar)/((2 lbar+1) (2 lbar+3)) my\[ScriptCapitalI]Num[\[Rho]][zc][nbar,lbar+2,lbar+2,m])/\[Rho]
myd\[ScriptCapitalI]Num[\[Rho]_?NumericQ][zc_?NumericQ][nbar_?NumericQ,lbar_?NumericQ,mbar_?NumericQ,m_?NumericQ]:= (myd\[ScriptCapitalI]Num[\[Rho]][zc][nbar,lbar,mbar,m] = (If[mbar+2 > lbar,0,-myd\[ScriptCapitalI]Num[\[Rho]][zc][nbar,lbar,mbar+2,m]] + If[mbar+2 > lbar-2,0,myd\[ScriptCapitalI]Num[\[Rho]][zc][nbar,lbar-2,mbar+2,m]] + 
If[mbar > lbar-2,0,(lbar+mbar) (lbar+mbar-1) myd\[ScriptCapitalI]Num[\[Rho]][zc][nbar,lbar-2,mbar,m]])/((lbar-mbar-1) (lbar-mbar))) /; mbar < lbar-1


(* ::Input::Initialization:: *)
Clear[myBox1\[ScriptCapitalI]pFactorNum]
myBox1\[ScriptCapitalI]pFactorNum[\[CapitalDelta]r_?NumericQ,\[Rho]_?NumericQ][r0_?NumericQ,f0_?NumericQ][n_?NumericQ,lbar_?NumericQ,mbar_?NumericQ,m_?NumericQ] := Block[{y = (\[Rho]^2-\[CapitalDelta]r^2/f0)/r0^2},((\[CapitalDelta]r (r0+f0 r0+2 \[CapitalDelta]r))/f0 + (1-3 y) r0^2)/(\[Rho] (r0+\[CapitalDelta]r)^2)]


(* ::Section::Closed:: *)
(*\[Square]^2\[ScriptCapitalI](\!\(\*OverscriptBox[\(n\), \(_\)]\),\!\(\*OverscriptBox[\(l\), \(_\)]\),\!\(\*OverscriptBox[\(m\), \(_\)]\),m)*)


(* ::Input::Initialization:: *)
Clear[myd2\[ScriptCapitalI]Num]
myd2\[ScriptCapitalI]Num[\[Rho]_?NumericQ][zc_?NumericQ][nbar_?NumericQ,1,0,m_?NumericQ] := 0;
myd2\[ScriptCapitalI]Num[\[Rho]_?NumericQ][zc_?NumericQ][nbar_?NumericQ,lbar_?NumericQ,lm1_?NumericQ,m_?NumericQ]  := 0 /; lm1 == lbar-1
myd2\[ScriptCapitalI]Num[\[Rho]_?NumericQ][zc_?NumericQ][nbar_?NumericQ,lbar_?NumericQ,lbar_?NumericQ,m_?NumericQ] :=myd2\[ScriptCapitalI]Num[\[Rho]][zc][nbar,lbar,lbar,m] =  
(lbar (lbar-1) my\[ScriptCapitalI]Num[\[Rho]][zc][nbar,lbar,lbar,m]+ (nbar-lbar)/(2 lbar+3) my\[ScriptCapitalI]Num[\[Rho]][zc][nbar,lbar+2,lbar+2,m] + ((nbar-lbar) (nbar-lbar-2))/((2 lbar+1) (2 lbar+3) (2 lbar+5) (2 lbar+7)) my\[ScriptCapitalI]Num[\[Rho]][zc][nbar,lbar+4,lbar+4,m])/\[Rho]^2
myd2\[ScriptCapitalI]Num[\[Rho]_?NumericQ][zc_?NumericQ][n_?NumericQ,lbar_?NumericQ,mbar_?NumericQ,m_?NumericQ] := (myd2\[ScriptCapitalI]Num[\[Rho]][zc][n,lbar,mbar,m] = (If[mbar+2 > lbar,0,-myd2\[ScriptCapitalI]Num[\[Rho]][zc][n,lbar,mbar+2,m]] + If[mbar+2 > lbar-2,0,myd2\[ScriptCapitalI]Num[\[Rho]][zc][n,lbar-2,mbar+2,m]] + 
If[mbar > lbar-2,0,(lbar+mbar) (lbar+mbar-1) myd2\[ScriptCapitalI]Num[\[Rho]][zc][n,lbar-2,mbar,m]])/((lbar-mbar-1) (lbar-mbar))) /; mbar < lbar-1


(* ::Input::Initialization:: *)
Clear[myBox2\[ScriptCapitalI]pFactorNum,myBox2\[ScriptCapitalI]ppFactorNum]
myBox2\[ScriptCapitalI]pFactorNum[\[CapitalDelta]r_?NumericQ,\[Rho]_?NumericQ][r0_?NumericQ,f0_?NumericQ][n_?NumericQ,lbar_?NumericQ,mbar_?NumericQ,m_?NumericQ] := Block[{y = (\[Rho]^2-\[CapitalDelta]r^2/f0)/r0^2},(r0^4 (-1+y) y+((r0+\[CapitalDelta]r) (f0 r0+\[CapitalDelta]r) (-\[CapitalDelta]r^2+f0 \[Rho]^2))/f0^2)/((r0+\[CapitalDelta]r)^2 \[Rho]^3) ]
myBox2\[ScriptCapitalI]ppFactorNum[\[CapitalDelta]r_?NumericQ,\[Rho]_?NumericQ][r0_?NumericQ,f0_?NumericQ][n_?NumericQ,lbar_?NumericQ,mbar_?NumericQ,m_?NumericQ] := Block[{y = (\[Rho]^2-\[CapitalDelta]r^2/f0)/r0^2},(-f0^2 r0^4 (-1+y) y+f0 r0 \[CapitalDelta]r^2 (r0+\[CapitalDelta]r)+\[CapitalDelta]r^3 (r0+\[CapitalDelta]r))/(f0^2 (r0+\[CapitalDelta]r)^2 \[Rho]^2) ]


(* ::Section::Closed:: *)
(*Cross-term*)


(* ::Input::Initialization:: *)
Clear[crossTerm\[ScriptCapitalI]pFactorNum]
crossTerm\[ScriptCapitalI]pFactorNum[\[CapitalDelta]r_?NumericQ,\[Rho]_?NumericQ][r0_?NumericQ,f0_?NumericQ][n_?NumericQ,lbar_?NumericQ,mbar_?NumericQ,m_?NumericQ] := 
Block[{y = (\[Rho]^2-\[CapitalDelta]r^2/f0)/r0^2},2 ((f0 r0+\[CapitalDelta]r)/(f0 r0+f0 \[CapitalDelta]r) - r0^2/(r0+\[CapitalDelta]r)^2 (1-y)) \[CapitalDelta]r d\[CapitalDelta]rCosmbar\[Phi]bar[\[CapitalDelta]r,y][r0,f0][mbar] /\[Rho]]


(* ::Section::Closed:: *)
(*The effective source *)


(* ::Text:: *)
(*Altogether, the effective source can be written in the form: \!\( *)
(*\*SubscriptBox[\(\[Sum]\), \( *)
(*\*OverscriptBox[\(l\), \(_\)] *)
(*\*OverscriptBox[\(m\), \(_\)] *)
(*\*OverscriptBox[\(n\), \(_\)]\)]\ \( *)
(*\*SubscriptBox[\(\[Phi]\), \( *)
(*\*OverscriptBox[\(l\), \(_\)] *)
(*\*OverscriptBox[\(m\), \(_\)] *)
(*\*OverscriptBox[\(n\), \(_\)]\)]\ \((\[ScriptCapitalI]Factor\ \[ScriptCapitalI]\  + \ \[ScriptCapitalI]pFactor\ \[ScriptCapitalI]'\  + \ \[ScriptCapitalI]ppFactor\ \[ScriptCapitalI]'')\)\)\), where \[ScriptCapitalI]' := (\[PartialD]/\[PartialD]\[Rho])\[ScriptCapitalI].*)


(* ::Input::Initialization:: *)
Clear[SeffmNumVec]
SeffmNumVec[\[CapitalDelta]rVec_?ListQ,\[Rho]Vec_?ListQ][r0para_?NumericQ,Mpara_?NumericQ,qpara_:1][nbar_Integer,m_Integer] := 
Block[{q,r0,F0,i,mbar,lbar,\[Phi]lmnMatrixTilde,f0Block = 1-2 Mpara/r0para, yVec = (\[Rho]Vec^2 - \[CapitalDelta]rVec^2/(1-2 Mpara/r0para))/r0para^2, zcBlock = -2 r0para Sqrt[(r0para-2 Mpara)/(r0para-3 Mpara)]},
\[Phi]lmnMatrixTilde =Table[2 \[Phi]lmn[lbar,mbar,nbar]/(1+KroneckerDelta[0,mbar]) /. {q -> qpara,r0 -> r0para,F0 -> Sqrt[f0Block]},{lbar,0,3 nbar+3},{mbar,0,lbar}];

ParallelTable[-Sum[If[\[Phi]lmnMatrixTilde[[lbar+1,mbar+1]] == 0,0,\[Phi]lmnMatrixTilde[[lbar+1,mbar+1]]  clm[lbar,mbar](
(box0Cosmbar\[Phi]bar[\[CapitalDelta]rVec[[i]],yVec[[i]]][r0para,f0Block][m,mbar]+ box1Cosmbar\[Phi]bar[\[CapitalDelta]rVec[[i]],yVec[[i]]][r0para,f0Block][mbar] + box2Cosmbar\[Phi]bar[\[CapitalDelta]rVec[[i]],yVec[[i]]][r0para,f0Block][mbar]) my\[ScriptCapitalI]Num[\[Rho]Vec[[i]]][zcBlock][nbar,lbar,mbar,m]+ (* Term proportional to \[ScriptCapitalI] *)
((myBox1\[ScriptCapitalI]pFactorNum[\[CapitalDelta]rVec[[i]],\[Rho]Vec[[i]]][r0para,f0Block] [nbar,lbar,mbar,m]+ myBox2\[ScriptCapitalI]pFactorNum[\[CapitalDelta]rVec[[i]],\[Rho]Vec[[i]]][r0para,f0Block] [nbar,lbar,mbar,m]) Cos[mbar ArcTan[\[CapitalDelta]rVec[[i]]/Sqrt[f0Block],- r0para Sqrt[yVec[[i]]]]]+crossTerm\[ScriptCapitalI]pFactorNum[\[CapitalDelta]rVec[[i]],\[Rho]Vec[[i]]][r0para,f0Block] [nbar,lbar,mbar,m]) myd\[ScriptCapitalI]Num[\[Rho]Vec[[i]]][zcBlock][nbar,lbar,mbar,m]+ (* Term proportional to \[ScriptCapitalI]' *)
Cos[mbar ArcTan[\[CapitalDelta]rVec[[i]]/Sqrt[f0Block],- r0para Sqrt[yVec[[i]]]]] myBox2\[ScriptCapitalI]ppFactorNum[\[CapitalDelta]rVec[[i]],\[Rho]Vec[[i]]][r0para,f0Block][nbar,lbar,mbar,m] myd2\[ScriptCapitalI]Num[\[Rho]Vec[[i]]][zcBlock][nbar,lbar,mbar,m] (* Term proportional to \[ScriptCapitalI]'' *) 
)]
,{lbar,0,3 nbar+3},{mbar,lbar,0,-1}]/(2 \[Pi])
,{i,1,Length[\[CapitalDelta]rVec]},Method->"CoarsestGrained"]]


(* ::Section:: *)
(*\!\( *)
(*\*SubscriptBox[\(\[PartialD]\), \(y\)]*)
(*\*SubscriptBox[*)
(*SubscriptBox[\(\[CapitalPhi]\), \(m\)], \(\n\)]\) and \!\( *)
(*\*SubscriptBox[\(\[PartialD]\), \(\[CapitalDelta]r\)]\ *)
(*\*SubscriptBox[\(\[CapitalPhi]\), \(m\)]\)*)


(* ::Text:: *)
(*These are needed because we will need to know the value of the normal derivative of \[CapitalPhi] on the shell/torus.*)


(* ::Input::Initialization:: *)
Clear[d\[CapitalPhi]mdyNumVec]
d\[CapitalPhi]mdyNumVec[\[CapitalDelta]rVec_?ListQ,\[Rho]Vec_?ListQ][r0para_?NumericQ,Mpara_?NumericQ,qpara_:1][nbar_Integer,m_Integer] := Block[{q,r0,F0,i,mbar,lbar,\[Phi]lmnMatrixTilde,f0Block = 1-2 Mpara/r0para,yVec = (\[Rho]Vec^2 - \[CapitalDelta]rVec^2/(1-2 Mpara/r0para))/r0para^2, zcBlock = -2 r0para Sqrt[(r0para-2 Mpara)/(r0para-3 Mpara)]},
\[Phi]lmnMatrixTilde =Table[2 \[Phi]lmn[lbar,mbar,nbar]/(1+KroneckerDelta[0,mbar]) /. {q -> qpara,r0 -> r0para,F0 -> Sqrt[f0Block]},{lbar,0,3 nbar+3},{mbar,0,lbar}];
ParallelTable[Sum[If[\[Phi]lmnMatrixTilde[[lbar+1,mbar+1]] == 0,0,\[Phi]lmnMatrixTilde[[lbar+1,mbar+1]] clm[lbar,mbar] (
-\[CapitalDelta]rVec[[i]]/2 If[yVec[[i]] < ymin,(f0Block mbar^2 r0para^2)/\[CapitalDelta]rVec[[i]]^3 If[\[CapitalDelta]rVec[[i]] < 0,Cos[mbar \[Pi]],1],1/yVec[[i]] d\[CapitalDelta]rCosmbar\[Phi]bar[\[CapitalDelta]rVec[[i]],yVec[[i]]][r0para,f0Block][mbar]] my\[ScriptCapitalI]Num[\[Rho]Vec[[i]]][zcBlock][nbar,lbar,mbar,m]+Cos[mbar ArcTan[\[CapitalDelta]rVec[[i]]/Sqrt[1-2 Mpara/r0para],- r0para Sqrt[yVec[[i]]]]] (r0para^2/(2 \[Rho]Vec[[i]]))myd\[ScriptCapitalI]Num[\[Rho]Vec[[i]]][zcBlock][nbar,lbar,mbar,m])
],{lbar,0,3 nbar+3},{mbar,lbar,0,-1}]/(2 \[Pi])
,{i,1,Length[\[CapitalDelta]rVec]},Method->"CoarsestGrained"]]


(* ::Input::Initialization:: *)
Clear[d\[CapitalPhi]md\[CapitalDelta]rNumVec]
d\[CapitalPhi]md\[CapitalDelta]rNumVec[\[CapitalDelta]rVec_?ListQ,\[Rho]Vec_?ListQ][r0para_?NumericQ,Mpara_?NumericQ,qpara_:1][nbar_Integer,m_Integer] := Block[{q,r0,F0,i,mbar,lbar,\[Phi]lmnMatrixTilde,f0Block = 1-2 Mpara/r0para,yVec = (\[Rho]Vec^2 - \[CapitalDelta]rVec^2/(1-2 Mpara/r0para))/r0para^2, zcBlock = -2 r0para Sqrt[(r0para-2 Mpara)/(r0para-3 Mpara)]},
\[Phi]lmnMatrixTilde =Table[2 \[Phi]lmn[lbar,mbar,nbar]/(1+KroneckerDelta[0,mbar]) /. {q -> qpara,r0 -> r0para,F0 -> Sqrt[f0Block]},{lbar,0,3 nbar+3},{mbar,0,lbar}];
ParallelTable[Sum[If[\[Phi]lmnMatrixTilde[[lbar+1,mbar+1]] == 0,0,\[Phi]lmnMatrixTilde[[lbar+1,mbar+1]] clm[lbar,mbar] (
d\[CapitalDelta]rCosmbar\[Phi]bar[\[CapitalDelta]rVec[[i]],yVec[[i]]][r0para,f0Block][mbar] my\[ScriptCapitalI]Num[\[Rho]Vec[[i]]][zcBlock][nbar,lbar,mbar,m]+Cos[mbar ArcTan[\[CapitalDelta]rVec[[i]]/Sqrt[1-2 Mpara/r0para],- r0para Sqrt[yVec[[i]]]]] (\[CapitalDelta]rVec[[i]]/(f0Block \[Rho]Vec[[i]]))myd\[ScriptCapitalI]Num[\[Rho]Vec[[i]]][zcBlock][nbar,lbar,mbar,m])
],{lbar,0,3 nbar+3},{mbar,lbar,0,-1}]/(2 \[Pi])
,{i,1,Length[\[CapitalDelta]rVec]},Method->"CoarsestGrained"]]


(* ::Section::Closed:: *)
(*\!\( *)
(*\*SubscriptBox[\(\[PartialD]\), \(yy\)]*)
(*\*SubscriptBox[*)
(*SubscriptBox[\(\[CapitalPhi]\), \(m\)], \(\n\)]\) and \!\( *)
(*\*SubscriptBox[\(\[PartialD]\), \(\[CapitalDelta]r\[CapitalDelta]r\)]\ *)
(*\*SubscriptBox[\(\[CapitalPhi]\), \(m\)]\)*)


(* ::Text:: *)
(*These are temporary to check that everything works fine.*)


(* ::Input::Initialization:: *)
Clear[d\[CapitalPhi]mdyyNumVec]
d\[CapitalPhi]mdyyNumVec[\[CapitalDelta]rVec_?ListQ,\[Rho]Vec_?ListQ][r0para_?NumericQ,Mpara_?NumericQ,qpara_:1][nbar_Integer,m_Integer] := Block[{q,r0,F0,i,mbar,lbar,\[Phi]lmnMatrixTilde,f0Block = 1-2 Mpara/r0para,yVec = (\[Rho]Vec^2 - \[CapitalDelta]rVec^2/(1-2 Mpara/r0para))/r0para^2, zcBlock = -2 r0para Sqrt[(r0para-2 Mpara)/(r0para-3 Mpara)]},
\[Phi]lmnMatrixTilde =Table[2 \[Phi]lmn[lbar,mbar,nbar]/(1+KroneckerDelta[0,mbar]) /. {q -> qpara,r0 -> r0para,F0 -> Sqrt[f0Block]},{lbar,0,3 nbar+3},{mbar,0,lbar}];
ParallelTable[Sum[If[\[Phi]lmnMatrixTilde[[lbar+1,mbar+1]] == 0,0,\[Phi]lmnMatrixTilde[[lbar+1,mbar+1]] clm[lbar,mbar] (
(If[yVec[[i]] < ymin,(f0Block^2 mbar^2 (8+mbar^2) r0para^4)/(12 \[CapitalDelta]rVec[[i]]^4) If[\[CapitalDelta]rVec[[i]] <0,(-1)^mbar,1],-(1/(4 yVec[[i]]^(3/2) (f0Block r0para^2 yVec[[i]]+\[CapitalDelta]rVec[[i]]^2)^2)) Sqrt[f0Block] mbar r0para \[CapitalDelta]rVec[[i]] (mbar r0para Sqrt[f0Block yVec[[i]]] \[CapitalDelta]rVec[[i]] Cos[mbar ArcTan[\[CapitalDelta]rVec[[i]],-r0para Sqrt[f0Block yVec[[i]]]]]+(3 f0Block r0para^2 yVec[[i]]+\[CapitalDelta]rVec[[i]]^2) Sin[mbar ArcTan[\[CapitalDelta]rVec[[i]],-r0para Sqrt[f0Block yVec[[i]]]]])]) my\[ScriptCapitalI]Num[\[Rho]Vec[[i]]][zcBlock][nbar,lbar,mbar,m]+
(Cos[mbar ArcTan[\[CapitalDelta]rVec[[i]]/Sqrt[1-2 Mpara/r0para],- r0para Sqrt[yVec[[i]]]]] ((r0para^4 (-1))/(4 \[Rho]Vec[[i]]^3))- \[CapitalDelta]rVec[[i]] r0para^2/(2 \[Rho]Vec[[i]]) If[yVec[[i]] < ymin,(f0Block mbar^2 r0para^2)/\[CapitalDelta]rVec[[i]]^3 If[\[CapitalDelta]rVec[[i]] < 0,Cos[mbar \[Pi]],1],1/yVec[[i]] d\[CapitalDelta]rCosmbar\[Phi]bar[\[CapitalDelta]rVec[[i]],yVec[[i]]][r0para,f0Block][mbar]]) myd\[ScriptCapitalI]Num[\[Rho]Vec[[i]]][zcBlock][nbar,lbar,mbar,m]+
Cos[mbar ArcTan[\[CapitalDelta]rVec[[i]]/Sqrt[1-2 Mpara/r0para],- r0para Sqrt[yVec[[i]]]]](r0para^4/(4 \[Rho]Vec[[i]]^2)) myd2\[ScriptCapitalI]Num[\[Rho]Vec[[i]]][zcBlock][nbar,lbar,mbar,m])
],{lbar,0,3 nbar+3},{mbar,lbar,0,-1}]/(2 \[Pi])
,{i,1,Length[\[CapitalDelta]rVec]},Method->"CoarsestGrained"]]


(* ::Input::Initialization:: *)
Clear[d\[CapitalPhi]md\[CapitalDelta]r\[CapitalDelta]rNumVec]
d\[CapitalPhi]md\[CapitalDelta]r\[CapitalDelta]rNumVec[\[CapitalDelta]rVec_?ListQ,\[Rho]Vec_?ListQ][r0para_?NumericQ,Mpara_?NumericQ,qpara_:1][nbar_Integer,m_Integer] := Block[{q,r0,F0,i,mbar,lbar,\[Phi]lmnMatrixTilde,f0Block = 1-2 Mpara/r0para,yVec = (\[Rho]Vec^2 - \[CapitalDelta]rVec^2/(1-2 Mpara/r0para))/r0para^2, zcBlock = -2 r0para Sqrt[(r0para-2 Mpara)/(r0para-3 Mpara)]},
\[Phi]lmnMatrixTilde =Table[2 \[Phi]lmn[lbar,mbar,nbar]/(1+KroneckerDelta[0,mbar]) /. {q -> qpara,r0 -> r0para,F0 -> Sqrt[f0Block]},{lbar,0,3 nbar+3},{mbar,0,lbar}];
ParallelTable[Sum[If[\[Phi]lmnMatrixTilde[[lbar+1,mbar+1]] == 0,0,\[Phi]lmnMatrixTilde[[lbar+1,mbar+1]] clm[lbar,mbar] (
(1/(f0Block r0para^2 yVec[[i]]+\[CapitalDelta]rVec[[i]]^2)^2 mbar r0para (-f0Block mbar r0para yVec[[i]] Cos[mbar ArcTan[\[CapitalDelta]rVec[[i]],-r0para Sqrt[f0Block yVec[[i]]]]]+2 Sqrt[f0Block yVec[[i]]] \[CapitalDelta]rVec[[i]] Sin[mbar ArcTan[\[CapitalDelta]rVec[[i]],-r0para Sqrt[f0Block yVec[[i]]]]])) my\[ScriptCapitalI]Num[\[Rho]Vec[[i]]][zcBlock][nbar,lbar,mbar,m]+
(Cos[mbar ArcTan[\[CapitalDelta]rVec[[i]]/Sqrt[1-2 Mpara/r0para],- r0para Sqrt[yVec[[i]]]]] ((-\[CapitalDelta]rVec[[i]]^2+f0Block \[Rho]Vec[[i]]^2)/(f0Block^2 \[Rho]Vec[[i]]^3))+2 \[CapitalDelta]rVec[[i]] /(f0Block \[Rho]Vec[[i]]) d\[CapitalDelta]rCosmbar\[Phi]bar[\[CapitalDelta]rVec[[i]],yVec[[i]]][r0para,f0Block][mbar]) myd\[ScriptCapitalI]Num[\[Rho]Vec[[i]]][zcBlock][nbar,lbar,mbar,m]+
Cos[mbar ArcTan[\[CapitalDelta]rVec[[i]]/Sqrt[1-2 Mpara/r0para],- r0para Sqrt[yVec[[i]]]]](\[CapitalDelta]rVec[[i]]^2/(f0Block^2 \[Rho]Vec[[i]]^2)) myd2\[ScriptCapitalI]Num[\[Rho]Vec[[i]]][zcBlock][nbar,lbar,mbar,m])
],{lbar,0,3 nbar+3},{mbar,lbar,0,-1}]/(2 \[Pi])
,{i,1,Length[\[CapitalDelta]rVec]},Method->"CoarsestGrained"]]


(* ::Chapter::Closed:: *)
(*(\[Sigma],y)-Domain Decomposition*)


CoordFlagValue


(* ::Section:: *)
(*Domain 0*)


(* ::Subsection:: *)
(*Coordinate \[Sigma](x1,x2)*)


\[Sigma]$Dom0=\[Sigma]m/2(1+x1);

\[Sigma]1$Dom0=D[\[Sigma]$Dom0,x1];
\[Sigma]11$Dom0=D[\[Sigma]$Dom0,{x1,2}];
\[Sigma]12$Dom0=D[\[Sigma]$Dom0,x1,x2];
\[Sigma]2$Dom0=D[\[Sigma]$Dom0,x2];
\[Sigma]22$Dom0=D[\[Sigma]$Dom0,{x2,2}];


(* ::Subsection:: *)
(*Coordinate y(x1,x2)*)


y$Dom0=(1+x2)/2;
y1$Dom0=D[y$Dom0,x1];
y11$Dom0=D[y$Dom0,{x1,2}];
y12$Dom0=D[y$Dom0,x1,x2];
y2$Dom0=D[y$Dom0,x2];
y22$Dom0=D[y$Dom0,{x2,2}];


(* ::Section:: *)
(*Domain 1*)


(* ::Subsection:: *)
(*Coordinate \[Sigma](x1,x2)*)


\[Sigma]$Dom1=\[Sigma]0/(1-\[Eta] \[Sigma]0 Sqrt[f0] x1);

\[Sigma]1$Dom1=D[\[Sigma]$Dom1,x1];
\[Sigma]11$Dom1=D[\[Sigma]$Dom1,{x1,2}];
\[Sigma]12$Dom1=D[\[Sigma]$Dom1,x1,x2];
\[Sigma]2$Dom1=D[\[Sigma]$Dom1,x2];
\[Sigma]22$Dom1=D[\[Sigma]$Dom1,{x2,2}];


(* ::Subsection::Closed:: *)
(*Coordinate y(x1,x2)*)


y$Dom1 = \[Eta]^2\[Sigma]0^2(1-x1^2)(1-x2)/2+(1+x2)/2;

y1=y1$Dom1=D[y$Dom1,x1];
y11=y11$Dom1=D[y$Dom1,{x1,2}];
y12=y12$Dom1=D[y$Dom1,x1,x2];
y2=y2$Dom1=D[y$Dom1,x2];
y22=y22$Dom1=D[y$Dom1,{x2,2}];


(* ::Section:: *)
(*Domain 2*)


(* ::Subsection:: *)
(*Coordinate \[Sigma](x1,x2)*)


\[Rho]=(1+x2)/2;
\[Sigma]$Dom2=\[Sigma]0/(1-\[Eta] \[Sigma]0 Sqrt[f0] \[Rho] x1);

\[Sigma]1$Dom2=D[\[Sigma]$Dom2,x1];
\[Sigma]11$Dom2=D[\[Sigma]$Dom2,{x1,2}];
\[Sigma]12$Dom2=D[\[Sigma]$Dom2,x1,x2];
\[Sigma]2$Dom2=D[\[Sigma]$Dom2,x2];
\[Sigma]22$Dom2=D[\[Sigma]$Dom2,{x2,2}];


(* ::Subsection:: *)
(*Coordinate y(x1,x2)*)


\[Rho]=(1+x2)/2;

y$Dom2=If[CoordFlagValue==1,\[Eta]^2 \[Sigma]0^2 \[Rho] (1-x1^2), \[Eta]^2 \[Sigma]0^2 \[Rho]^2(1-x1^2)];
(*y$Dom2=\[Eta]^2 \[Sigma]0^2 \[Rho]^2(1-x1^2);*)

y1$Dom2=D[y$Dom2,x1];
y11$Dom2=D[y$Dom2,{x1,2}];
y12$Dom2=D[y$Dom2,x1,x2];
y2$Dom2=D[y$Dom2,x2];
y22$Dom2=D[y$Dom2,{x2,2}];


(* ::Section::Closed:: *)
(*Domain 3*)


(* ::Subsection:: *)
(*Coordinate \[Sigma](x1,x2)*)


\[Sigma]$Dom3=(1-x1)\[Sigma]p/2+(1+x1)/2;

\[Sigma]1$Dom3=D[\[Sigma]$Dom3,x1];
\[Sigma]11$Dom3=D[\[Sigma]$Dom3,{x1,2}];
\[Sigma]12$Dom3=D[\[Sigma]$Dom3,x1,x2];
\[Sigma]2$Dom3=D[\[Sigma]$Dom3,x2];
\[Sigma]22$Dom3=D[\[Sigma]$Dom3,{x2,2}];


(* ::Subsection::Closed:: *)
(*Coordinate y(x1,x2)*)


y$Dom3=(1+x2)/2;
y1$Dom3=D[y$Dom3,x1];
y11$Dom3=D[y$Dom3,{x1,2}];
y12$Dom3=D[y$Dom3,x1,x2];
y2$Dom3=D[y$Dom3,x2];
y22$Dom3=D[y$Dom3,{x2,2}];


(* ::Chapter::Closed:: *)
(*Hyperboloidal Puncture and Effective Source*)


(*Print["Set Hyp"]*)


(* ::Section:: *)
(*Coordinate Transformation*)


r$of$\[Sigma][\[Sigma]_]:=rh/\[Sigma];
\[Sigma]$of$r[r_]:=rh/r;

\[CapitalDelta]r[\[Sigma]_]:=r$of$\[Sigma][\[Sigma]]-r0;

\[Sigma]0=\[Sigma]$of$r[r0];


(* ::Subsection:: *)
(*Hyperboloidal Transformation*)


(*\[Lambda]=2rh;*)
s0=-I*\[Omega]0*\[Lambda]
r$toroise[r_]:= r + rh Log[r/rh-1]
x[\[Sigma]_]:=r$toroise[r$of$\[Sigma][\[Sigma]]]/\[Lambda]
H[\[Sigma]_]:=rh/\[Lambda] (Log[1-\[Sigma]]-1/\[Sigma]+Log[\[Sigma]]);

\[CapitalOmega][\[Sigma]_]:=\[Sigma]/\[Lambda]


(* ::Section:: *)
(*Re-scalling factor for Field*)


s=s0;
Z[\[Sigma]_,y_]:=(1-y)^(-m/2) * \[CapitalOmega][\[Sigma]] * Exp[s*H[\[Sigma]]]

dZdy[\[Sigma]_,y_]:= m/2 (1-y)^(-m/2-1) * \[CapitalOmega][\[Sigma]] * Exp[s*H[\[Sigma]]]

d2Zdy2[\[Sigma]_,y_]:=\!\(\*SuperscriptBox[\(dZdy\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Sigma],y]
dZd\[Sigma][\[Sigma]_,y_]:= (1-y)^(-m/2) * Exp[s*H[\[Sigma]]] ( 1/\[Lambda] + \[CapitalOmega][\[Sigma]] * s * (rh (-1+2 \[Sigma]^2))/(\[Lambda] (-1+\[Sigma]) \[Sigma]^2))
d2Zd\[Sigma]2[\[Sigma]_,y_]:=\!\(\*SuperscriptBox[\(dZd\[Sigma]\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[\[Sigma],y]


(* ::Section:: *)
(*Re-scalling factor for Effective Source*)


F[\[Sigma]_,y_]:=Z[\[Sigma],y]/r$of$\[Sigma][\[Sigma]]^2


(* ::Section:: *)
(*Hyperboloidal Puncture*)


\[CapitalPhi]mNumVec$hyp[\[Sigma]Vec_,\[Rho]Vec_][nbar_]:= With[{yVec = (\[Rho]Vec^2-\[CapitalDelta]r[\[Sigma]Vec]^2/f0)/r0^2},								
								\[CapitalPhi]mNumVec[\[CapitalDelta]r[\[Sigma]Vec],\[Rho]Vec][r0,M][nbar,m] / Z[\[Sigma]Vec,yVec]
								];

d\[CapitalPhi]mdyNumVec$hyp[\[Sigma]Vec_,\[Rho]Vec_][nbar_]:= With[{yVec = (\[Rho]Vec^2-\[CapitalDelta]r[\[Sigma]Vec]^2/f0)/r0^2}, 
								d\[CapitalPhi]mdyNumVec[\[CapitalDelta]r[\[Sigma]Vec],\[Rho]Vec][r0,M][nbar,m]/Z[\[Sigma]Vec,yVec] + \[CapitalPhi]mNumVec[\[CapitalDelta]r[\[Sigma]Vec],\[Rho]Vec][r0,M][nbar,m] * (-dZdy[\[Sigma]Vec,yVec]/Z[\[Sigma]Vec,yVec]^2)
								];

d\[CapitalPhi]md\[Sigma]NumVec$hyp[\[Sigma]Vec_,\[Rho]Vec_][nbar_]:= With[{yVec = (\[Rho]Vec^2-\[CapitalDelta]r[\[Sigma]Vec]^2/f0)/r0^2},
								r$of$\[Sigma]'[\[Sigma]Vec]*d\[CapitalPhi]md\[CapitalDelta]rNumVec[\[CapitalDelta]r[\[Sigma]Vec],\[Rho]Vec][r0,M][nbar,m]/Z[\[Sigma]Vec,yVec] + \[CapitalPhi]mNumVec[\[CapitalDelta]r[\[Sigma]Vec],\[Rho]Vec][r0,M][nbar,m]*(-dZd\[Sigma][\[Sigma]Vec,yVec]/Z[\[Sigma]Vec,yVec]^2)
								];
								
								
d2\[CapitalPhi]md\[Sigma]2NumVec$hyp[\[Sigma]Vec_,\[Rho]Vec_][nbar_]:= With[{yVec = (\[Rho]Vec^2-\[CapitalDelta]r[\[Sigma]Vec]^2/f0)/r0^2},

r$of$\[Sigma]''[\[Sigma]Vec]  *    d\[CapitalPhi]md\[CapitalDelta]rNumVec[\[CapitalDelta]r[\[Sigma]Vec],\[Rho]Vec][r0,M][nbar,m]   /    Z[\[Sigma]Vec,yVec] + (*((r^\[Prime]\[Prime])[\[Sigma]] (\[CapitalPhi]^(1,0))[r[\[Sigma]],y])/ZZ[\[Sigma],y]*)

r$of$\[Sigma]'[\[Sigma]Vec]  * r$of$\[Sigma]'[\[Sigma]Vec]  *   d\[CapitalPhi]md\[CapitalDelta]r\[CapitalDelta]rNumVec[\[CapitalDelta]r[\[Sigma]Vec],\[Rho]Vec][r0,M][nbar,m]   /    Z[\[Sigma]Vec,yVec] (*((r^\[Prime])[\[Sigma]]^2 (\[CapitalPhi]^(2,0))[r[\[Sigma]],y])/ZZ[\[Sigma],y]*)

+ 2*d\[CapitalPhi]md\[CapitalDelta]rNumVec[\[CapitalDelta]r[\[Sigma]Vec],\[Rho]Vec][r0,M][nbar,m]*(-dZd\[Sigma][\[Sigma]Vec,yVec]/Z[\[Sigma]Vec,yVec]^2) (*-((2 (r^\[Prime])[\[Sigma]] (ZZ^(1,0))[\[Sigma],y] (\[CapitalPhi]^(1,0))[r[\[Sigma]],y])/ZZ[\[Sigma],y]^2)*)

+ \[CapitalPhi]mNumVec[\[CapitalDelta]r[\[Sigma]Vec],\[Rho]Vec][r0,M][nbar,m]*(-d2Zd\[Sigma]2[\[Sigma]Vec,yVec]/Z[\[Sigma]Vec,yVec]^2) (*-((\[CapitalPhi][r[\[Sigma]],y] (ZZ^(2,0))[\[Sigma],y])/ZZ[\[Sigma],y]^2)*)

+ \[CapitalPhi]mNumVec[\[CapitalDelta]r[\[Sigma]Vec],\[Rho]Vec][r0,M][nbar,m]*(2* dZd\[Sigma][\[Sigma]Vec,yVec]*dZd\[Sigma][\[Sigma]Vec,yVec]/Z[\[Sigma]Vec,yVec]^3) (*(2 \[CapitalPhi][r[\[Sigma]],y] (ZZ^(1,0))[\[Sigma],y]^2)/ZZ[\[Sigma],y]^3*)
									];

d2\[CapitalPhi]mdy2NumVec$hyp[\[Sigma]Vec_,\[Rho]Vec_][nbar_]:= With[{yVec = (\[Rho]Vec^2-\[CapitalDelta]r[\[Sigma]Vec]^2/f0)/r0^2}, 
d\[CapitalPhi]mdyyNumVec[\[CapitalDelta]r[\[Sigma]Vec],\[Rho]Vec][r0,M][nbar,m]/Z[\[Sigma]Vec,yVec] (*(\[CapitalPhi]^(0,2))[r[\[Sigma]],y]/ZZ[\[Sigma],y]*)
 + 2*d\[CapitalPhi]mdyNumVec[\[CapitalDelta]r[\[Sigma]Vec],\[Rho]Vec][r0,M][nbar,m] * (-dZdy[\[Sigma]Vec,yVec]/Z[\[Sigma]Vec,yVec]^2) (*-((2 (ZZ^(0,1))[\[Sigma],y] (\[CapitalPhi]^(0,1))[r[\[Sigma]],y])/ZZ[\[Sigma],y]^2)*)
 + \[CapitalPhi]mNumVec[\[CapitalDelta]r[\[Sigma]Vec],\[Rho]Vec][r0,M][nbar,m] * (-d2Zdy2[\[Sigma]Vec,yVec]/Z[\[Sigma]Vec,yVec]^2 + 2*(dZdy[\[Sigma]Vec,yVec])^2/Z[\[Sigma]Vec,yVec]^3) (*\[CapitalPhi][r[\[Sigma]],y] ((2 (ZZ^(0,1))[\[Sigma],y]^2)/ZZ[\[Sigma],y]^3-(ZZ^(0,2))[\[Sigma],y]/ZZ[\[Sigma],y]^2)*)
								];


(* ::Section:: *)
(*Hyperboloidal Effective Source*)


SeffmNum$hyp[\[Sigma]Vec_,\[Rho]Vec_][nbar_]:=With[{yVec = (\[Rho]Vec^2-\[CapitalDelta]r[\[Sigma]Vec]^2/f0)/r0^2},								
									SeffmNumVec[\[CapitalDelta]r[\[Sigma]Vec],\[Rho]Vec][r0,M][nbar,m]/F[\[Sigma]Vec,yVec]
								];


(* ::Chapter:: *)
(*Calculate Puncture and Effective Source*)


(* ::Section::Closed:: *)
(*Derived Parameters*)


(*Print["Set Par"]*)


(*ParSpaceDirName="ParameterSpace"
InputParFile=ToString[$CommandLine[[Length@$CommandLine]]];
Get[InputParFile];*)

(*q=1;
M=1;

r0=10 M;*)
f0 = 1-rh/r0;
rh = 2*M;




(*m=0;
nbar=0;
nbarmax=nbar;
NN=100;
Prec=4 MachinePrecision;*)

\[Eta]=N[(r0/rh)^2*Sqrt[1-rh/r0]/(1+r0/rh),Prec];

Print["\t\tM=",M," r0=",r0," eta=",\[Eta]," m=",m, " nbar=",nbarmax," N=",NN]


\[Delta]r=N[Sqrt[f0]*\[Eta]*rh,Prec]

rplus=r0 + \[Delta]r;
\[Sigma]minus=\[Sigma]$of$r[rplus]//FullSimplify

rminus=r0 - \[Delta]r;
\[Sigma]plus=\[Sigma]$of$r[rminus]//FullSimplify


(* ::Subsubsection:: *)
(*Orbital Parameters*)


E0 = f0/Sqrt[1 - 3M/r0];

L0 = Sqrt[r0 M]/Sqrt[1-3M/r0];

\[CapitalOmega]0=Sqrt[M/r0^3];
\[Omega]0=m*\[CapitalOmega]0;


(* ::Section:: *)
(*Numerical Data*)


(* ::Subsection:: *)
(*Numerical Grid*)


(*Print["Set Num Data"]*)


(*Prec=2*MachinePrecision;*)
(*NN=100*)
N1=NN; 

n1=N1+1;
N2=NN; 

n2=N2+1;
x1$grid=Table[Grid$Lobatto[i1,N1], {i1,0,N1}];
x2$gridLobatto=Table[Grid$Lobatto[i2,N2], {i2,0,N2}];
x2$grid=Table[Grid$RadauRHS[i2,N2], {i2,0,N2}]; 


(* ::Subsubsection:: *)
(*1D Grids Domain 2*)


(*Print["Set 1D Dom 2"]*)


\[Sigma]$Surface$dom2=N[\[Sigma]$Dom2/.{x1->x1$grid,x2->1},Prec];
y$Surface$dom2=N[y$Dom2/.{x1->x1$grid,x2->1},Prec];
\[Rho]$Surface$dom2=N[Sqrt[\[CapitalDelta]r[\[Sigma]$Surface$dom2]^2/(1-2 M/r0) + r0^2 y$Surface$dom2],Prec];
\[Sigma]y$Surface$dom2=Transpose@{\[Sigma]$Surface$dom2,y$Surface$dom2};

\[Sigma]$Bulk$x1$dom2=N[\[Sigma]$Dom2/.{x1->x1$grid,x2->0},Prec];
y$Bulk$x1$dom2=N[y$Dom2/.{x1->x1$grid,x2->0},Prec];
\[Rho]$Bulk$x1$dom2=N[Sqrt[\[CapitalDelta]r[\[Sigma]$Bulk$x1$dom2]^2/(1-2 M/r0) + r0^2 y$Bulk$x1$dom2],Prec];
\[Sigma]y$Bulk$x1$dom2=Transpose@{\[Sigma]$Bulk$x1$dom2,y$Bulk$x1$dom2};

\[Sigma]$Particle=N[\[Sigma]$Dom2/.{x1->x1$grid,x2->x2$grid[[n2]]},Prec];
y$Particle=N[y$Dom2/.{x1->x1$grid,x2->x2$grid[[n2]]},Prec];
\[Rho]$Particle=N[Sqrt[\[CapitalDelta]r[\[Sigma]$Particle]^2/(1-2 M/r0) + r0^2 y$Particle],Prec];
\[Sigma]y$Particle=Transpose@{\[Sigma]$Particle,y$Particle};

\[Sigma]$BoundRight=N[\[Sigma]$Dom2/.{x1->1,x2->x2$grid},Prec];
y$BoundRight=N[y$Dom2/.{x1->1,x2->x2$grid},Prec];
\[Rho]$BoundRight=N[Sqrt[\[CapitalDelta]r[\[Sigma]$BoundRight]^2/(1-2 M/r0) + r0^2 y$BoundRight],Prec];
\[Sigma]y$BoundRight=Transpose@{\[Sigma]$BoundRight,y$BoundRight};


\[Sigma]$Bulk$x2=N[\[Sigma]$Dom2/.{x1->0,x2->x2$grid},Prec];
y$Bulk$x2=N[y$Dom2/.{x1->0,x2->x2$grid},Prec];
\[Rho]$Bulk$x2=N[Sqrt[\[CapitalDelta]r[\[Sigma]$Bulk$x2]^2/(1-2 M/r0) + r0^2 y$Bulk$x2],Prec];
\[Sigma]y$Bulk$x2=Transpose@{\[Sigma]$Bulk$x2,y$Bulk$x2};

\[Sigma]$BoundLeft=N[\[Sigma]$Dom2/.{x1->-1,x2->x2$grid},Prec];
y$BoundLeft=N[y$Dom2/.{x1->-1,x2->x2$grid},Prec];
\[Rho]$BoundLeft=N[Sqrt[\[CapitalDelta]r[\[Sigma]$BoundLeft]^2/(1-2 M/r0) + r0^2 y$BoundLeft],Prec];
\[Sigma]y$BoundLeft=Transpose@{\[Sigma]$BoundLeft,y$BoundLeft};


(* ::Subsubsection:: *)
(*2D Bulk Grid Domain 2*)


(*Print["Set 2D Dom 2"]*)


\[Sigma]$Bulk$dom2=Table[N[\[Sigma]$Dom2/.{x1->x1$grid[[i1]],x2->x2$grid[[i2]]},Prec],{i1,1,n1},{i2,1,n2}];
y$Bulk$dom2=Table[N[y$Dom2/.{x1->x1$grid[[i1]],x2->x2$grid[[i2]]},Prec],{i1,1,n1},{i2,1,n2}];
\[Rho]$Bulk$dom2=N[Sqrt[\[CapitalDelta]r[\[Sigma]$Bulk$dom2]^2/(1-2 M/r0) + r0^2 y$Bulk$dom2],Prec];


(* ::Section:: *)
(*Surface Boundaries in Domain 2*)


(* ::Subsection:: *)
(*Puncture Field at Boundary*)


(*nKernels=1;*)
(*LaunchKernels[nKernels];*)
(*Print[$KernelCount," ", $ProcessorCount];*)


(*Print["Calculate Punc"]*)


\[CapitalPhi]m$Punc$Surface$BENCHMARK=AbsoluteTiming[Sum[\[CapitalPhi]mNumVec$hyp[\[Sigma]$Surface$dom2,\[Rho]$Surface$dom2][nbar],{nbar,-1,nbarmax}]];
d\[CapitalPhi]m$Punc$dy$Surface$BENCHMARK=AbsoluteTiming[Sum[d\[CapitalPhi]mdyNumVec$hyp[\[Sigma]$Surface$dom2,\[Rho]$Surface$dom2][nbar],{nbar,-1,nbarmax}]];
d\[CapitalPhi]m$Punc$d\[Sigma]$Surface$BENCHMARK=AbsoluteTiming[Sum[d\[CapitalPhi]md\[Sigma]NumVec$hyp[\[Sigma]$Surface$dom2,\[Rho]$Surface$dom2][nbar],{nbar,-1,nbarmax}]];

\[CapitalPhi]m$Punc$Surface$Time=\[CapitalPhi]m$Punc$Surface$BENCHMARK[[1]];
\[CapitalPhi]m$Punc$Surface=\[CapitalPhi]m$Punc$Surface$BENCHMARK[[2]];


d\[CapitalPhi]m$Punc$dy$Surface$Time=d\[CapitalPhi]m$Punc$dy$Surface$BENCHMARK[[1]];
d\[CapitalPhi]m$Punc$dy$Surface=d\[CapitalPhi]m$Punc$dy$Surface$BENCHMARK[[2]];

d\[CapitalPhi]m$Punc$d\[Sigma]$Surface$Time=d\[CapitalPhi]m$Punc$d\[Sigma]$Surface$BENCHMARK[[1]];
d\[CapitalPhi]m$Punc$d\[Sigma]$Surface=d\[CapitalPhi]m$Punc$d\[Sigma]$Surface$BENCHMARK[[2]];
(*CloseKernels[];*)


(* ::Subsubsection:: *)
(*Chebyshev Surface*)


(*Print["Calculate Cheb Punc"]*)


\[CapitalPhi]m$Punc$Surface$Cheb=ChebLobattoCoefficients[\[CapitalPhi]m$Punc$Surface, Prec];
d\[CapitalPhi]m$Puncdy$Surface$Cheb=ChebLobattoCoefficients[d\[CapitalPhi]m$Punc$dy$Surface, Prec];
d\[CapitalPhi]m$Puncd\[Sigma]$Surface$Cheb=ChebLobattoCoefficients[d\[CapitalPhi]m$Punc$d\[Sigma]$Surface, Prec];


(* ::Section:: *)
(*Bulk Data Dom 2*)


(* ::Subsection:: *)
(*Effective Source in whole Domain*)


(*nKernels=8;*)
(*LaunchKernels[nKernels];*)
(*Print["Get Seff"]*)

Seffm$Flat$BENCHMARK=AbsoluteTiming[Sum[SeffmNum$hyp[Flatten[\[Sigma]$Bulk$dom2],Flatten[\[Rho]$Bulk$dom2]][nbar],{nbar,-1,nbarmax}]];

Seffm$Flat$Time = Seffm$Flat$BENCHMARK[[1]];
Seffm$Flat=Seffm$Flat$BENCHMARK[[2]];

Seffm$BENCHMARK=AbsoluteTiming[Partition[Seffm$Flat,n2]];

Seffm$Time=Seffm$Flat$Time+Seffm$BENCHMARK[[1]];
Seffm=Seffm$BENCHMARK[[2]];
(*CloseKernels[];*)


(* ::Subsubsection:: *)
(*Chebyschev Coefficients*)


(*Print["Get Cheb Seff"]*)
Seffm$ChebCoeff1=Table[ChebLobattoCoefficients[Seffm[[All,i2]], Prec],{i2,1,n2}];
Seffm$ChebCoeff=Table[ChebRadauRHSCoefficients[Seffm$ChebCoeff1[[All,i1]], Prec],{i1,1,n1}];


(* ::Section:: *)
(*Export Data*)


(* ::Subsection:: *)
(*Export Boundary Data*)


Print["Export External Data Punc"]

DirName="../../InputData/Coord"<>ToString[CoordFlagValue]<>"/r0_over_M_"<>ToString[NumberForm[N[r0],{20,5}]]<>"/eta_"<>ToString[NumberForm[N[\[Eta]],{20,5}]]<>"/m_"<>ToString[m]<>"/nbarMax_"<>ToString[NumberForm[nbarmax]]<>"/N_"<>ToString[N1]<>"/Prec"<>ToString[NumberForm[N[Prec],{20,5}]]
If[!DirectoryQ[DirName],CreateDirectory[DirName]];

fn$\[CapitalPhi]P$Surface=DirName<>"/phi_punc_boundary_nbar"<>ToString[nbarmax]<>".dat";

Print[fn$\[CapitalPhi]P$Surface]

PhysicalData={
"r0\t"<>ToString[r0]<>"\n",
"m\t"<>ToString[m]<>"\n",
"eta\t"<>ToString[\[Eta]]<>"\n"
};

fp=OpenWrite[fn$\[CapitalPhi]P$Surface];
Table[ WriteString[fp,PhysicalData[[i]]], {i,1,Length@PhysicalData}];
Close[fp];


NumericalData={
"nbarmax\t"<>ToString[nbarmax]<>"\n\n",
"N1\t"<>ToString[N1]<>"\n",
"N2\t"<>ToString[N2]<>"\n\n"
};

fp=OpenAppend[fn$\[CapitalPhi]P$Surface];
Table[ WriteString[fp,NumericalData[[i]]], {i,1,Length@NumericalData}];
Close[fp];

fp=OpenAppend[fn$\[CapitalPhi]P$Surface];
WriteString[fp,"#1:i1\t #2:Real(Cheb_PhiP)\t #3:Imag(Cheb_PhiP)\t #4:Real(Cheb_PhiP,sigma)\t #5:Imag(Cheb_PhiP,sigma)\t #6:Real(Cheb_PhiP,y)\t #7:Imag(Cheb_PhiP,y)\n"]
For[i1=0, i1<=N1, i1++,
WriteString[fp,ToString[i1],"\t",
ToString[Re@\[CapitalPhi]m$Punc$Surface$Cheb[[i1+1]],CForm],"\t", ToString[Im@\[CapitalPhi]m$Punc$Surface$Cheb[[i1+1]],CForm],"\t", 
ToString[Re@d\[CapitalPhi]m$Puncd\[Sigma]$Surface$Cheb[[i1+1]],CForm],"\t",ToString[Im@d\[CapitalPhi]m$Puncd\[Sigma]$Surface$Cheb[[i1+1]],CForm],"\t",
ToString[Re@d\[CapitalPhi]m$Puncdy$Surface$Cheb[[i1+1]],CForm],"\t",ToString[Im@d\[CapitalPhi]m$Puncdy$Surface$Cheb[[i1+1]],CForm],"\n" ];]
Close[fp];


(* ::Subsection:: *)
(*Export Bulk Data*)


(*Print["Export External Seff"]*)
fn$Seff=DirName<>"/S_eff_nbar"<>ToString[nbarmax]<>".dat"

fp=OpenWrite[fn$Seff];
Table[ WriteString[fp,PhysicalData[[i]]], {i,1,Length@PhysicalData}];
Close[fp];

fp=OpenAppend[fn$Seff];
Table[ WriteString[fp,NumericalData[[i]]], {i,1,Length@NumericalData}];
Close[fp];

fp=OpenAppend[fn$Seff];
WriteString[fp,"#1:i1\t #2:i2\t #3:Real(Cheb S_eff)\t #4:Imag(Cheb Seff)\n"];
For[i1=0, i1<=N1, i1++,
For[i2=0, i2<=N2, i2++,
WriteString[fp,ToString[i1],"\t",ToString[i2],"\t",
ToString[Re@Seffm$ChebCoeff[[i1+1,i2+1]],CForm],"\t", ToString[Im@Seffm$ChebCoeff[[i1+1,i2+1]],CForm],"\n"
];
];
WriteString[fp,"\n"];
]
Close[fp];


CloseKernels[];
Print["Done"];





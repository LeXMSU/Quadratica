(* ::Package:: *)

BeginPackage["Quadratica`"];


(* ::Subchapter:: *)
(* Public section*)


(* ::Text:: *)
(*ONLY functions inside this section will be seen outside.*)
(**)
(*    	Usage description and error messages*)


(* ::Section:: *)
(*Bosonic case*)


(* ::Subsection:: *)
(*Conversion*)


QuadraticForm::usage =
 "Converts the general case quadratic form to the simmetric one.";
QuadraticForm::wrongsize =
 "The size of the quadratic form and the linear form have to be the same.";

FromNormalForm::usage =
 "Converts the normally ordered quadratic form to the simmetrically ordered one.";

ToNormalForm::usage =
 "Converts the general case quadratic form to the normally ordered one.";

FromAntinormalForm::usage =
 "Converts the antinormally ordered quadratic form to the simmetrically ordered one.";

ToAntinormalForm::usage =
 "Converts the general case quadratic form to the antinormally ordered one.";

FromPqForm::usage =
 "Converts the pq ordered quadratic form to the simmetrically ordered one.";

ToPqForm::usage =
 "Converts the general case quadratic form to the pq ordered one.";

FromQpForm::usage =
 "Converts the qp ordered quadratic form to the simmetrically ordered one.";

ToQpForm::usage =
 "Converts the general case quadratic form to the qp ordered one.";

FromGeneralForm::usage =
 "Converts the generally ordered quadratic form to the simmetrically ordered one.";

ToGeneralForm::usage =
 "Converts the generally case quadratic form to the qp ordered one.";


(* ::Subsection:: *)
(*General symbol evaluation*)


NormalSymbol::usage =
 "Obtains a nromal symbol of operator exponent of quadratic operator.";

ExponentFromNormalSymbol::usage =
  "Calculates the exponent from the normal form.";

AntinormalSymbol::usage =
  "Obtains a nromal symbol of a quadratic form";

ExponentFromAntinormalSymbol::usage =
  "Calculates the exponent from the antinormal form.";


(* ::Subsection:: *)
(*Composition evaluation*)


NormalSymbolsComposition::usage =
  "Evaluates a composition of two normal symbols, i.e. the normal symbol of the product of two operators with known normal symbols.";

ExponentsComposition::usage = 
  "Evaluates a composition of two exponents with quadratic generator as exponent of the third exponent.";


(* ::Subsection:: *)
(*Lindblad equation*)


LindbdladEvolution::usage =
  "Lindblad evolution in terms of covariance matrices and mean values.";

QuantumChannel::usage =
  "Quantum parameters from Lindblad equaiton.";

LindbdladStationaryState::usage =
  "Calcualte a stationary state if the Lindblad equation with a given generator.";


(* ::Subchapter:: *)
(*Private Section*)


Begin["`Private`"];


(* ::Section:: *)
(*Standard matrices definitions*)


(* The matrix of the simplectic form *)
J[n_]:=KroneckerProduct[({
 {0, -1},
 {1, 0}
}),IdentityMatrix[n]];

(* The matrix E changes Subscript[a, i] with Subscript[a, i]^+ , i.e. E a = ( Subscript[a, 1]^+, ..., Subscript[a, n]^+, Subscript[a, 1], ..., Subscript[a, n])^T*)
EE[n_]:=KroneckerProduct[({
 {0, 1},
 {1, 0}
}),IdentityMatrix[n]];

(* Canonical isomorphism *)
UU[n_]:=KroneckerProduct[1/Sqrt[2] ({
 {1, I},
 {1, -I}
}),IdentityMatrix[n]];


(* ::Section:: *)
(*Bosonic case*)


(* ::Subsection::Closed:: *)
(*Simmetrization*)


(* Symmetrize the general quadratic form *)
QuadraticToSymmetricForm[x_]:=
  Module[{K=x[[1]], g=x[[2]], c=x[[3]], n, Jn}, 
	n=Length[K]/2 ; Jn=J[n];
	{1/2 (K+Transpose[K]), g, c+1/2 Tr[K.Jn]}
  ];

(* Symmetrize the form in the exponent of the normal symbol *)
SymmetrizeNormalForm[x_]:= 
  Module[{R, q, s, n}, 
	{R, q, s} =  x;
	n = Length[R]/2;
	{1/2 (R+Transpose[R]), q, s}
  ];


(* ::Subsection::Closed:: *)
(*Conversion*)


(* Trivial constructor of the quadratic form *)
QuadraticForm[x_]/;
  If[RightSizeQ[x],True,Message[QuadraticForm::wrongsize];False] :=
	QuadraticToSymmetricForm[x];

(* Converts quadratic form to normally ordered form *)
ToNormalForm[x_]:= 
  Module[{K, g, c, n, En}, 
	{K, g, c} = QuadraticToSymmetricForm[x];
	n = Length[K]/2 ; En = EE[n];
    {K, g, c+1/2 Tr[En.K]}
  ];

(* Obtains quadratic form from normal form *)
FromNormalForm[x_]:=
  Module[{K, g, c, n, En}, 
	{K, g, c} = SymmetrizeNormalForm[x];
	n = Length[K]/2 ; En = EE[n];
	{K, g, c-1/2 Tr[En.K]}
  ];

(* Converts quadratic form to antinormally ordered form *)
ToAntinormalForm[x_]:= 
  Module[{K, g, c, n, En}, 
	{K, g, c} = QuadraticToSymmetricForm[x];
	n = Length[K]/2 ; En = -EE[n];
    {K, g, c+1/2 Tr[En.K]}
  ];

(* Obtains quadratic form from antinormal form *)
FromAntinormalForm[x_]:=
  Module[{K, g, c, n, En}, 
	{K, g, c} = SymmetrizeNormalForm[x];
	n = Length[K]/2 ; En = -EE[n];
	{K, g, c-1/2 Tr[En.K]}
  ];

(* Converts quadratic form to pq ordered form *)
ToPqForm[x_]:= 
  Module[{K, g, c, n, En}, 
	{K, g, c} = QuadraticToSymmetricForm[x];
	n = Length[K]/2 ; En = EE[n].J[n];
    {K, g, c+1/2 Tr[En.K]}
  ];

(* Obtains quadratic form from pq form *)
FromPqForm[x_]:=
  Module[{K, g, c, n, En}, 
	{K, g, c} = SymmetrizeNormalForm[x];
	n = Length[K]/2 ; En = EE[n].J[n];
	{K, g, c-1/2 Tr[En.K]}
  ];

(* Converts quadratic form to qp ordered form *)
ToQpForm[x_]:= 
  Module[{K, g, c, n, En}, 
	{K, g, c} = QuadraticToSymmetricForm[x];
	n = Length[K]/2 ; En = J[n].EE[n];
    {K, g, c+1/2 Tr[En.K]}
  ];

(* Obtains quadratic form from qp form *)
FromQpForm[x_]:=
  Module[{K, g, c, n, En}, 
	{K, g, c} = SymmetrizeNormalForm[x];
	n = Length[K]/2 ; En = J[n].EE[n];
	{K, g, c-1/2 Tr[En.K]}
  ];

(* Converts quadratic form to general ordered form *)
ToGeneralForm[x_, S_]:= 
  Module[{K, g, c, n, En}, 
	{K, g, c} = QuadraticToSymmetricForm[x];
	n = Length[K]/2 ; En = S.EE[n].Transpose[S];
    {K, g, c+1/2 Tr[En.K]}
  ];

(* Obtains quadratic form from general form *)
FromGeneralForm[x_, S_]:=
  Module[{K, g, c, n, En}, 
	{K, g, c} = SymmetrizeNormalForm[x];
	n = Length[K]/2 ; En = S.EE[n].Transpose[S];
	{K, g, c-1/2 Tr[En.K]}
  ];


(* ::Subsection::Closed:: *)
(*Tests*)


(* Checks if quadratic and linear forms are of the same size *)
RightSizeQ[x_]:= 
  (SquareMatrixQ[x[[1]]]
	&& (Length[x[[1]]] == Length[x[[2]]]));

(* Checks if the qadratic form is hermitian *)
QuadraticHermitianQ[x_]:=
  Module[{K, g, c, n, En},
	{K, g, c} = QuadraticToSymmetricForm[x];
	n = Length[K]/2;
	En = EE[n];
	(c \[Element] Reals) 
		&& (Chop[Norm[En.Conjugate[K].En-K,"Frobenius"]] == 0) 
		&& (Chop[Norm[En.Conjugate[g]-g,"Frobenius"]] == 0)
  ];

(* Checks if the qadratic form is positive hermitian one *)
QuadraticPositiveQ[x_]:=
  Module[{K, g, c, n, En},
	{K, g, c} = QuadraticToSymmetricForm[x];
	n = Length[K]/2;
	En = EE[n];
	(QuadraticHermitianQ[x]) 
		&& (PositiveDefiniteMatrixQ[K.En])
  ];

(*
UnitarySymbolQ[x_]:={};
*)


(* ::Subsection:: *)
(*Order of variables*)


(* Order of tensor product *)
JTransform[n_]:=UnitVector[2 n,#]&/@Flatten[Table[{i,n+i},{i,1, n}]];

(* Order of Hilbert spaces *)
HilbertSpacesOrder[traceOff_,n_]:=KroneckerProduct[UnitVector[n,#]&/@(traceOff~Join~Complement[Range[n],traceOff]),IdentityMatrix[2]];

(* Linear transform to the variables better ordered for trace operation *)
TraceTransform[traceOff_, n_]:=ArrayFlatten[({
 {Inverse[JTransform[Length[traceOff]]], 0},
 {0, Inverse[JTransform[n - Length[traceOff]]]}
})].HilbertSpacesOrder[traceOff,n].JTransform[n];

(* Transfrom of matrix symbol exponent better ordered for trace *)
NormalSymbolTraceTransform[x_, traceOff_]:=
  Module[{R, q, s, n,  T},
	{R, q, s} = x;
	n = Length[R]/2 ;
	T = Inverse[TraceTransform[traceOff, n]];
	{Transpose[T]R.T, T.q, s} 
  ];


(* ::Subsection:: *)
(*General symbol evaluation*)


(* ::Subsubsection:: *)
(*Normal symbol evaluation*)


(* The quadratic form in the exponet of the normal symbol *)
RNormalSymbol[K_]:=
  Module[{n, Jn, En, G},
	n = Length[K]/2 ;
	Jn = J[n]; En = EE[n];
	G = MatrixExp[K.Jn];
	(G-IdentityMatrix[2n]).Inverse[1/2 (En+Jn)-1/2 (En-Jn).G]
  ];

(* The linear form in the exponet of the normal symbol *)
qNormalSymbol[K_,g_]:=
  Module[{n, Jn, En},
	n=Length[K]/2 ;
	Jn=J[n]; En=EE[n];
	Inverse[MatrixFunction[
		\[Piecewise]{
		  {1, #==0},
		  {((# E^#)/(E^#-1)), True}
         }&
	, K.Jn]-1/2 K.(En + Jn)].g
  ];

(* The scalar term in the exponet of the normal symbol if the initial scalar term is zero *)
s0NormalSymbol[K_, g_]:= 
  Module[{n, Jn, En, \[CapitalSigma]m, \[CapitalSigma]p},
	n = Length[K]/2 ;
	Jn=J[n]; En=EE[n];
	\[CapitalSigma]m = 1/2 (En +Jn); \[CapitalSigma]p=1/2 (En -Jn);
	Integrate[1/2 Transpose[qNormalSymbol[K t,g t]].\[CapitalSigma]m.K.\[CapitalSigma]m.qNormalSymbol[K t,g t]+ 1/2 Tr[\[CapitalSigma]m.K.(IdentityMatrix[2 n] + \[CapitalSigma]p.RNormalSymbol[K t])] + Transpose[g].\[CapitalSigma]p.qNormalSymbol[K t,g t],{t,0,1}][[1]][[1]]
  ];

(* Evaluating the normal symbol of  general qaudratic form exponent *)
NormalSymbol[x_]:=
  Module[{K, g, c, n, Jn, En, s},
	{K, g, c} = QuadraticToSymmetricForm[x];
	n = Length[K]/2 ;
	Jn = J[n]; En=EE[n];
	s = c + s0NormalSymbol[K, g];
	{RNormalSymbol[K], qNormalSymbol[K, g], s}
  ];

(* Evaluating exponent by the normal symbol *)
ExponentFromNormalSymbol[x_]:=
  Module[{R, q, s, n, Jn, En, \[CapitalSigma]m, \[CapitalSigma]p, K, g, c},
	{R, q, s} = SymmetrizeNormalForm[x];
	n = Length[R]/2 ;
	Jn = J[n]; En=EE[n];
	\[CapitalSigma]m = 1/2 (En +Jn); \[CapitalSigma]p=1/2 (En -Jn);
	K = -MatrixFunction[Log, Inverse[IdentityMatrix[2 n]+R.\[CapitalSigma]p].(IdentityMatrix[2 n]+R.\[CapitalSigma]m)].Jn;
	g =  K.Inverse[R].q;
	c = s - s0NormalSymbol[K, g];
	{K, q, c}
  ];


(* ::Subsubsection::Closed:: *)
(*Antinormal symbol evaluation*)


(* The quadratic form in the exponet of the antinormal symbol *)
RAntinormalSymbol[K_]:=
  Module[{n, Jn, En, G},
	n = Length[K]/2 ;
	Jn = J[n]; En = -EE[n];
	G = MatrixExp[K.Jn];
	(G-IdentityMatrix[2n]).Inverse[1/2 (En+Jn)-1/2 (En-Jn).G]
  ];

(* The linear form in the exponet of the antinormal symbol *)
qAntinormalSymbol[K_,g_]:=
  Module[{n, Jn, En},
	n=Length[K]/2 ;
	Jn=J[n]; En=-EE[n];
	Inverse[MatrixFunction[
		\[Piecewise]{
		  {1, #==0},
		  {((# E^#)/(E^#-1)), True}
         }&
	, K.Jn]-1/2 K.(En + Jn)].g
  ];

(* The scalar term in the exponet of the antinormal symbol *)
s0AntinormalSymbol[K_, g_]:= 
  Module[{n, Jn, En, \[CapitalSigma]m, \[CapitalSigma]p},
	n = Length[K]/2 ;
	Jn=J[n]; En = -EE[n];
	\[CapitalSigma]m = 1/2 (En +Jn); \[CapitalSigma]p=1/2 (En -Jn);
	Integrate[1/2 Transpose[qAntinormalSymbol[K t,g t]].\[CapitalSigma]m.K.\[CapitalSigma]m.qAntinormalSymbol[K t,g t]+ 1/2 Tr[\[CapitalSigma]m.K.(IdentityMatrix[2 n] + \[CapitalSigma]p.RAntinormalSymbol[K t])] + Transpose[g].\[CapitalSigma]p.qAntinormalSymbol[K t,g t],{t,0,1}][[1]][[1]]
  ];

(* Evaluating the antinormal symbol of  general qaudratic form exponent *)
AntinormalSymbol[x_] :=
  Module[{K, g, c, n, Jn, En, s},
	{K, g, c} = QuadraticToSymmetricForm[x];
	n = Length[K]/2 ;
	Jn = J[n]; En=-EE[n];
	s = c + s0AntinormalSymbol[K, g];
	{RAntinormalSymbol[K], qAntinormalSymbol[K, g], s}
  ];

(* Evaluating exponent by the antinormal symbol *)
ExponentFromAntinormalSymbol[x_]:=
  Module[{R, q, s, n, Jn, En, \[CapitalSigma]m, \[CapitalSigma]p, K, g, c},
	{R, q, s} = SymmetrizeNormalForm[x];
	n = Length[R]/2 ;
	Jn = J[n]; En=-EE[n];
	\[CapitalSigma]m = 1/2 (En +Jn); \[CapitalSigma]p=1/2 (En -Jn);
	K = -MatrixFunction[Log, Inverse[IdentityMatrix[2 n]+R.\[CapitalSigma]p].(IdentityMatrix[2 n]+ R.\[CapitalSigma]m)].Jn;
	q = R.Inverse[K].g;
	c = s - s0AntinormalSymbol[K, g];
	{K, q, c}
  ];


(* ::Subsubsection::Closed:: *)
(*pq - symbol evaluation*)


(* The quadratic form in the exponet of the pq-symbol *)
RPqSymbol[K_]:=
  Module[{n, Jn, En, G},
	n = Length[K]/2 ;
	Jn = J[n]; En = EE[n].Jn;
	G = MatrixExp[K.Jn];
	(G-IdentityMatrix[2n]).Inverse[1/2 (En+Jn)-1/2 (En-Jn).G]
  ];

(* The linear form in the exponet of the pq-symbol *)
qPqSymbol[K_,g_]:=
  Module[{n, Jn, En},
	n=Length[K]/2 ;
	Jn=J[n]; En=EE[n].Jn;
	Inverse[MatrixFunction[
		\[Piecewise]{
		  {1, #==0},
		  {((# E^#)/(E^#-1)), True}
         }&
	, K.Jn]-1/2 K.(En + Jn)].g
  ];

(* The scalar term in the exponet of the pq-symbol *)
s0PqSymbol[K_, g_]:= 
  Module[{n, Jn, En, \[CapitalSigma]m, \[CapitalSigma]p},
	n = Length[K]/2 ;
	Jn=J[n]; En = EE[n].Jn;
	\[CapitalSigma]m = 1/2 (En +Jn); \[CapitalSigma]p=1/2 (En -Jn);
	Integrate[1/2 Transpose[qPqSymbol[K t,g t]].\[CapitalSigma]m.K.\[CapitalSigma]m.qPqSymbol[K t,g t]+ 1/2 Tr[\[CapitalSigma]m.K.(IdentityMatrix[2 n] + \[CapitalSigma]p.RPqSymbol[K t])] + Transpose[g].\[CapitalSigma]p.qPqSymbol[K t,g t],{t,0,1}][[1]][[1]]
  ];

(* Evaluating the pq-symbol of  general qaudratic form exponent *)
PqSymbol[x_] :=
  Module[{K, g, c, n, Jn, En, s},
	{K, g, c} = QuadraticToSymmetricForm[x];
	n = Length[K]/2 ;
	Jn = J[n]; En=EE[n].Jn;
	s = c + s0PqSymbol[K, g];
	{RPqSymbol[K], qPqSymbol[K, g], s}
  ];

(* Evaluating exponent by the pq-symbol *)
ExponentFromPqSymbol[x_]:=
  Module[{R, q, s, n, Jn, En, \[CapitalSigma]m, \[CapitalSigma]p, K, g, c},
	{R, q, s} = SymmetrizeNormalForm[x];
	n = Length[R]/2 ;
	Jn = J[n]; En = EE[n].Jn;
	\[CapitalSigma]m = 1/2 (En +Jn); \[CapitalSigma]p=1/2 (En -Jn);
	K = -MatrixFunction[Log, Inverse[IdentityMatrix[2 n]+R.\[CapitalSigma]p].(IdentityMatrix[2 n]+ R.\[CapitalSigma]m)].Jn;
	q = R.Inverse[K].g;
	c = s - s0PqSymbol[K, g];
	{K, q, c}
  ];


(* ::Subsubsection::Closed:: *)
(*qp - symbol evaluation*)


(* The quadratic form in the exponet of the qp-symbol *)
RQpSymbol[K_]:=
  Module[{n, Jn, En, G},
	n = Length[K]/2 ;
	Jn = J[n]; En = -EE[n].Jn;
	G = MatrixExp[K.Jn];
	(G-IdentityMatrix[2n]).Inverse[1/2 (En+Jn)-1/2 (En-Jn).G]
  ];

(* The linear form in the exponet of the qp-symbol *)
qQpSymbol[K_,g_]:=
  Module[{n, Jn, En},
	n=Length[K]/2 ;
	Jn=J[n]; En = -EE[n].Jn;
	Inverse[MatrixFunction[
		\[Piecewise]{
		  {1, #==0},
		  {((# E^#)/(E^#-1)), True}
         }&
	, K.Jn]-1/2 K.(En + Jn)].g
  ];

(* The scalar term in the exponet of the qp-symbol *)
s0QpSymbol[K_, g_]:= 
  Module[{n, Jn, En, \[CapitalSigma]m, \[CapitalSigma]p},
	n = Length[K]/2 ;
	Jn=J[n]; En = -EE[n].Jn;
	\[CapitalSigma]m = 1/2 (En +Jn); \[CapitalSigma]p=1/2 (En -Jn);
	Integrate[1/2 Transpose[qQpSymbol[K t,g t]].\[CapitalSigma]m.K.\[CapitalSigma]m.qQpSymbol[K t,g t]+ 1/2 Tr[\[CapitalSigma]m.K.(IdentityMatrix[2 n] + \[CapitalSigma]p.RQpSymbol[K t])] + Transpose[g].\[CapitalSigma]p.qQpSymbol[K t,g t],{t,0,1}][[1]][[1]]
  ];

(* Evaluating the qp-symbol of  general qaudratic form exponent *)
QpSymbol[x_] :=
  Module[{K, g, c, n, Jn, En, s},
	{K, g, c} = QuadraticToSymmetricForm[x];
	n = Length[K]/2 ;
	Jn = J[n]; En = -EE[n].Jn;
	s = c + s0QpSymbol[K, g];
	{RQpSymbol[K], qQpSymbol[K, g], s}
  ];

(* Evaluating exponent by the qp-symbol *)
ExponentFromQpSymbol[x_]:=
  Module[{R, q, s, n, Jn, En, \[CapitalSigma]m, \[CapitalSigma]p, K, g, c},
	{R, q, s} = SymmetrizeNormalForm[x];
	n = Length[R]/2 ;
	Jn = J[n]; En = -EE[n].Jn;
	\[CapitalSigma]m = 1/2 (En +Jn); \[CapitalSigma]p=1/2 (En -Jn);
	K = -MatrixFunction[Log, Inverse[IdentityMatrix[2 n]+R.\[CapitalSigma]p].(IdentityMatrix[2 n]+ R.\[CapitalSigma]m)].Jn;
	q = R.Inverse[K].g;
	c = s - s0QpSymbol[K, g];
	{K, q, c}
  ];


(* ::Subsubsection::Closed:: *)
(*Linear normal symbol evaluation*)


(* The quadratic form in the exponet of the linear normal symbol *)
RLinearNormalSymbol[K_, S_]:=
  Module[{n, Jn, En, G},
	n = Length[K]/2 ;
	Jn = J[n]; En = S.EE[n].Transpose[S];
	G = MatrixExp[K.Jn];
	(G-IdentityMatrix[2n]).Inverse[1/2 (En+Jn)-1/2 (En-Jn).G]
  ];

(* The linear form in the exponet of the linear normal symbol *)
qLinearNormalSymbol[K_,g_, S_]:=
  Module[{n, Jn, En},
	n=Length[K]/2 ;
	Jn = J[n]; En = S.EE[n].Transpose[S];
	Inverse[MatrixFunction[
		\[Piecewise]{
		  {1, #==0},
		  {((# E^#)/(E^#-1)), True}
         }&
	, K.Jn]-1/2 K.(En + Jn)].g
  ];

(* The scalar term in the exponet of the linear normal symbol *)
s0LinearNormalSymbol[K_, g_, S_]:= 
  Module[{n, Jn, En, \[CapitalSigma]m, \[CapitalSigma]p},
	n = Length[K]/2 ;
	Jn = J[n]; En = S.EE[n].Transpose[S];
	\[CapitalSigma]m = 1/2 (En +Jn); \[CapitalSigma]p=1/2 (En -Jn);
	Integrate[1/2 Transpose[qLinearNormalSymbol[K t, g t, S]].\[CapitalSigma]m.K.\[CapitalSigma]m.qLinearNormalSymbol[K t,g t, S]+ 1/2 Tr[\[CapitalSigma]m.K.(IdentityMatrix[2 n] + \[CapitalSigma]p.RLinearNormalSymbol[K t, S])] + Transpose[g].\[CapitalSigma]p.qLinearNormalSymbol[K t, g t, S],{t,0,1}][[1]][[1]]
  ];

(* Evaluating the linear normal symbol of  general qaudratic form exponent *)
LinearNormalSymbol[x_, S_] :=
  Module[{K, g, c, n, Jn, En, s},
	{K, g, c} = QuadraticToSymmetricForm[x];
	n = Length[K]/2 ;
	Jn = J[n]; En = S.EE[n].Transpose[S];
	s = c + s0LinearNormalSymbol[K, g, S];
	{RLinearNormalSymbol[K, S], qLinearNormalSymbol[K, g, S], s}
  ];

(* Evaluating exponent by the linear normal symbol *)
ExponentFromLinearNormalSymbol[x_, S_]:=
  Module[{R, q, s, n, Jn, En, \[CapitalSigma]m, \[CapitalSigma]p, K, g, c},
	{R, q, s} = SymmetrizeNormalForm[x];
	n = Length[R]/2 ;
	Jn = J[n]; En = S.EE[n].Transpose[S];
	\[CapitalSigma]m = 1/2 (En +Jn); \[CapitalSigma]p=1/2 (En -Jn);
	K = -MatrixFunction[Log, Inverse[IdentityMatrix[2 n]+R.\[CapitalSigma]p].(IdentityMatrix[2 n]+ R.\[CapitalSigma]m)].Jn;
	q = R.Inverse[K].g;
	c = s - s0LinearNormalSymbol[K, g, S];
	{K, q, c}
  ];

(* General linear normal simbol to linear normal simbol transform *)
LinearNormalToLinearNormalSymbol[x_, S1_, S2_] :=
  Module[{R1, q1, s1, n, Jn, En, \[CapitalSigma]p1, \[CapitalSigma]p2, \[CapitalSigma]m1, \[CapitalSigma]m2, G, R2, q2, s2},
	{R1, q1, s1} = SymmetrizeNormalForm[x];
	n = Length[K]/2 ;
	Jn = J[n]; En = EE[n];
	\[CapitalSigma]p1 = 1/2 S1.(En - Jn).Transpose[S1]; \[CapitalSigma]m1 = 1/2 S1.(En + Jn).Transpose[S1];
	\[CapitalSigma]p2 = 1/2 S2.(En - Jn).Transpose[S2]; \[CapitalSigma]m2 = 1/2 S2.(En + Jn).Transpose[S2];
	G = Inverse[IdentityMatrix[2 n] + R1.\[CapitalSigma]p1].(IdentityMatrix[2 n] + R1.\[CapitalSigma]m1);
	R2 = (G - IdentityMatrix[2 n]).Inverse[\[CapitalSigma]m2 - \[CapitalSigma]p2.G]; 
	q2 = R2.Inverse[R1].q1;
	{R2, q2, s2} (*s2 is not yet defined*)
  ];



(* ::Subsubsection::Closed:: *)
(*Composition of exponents and their symbols *)


(* Composition of two symbols of quadratic exponents *)
NormalSymbolsComposition[x_,y_]:=
  Module[{R1, q1, s1, R2, q2, s2, n, Jn, En, Ip,Im,X,Y,R12,q12,s12},
	{R1, q1, s1} = SymmetrizeNormalForm[x];
	{R2, q2, s2} = SymmetrizeNormalForm[y];
	n = Length[R1]/2;
	Jn = J[n]; En = EE[n];
	Ip = 1/2 (IdentityMatrix[2 n] - Jn.En); Im = 1/2 (IdentityMatrix[2 n] + Jn.En);
	X = Im.R1.Ip + Ip.R2.Im + En;
	Y = Ip.R1.Ip + Im.R2.Im - En;
	R12 = Im.R1.Im+Ip.R2.Ip-En-X.Inverse[Y].(Ip.R1.Im + Im.R2.Ip + En);
	q12 = Im.q1+Ip.q2 - X.Inverse[Y].(Ip.q1+Im.q2);
	s12 = s1 + s2 - 1/2 Log[Det[En.Y]];
	{R12, q12, s12}
  ];

(* Composition of two quadratic exponents *)
ExponentsComposition[x_,y_]:=
  ExponentFromNormalSymbol[NormalSymbolsComposition[NormalSymbol[x], NormalSymbol[y]]];


(* ::Subsubsection::Closed:: *)
(*Weyl symbol evaluation*)


(* The quadratic form in the exponet of the Weyl symbol *)
RWeylSymbol[K_]:=
  Module[{n, Jn, G},
	n = Length[K]/2 ;
	Jn = J[n];
	G = MatrixExp[K.Jn];
	-2 (G - IdentityMatrix[2 n]).Inverse[G + IdentityMatrix[2 n]].Jn
  ];

(* The linear form in the exponet of the Weyl symbol *)
qWeylSymbol[K_,g_]:=
  Module[{n, Jn, En},
	n=Length[K]/2 ;
	Jn=J[n];
	MatrixFunction[
		\[Piecewise]{
		  {1, #==0},
		  {2 (E^#-1)/(E^#+1) #^-1, True}
         }&
	, K.Jn].g
  ];

(* The scalar term in the exponet of the Weyl symbol if the initial scalar term is zero *)
s0WeylSymbol[K_, g_]:= 
  Module[{n, En},
	n = Length[K]/2 ;
	En = E[n];
	1/2 Transpose[qWeylSymbol[K,g]].Inverse[RWeylSymbol[K]].qWeylSymbol[K,g] + 1/2 Log[Det[En.RWeylSymbol[K]]]
  ];

(* Evaluating the Weyl symbol of  general qaudratic form exponent *)
WeylSymbol[x_]:=
  Module[{K, g, c, n, Jn, s},
	{K, g, c} = QuadraticToSymmetricForm[x];
	n = Length[K]/2 ;
	Jn = J[n];
	s = c + s0WeylSymbol[K, g];
	{RWeylSymbol[K], qWeylSymbol[K, g], s}
  ];

(* Evaluating exponent by the Weyl symbol *)
ExponentFromWeylSymbol[x_]:=
  Module[{R, q, s, n, Jn, En, K, g, c},
	{R, q, s} = SymmetrizeNormalForm[x];
	n = Length[R]/2 ;
	Jn = J[n]; En = EE[n];
	K = -MatrixFunction[Log, Inverse[IdentityMatrix[2 n] - 1/2 R.Jn].(IdentityMatrix[2 n]+1/2 R.Jn)].Jn;
	g =  K.Inverse[R].q;
	c = s - s0WeylSymbol[K, g];
	{K, q, c}
  ];


(* ::Subsection:: *)
(*Trace - class and density matrix calculations*)


(* ::Subsubsection::Closed:: *)
(*Moments*)


(* Simmetric covariance matrix \[Alpha] *)
CovarianceMatrixOfExponent[x_]:=
  Module[{K, g, c, n, Jn},
	{K, g, c} = QuadraticToSymmetricForm[x];
	n = Length[K]/2 ;
	Jn = J[n]; 
	-(1/2)Jn.MatrixFunction[Coth[#/2]&, K.Jn]
  ];

(* Mean value m *)
MeanValueOfExponent[x_]:=
  Module[{K, g, c, n, Jn},
	{K, g, c} = QuadraticToSymmetricForm[x];
	- Inverse[K].g
  ];


(* ::Subsubsection:: *)
(*Wigner characteristic function*)


(* ::Subsubsection:: *)
(*Trace and partial trace*)


(* Trace of the normal symbol *)
TraceOfNormalSymbol[x_]:=
  Module[{R, q, s, n, En},
	{R, q, s} = SymmetrizeNormalForm[x];
	n = Length[R]/2 ;
	En=EE[n];
	1/Sqrt[Det[En.R]] E^(-(1/2)Transpose[q].Inverse[R].q)
  ];

(* Trace of the exponent *)
TraceOfExponent[x_]:=
  Module[{y},
	y = NormalSymbol[x];
	TraceOfNormalSymbol[y]
  ];

(* Partial trace *)
PartialTraceCalc[R11_, R12_, R22_, q1_, q2_]:=
  Module[{R, q, s, n, En},
	R = R22 - Transpose[R12].Inverse[R11].R12;
	q = q2 - Transpose[R12].Inverse[R11].q1;
	n = Length[R]/2 ;
	En=EE[n];
	{R, q, -(1/2)Transpose[q].Inverse[R].q - 1/2 Log[Det[En.R]]}
  ];

(* Evaluate partial trace *)
PartialTrace[x_,traceOff_]:=
  Module[{R, q, s, n, d, R11, R12, R22, q1, q2},
	d = Length[traceOff];
	{R, q, s} = x;
	n = Length[R]/2 ;
	R11 = R[[1;;d]][[1;;d]];  R12 = R[[1;;d]][[d+1;;n]]; R11 = R[[d+1;;n]][[d+1;;n]];
	q1 = q[[1;;d]]; q2  = q[[d+1;;n]];
	PartialTraceCalc[R11, R12, R22, q1, q2]
  ];


(* ::Subsection:: *)
(*Lindblad equation*)


(* ::Subsubsection:: *)
(*Covariance and mean representation*)


LindbdladEvolution[Gen_, initCond_, t_]:=
  Module[{m0, C0, m, C, H, \[CapitalGamma], f, n, Jn},
	C0 = initCond[[1]]; (* Initial covariance matrix *)
	m0 = initCond[[2]]; (* Initial mean value *)
	(* Parameters of Lindbdal generator *)
	H = Gen[[1]]; \[CapitalGamma] = Gen[[2]]; f = Gen[[3]];
	n = Length[H]/2;
	Jn = J[n]; (* Simplectic form *)
	(* Mean value evolution *)
	m = MatrixExp[ Jn.(I H + 1/2 (Transpose[\[CapitalGamma]]-\[CapitalGamma])) t].m0 + I MatrixFunction[
		\[Piecewise]{
		  {1, #==0},
		  {((E^#-1)/#), True}
         }&, Jn.(I H + 1/2 (Transpose[\[CapitalGamma]]-\[CapitalGamma]))t].Jn.f;
	(* Cocariance matrix evolution *)
	 C = (#[[1;;2 n]].Inverse[#[[2n+1;;4 n]]])&[MatrixExp[ArrayFlatten[({
 {Jn.(I H + 1/2 (Transpose[\[CapitalGamma]]-\[CapitalGamma])), 1/2 Jn.(Transpose[\[CapitalGamma]]+\[CapitalGamma]).Jn},
 {0, (I H - 1/2 (Transpose[\[CapitalGamma]]-\[CapitalGamma])).Jn}
})]t].ArrayFlatten[({
 {C0},
 {IdentityMatrix[2 n]}
})]];
	(* Return covariance matrix and mean value pair: *)	
	{C, m}
  ];


(* ::Subsubsection:: *)
(*Quantum channel representation*)


QuantumChannel[Gen_, t_]:=
  Module[{H, \[CapitalGamma], f, n, Jn, G, \[Alpha], l},
	(* Parameters of Lindbdal generator *)
	H = Gen[[1]]; \[CapitalGamma] = Gen[[2]]; f = Gen[[3]];
	n = Length[H]/2;
	G = MatrixExp[ J[n].(I H + 1/2 (Transpose[\[CapitalGamma]]-\[CapitalGamma])) t];
	{\[Alpha], l} = LindbdladEvolution[Gen, {1/2 EE[n], 0}, t];
  ];


(* ::Subsubsection:: *)
(*Stationary state of a generator*)


(* Calculation of the stationary state by the help of Lyapunov equation solution *)
LindbdladStationaryState[Gen_]:=
  Module[{H, \[CapitalGamma], f, n, Jn, A, B, C, m},
	(* Parameters of Lindbdal generator *)
	H = Gen[[1]]; \[CapitalGamma] = Gen[[2]]; f = Gen[[3]];
	n = Length[H]/2;
	Jn = J[n];
	A= Jn.(I H + 1/2 (Transpose[\[CapitalGamma]]-\[CapitalGamma])) ;
	B = 1/2 Jn.( Transpose[\[CapitalGamma]]+\[CapitalGamma]).Jn;
	C = LyapunovSolve[{A, IdentityMatrix[2 n]},{Transpose[A],IdentityMatrix[2 n]},-B];
	m = - I Inverse[A].Jn.f;
	{C, m}
  ];




(* ::Section:: *)
(*Fermionic case*)


(* ::Subsection:: *)
(*Conversion*)


(* ::Subsection:: *)
(*Normal Symbol Evaluation*)


(* ::Subchapter::Closed:: *)
(*The End*)


End[];


EndPackage[];

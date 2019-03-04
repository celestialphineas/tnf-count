(* ::Package:: *)


(* ::Title:: *)
(*TNF Count*)


(* ::Subtitle:: *)
(*Counting the tetranucleotide frequency of a list of DNA sequence*)


(* :Context: TNFCount` *)
(* :Author: YIN Yehang @ Zhejiang University *)
(* :Summary: Counting the tetranucleotide frequency of a list of DNA sequence *)
(* :Package Version: 1.0 *)
(* :Mathematica Version: 11.1 *)
(* :Copyright: Copyright (c) 2019 YIN Yehang. All rights reserved. *)
(* :History:
	Version 1.0 Mar 4, 2019, YIN Yehang
*)
(* :Keywords:
  FASTA, Tetranucleotide Frequency
*)


(* ::Section:: *)
(*LICENSE - MIT*)


(*
Copyright (c) 2019 YIN Yehang
All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*)


BeginPackage["TNFCount`"]

Unprotect[
  TNF, BidirectionalTNF, CanonicalTNF,
  ContigsTNF, BidirectionalContigsTNF, CanonialContigsTNF,
  TNFSymbols, CanonicalSymbols
]

Begin["Private`"]


PairedSequence[sequence_String] := 
StringReplace[
  StringReverse[sequence],
  {
    "A" -> "T", "C" -> "G",
    "G" -> "C", "T" -> "A"
  }
];

(* Define the canonical form of a sequence *)
CanonicalForm[sequence_String] := 
Sort[{sequence, PairedSequence[sequence]}][[1]];

(* All tetranucleotide possibilities *)
symbols = Apply[StringJoin, Tuples[{"A", "C", "G", "T"}, 4], {1}];

(* Canonical symbols *)
canonicalSymbols = Union[CanonicalForm /@ symbols];

(* Gather the symbols with the same canonical froms *)
canonicalGathering = 
Gather[symbols, CanonicalForm[#1] == CanonicalForm[#2] &] /. 
  Apply[Rule, Transpose[{symbols, Range[Length[symbols]]}], {1}];

TNF[sequence_String] :=
Normal[
  Counts[
    StringPartition[sequence, 4, 1]
  ][[symbols]]
][[All, 2]] /. {_Missing -> 0};

BidirectionalTNF[sequence_String] :=
Plus @@ {TNF[sequence], TNF[PairedSequence[sequence]]};

CanonicalTNF[sequence_String] :=
Total /@ (canonicalGathering /. Apply[Rule, Transpose[{Range[256], TNF[sequence]}], {1}]);

ContigsTNF[sequences_List] :=
Plus @@ (TNF /@ sequences);

BidirectionalContigsTNF[sequences_List] :=
Plus @@ (BidirectionalTNF /@ sequences);

CanonialContigsTNF[sequences_List] :=
Plus @@ (CanonicalTNF /@ sequences);


TNFSymbols[___] := symbols;

CanonicalSymbols[___] := canonicalSymbols;


End[];

Protect[
  TNF, BidirectionalTNF, CanonicalTNF,
  ContigsTNF, BidirectionalContigsTNF, CanonialContigsTNF,
  TNFSymbols, CanonicalSymbols
]

EndPackage[];

SetDirectory[ParentDirectory @ Directory[]];
<< TNFCount`

sequence = "AAAACAA";
sequences = Import["tests/BS.fasta"];

Print @ Timing @ TNF[sequence]
Print @ Timing @ BidirectionalTNF[sequence]
Print @ Timing @ CanonicalTNF[sequence]

Print @ ContigsTNF[sequences]
Print @ BidirectionalContigsTNF[sequences]
Print @ CanonicalContigsTNF[sequences]

Print @ TNFSymbols[]
Print @ CanonicalSymbols[]

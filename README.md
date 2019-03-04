# TNF Count

Tetranucleotide frequency (TNF) counting tool for Mathematica

## Usage

```mathematica
<< TNFCount`

sequence = "AAAACAA";
sequences = Import["BS.fasta"];

Print @ TNF[sequence]
Print @ BidirectionalTNF[sequence]
Print @ CanonicalTNF[sequence]

Print @ ContigsTNF[sequences]
Print @ BidirectionalContigsTNF[sequences]
Print @ CanonialContigsTNF[sequences]

Print @ TNFSymbols[]
Print @ CanonicalSymbols[]
```

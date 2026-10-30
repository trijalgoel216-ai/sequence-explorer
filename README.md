# sequence-explorer

A small Python tool that reads a DNA sequence, analyses it, and translates it into a protein. It is a learning project, built step by step and tested on one gene: human insulin.

## What it does

- `gc_content(seq)`: the percentage of letters that are G or C.
- `complement(seq)` and `reverse_complement(seq)`: the partner DNA strand.
- `translate(dna)`: cuts DNA into codons (groups of 3 letters) and translates them into a protein, stopping at the first stop codon.
- `count_bases(seq)`: counts the A, T, C and G letters.
- Splits the insulin protein into its parts: signal peptide, B chain, C-peptide and A chain.
- Accepts upper-case and lower-case letters, and raises a clear error for invalid input.

## The data

- **Gene:** human insulin mRNA, NCBI accession `NM_000207.3` (saved as `insulin.fasta`).
- **Length:** 465 letters. GC content is 63.9% overall and 56.0% in the first 50 letters.
- **Base counts:** A 91, T 77, C 156, G 141.
- **Protein:** the translation starts at the first `ATG` (letter 60) and ends at a `TAG` stop codon (letters 390 to 392). The protein, preproinsulin, is 110 amino acids long.

| Piece | Length |
|---|---|
| Signal peptide | 24 |
| B chain | 30 |
| RR cut site | 2 |
| C-peptide | 31 |
| KR cut site | 2 |
| A chain | 21 |

## How to run

You need Python 3 (tested with Python 3.11). In a terminal, in this folder:

```
python sequence_tools.py
```

It prints 37 lines of test results and ends with `All tests finished`. Each line was checked against an answer worked out by hand beforehand.

## How it was checked

- The GC content of the first 50 letters was counted by hand (28 of 50, 56%) before it was coded.
- The first 10 and last 6 codons were translated by hand (`MALWMRLLPL` and `LENYCN`) before the translation function was written.
- The reverse complement of the reverse complement of the whole gene equals the original sequence.
- `GAATTC`, the EcoRI cutting site, is its own reverse complement, and `ATGC` is not.
- The B chain has 2 cysteines and the A chain has 4, matching the known structure of insulin.

## What I learned

- I changed one line in `reverse_complement` by mistake, writing `=` where I needed `+`. The function silently returned an empty string, and I only found out in a code review, because I had switched my tests off. Now I keep all my tests running.
- I write down what I expect a program to print before I run it. When the output did not match, for example when my reverse-complement tests were accidentally calling `complement`, the mismatch told me where to look.
- I can explain why the insulin protein is 110 amino acids and not 111: the last codon, TAG, is a stop signal and does not code for an amino acid. I also learned that the cell builds insulin as one long chain and then cuts it into pieces: signal peptide, B chain, C-peptide and A chain.

## Limits

- Only the letters A, C, G and T are accepted; anything else raises an error.
- It is tested on one gene only.
- It assumes the protein starts at the first `ATG`. It does not yet search a sequence for genes by itself.
- The piece positions are written for insulin and are not general.
- The tests print their results and are checked by eye. There is no automated test framework yet.

## Next

- An ORF finder that searches for every possible protein, in all three reading frames and on the reverse strand.
- Automated tests with `assert`.
- Reading FASTA files directly instead of typing the sequence into the code.

## Credits and honesty

- The gene sequence comes from NCBI, a public database.
- The codon table is the standard genetic code. It was provided to me, not written by me.
- This is a learning project. I wrote the functions and the tests myself, with an AI tutor (Claude) that explained concepts, set exercises and reviewed my code.

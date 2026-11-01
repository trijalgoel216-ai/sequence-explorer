# sequence-explorer

A small Python tool that reads a DNA sequence, analyses it, translates it into a protein, and searches both DNA strands for possible genes. It is a learning project, built step by step and tested on two genes: human insulin and human hemoglobin beta.

## What it does

- `gc_content(seq)`: the percentage of letters that are G or C.
- `complement(seq)` and `reverse_complement(seq)`: the partner DNA strand.
- `translate(dna)`: cuts DNA into codons (groups of 3 letters) and translates them into a protein, stopping at the first stop codon.
- `count_bases(seq)`: counts the A, T, C and G letters.
- `read_fasta(path)`: reads a DNA sequence from a FASTA file.
- `orf_from(dna, start)`: checks whether a complete open reading frame (ORF) starts at a given position. An ORF starts with `ATG`, is read in codons, and ends at a stop codon.
- `find_orfs(dna)`: runs `orf_from` at every position and keeps every complete ORF.
- `merge_nested(orfs)`: when several ORFs end at the same stop codon, keeps only the longest (the earliest `ATG`).
- `search_both_strands(dna)`: searches the DNA as given (strand 1) and its partner strand (strand 2), keeping only the longest ORF for each stop codon.
- `list_orfs(dna, min_length)`: keeps only ORFs of at least `min_length` amino acids and sorts them longest first.
- `print_report(dna, min_length)`: prints a readable report.
- Splits the insulin protein into its parts: signal peptide, B chain, C-peptide and A chain.
- Accepts upper-case and lower-case letters. `check_dna(seq)` checks the DNA once at the entrance of the ORF search, so any invalid letter gives one clear error message.

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

## The ORF finder

On the insulin DNA, with a minimum length of 30 amino acids, the report is:

```
ORF report: minimum length 30 amino acids
strand 1 | start 59 | 110 amino acids | MALWMRLLPL...
strand 2 | start 110 | 49 amino acids | MLVQHCSTMP...
Total: 2
```

- The 110-letter ORF on strand 1 is the real insulin protein.
- ORFs that start at a later `ATG` in the same reading frame end at the same stop codon, so they are the same ORF. They are merged: a 106-letter one on strand 1 and a 41-letter one on strand 2 are not listed separately.
- The ORF on strand 2 is shorter. A random stretch of DNA reaches 49 codons without a stop codon about 10% of the time, but 110 codons only about 0.5% of the time, so the 110-letter ORF is far more likely to be a real gene.
- An `ATG` at position 340 is not reported: the sequence ends before any stop codon is reached.
- Positions on strand 2 are counted from the start of strand 2's own reading, not from strand 1.

## A second gene: hemoglobin beta

I ran the tool on human hemoglobin beta (`NM_000518.5`, saved as `hbb.fasta`, 628 letters). NCBI's own annotation says the coding sequence runs from letter 51 to letter 494 and makes a protein of 147 amino acids. The tool found the same ORF, at start 50 (counting from 0), with 147 amino acids, starting `MVHLTPEEKS`.

A 92-amino-acid ORF in the same reading frame (a later `ATG` inside the same ORF) is merged into it. The report also lists an ORF of 39 amino acids in a different frame, which is most likely chance: a random stretch of DNA reaches 39 codons without a stop codon about 15% of the time, but 147 codons only about 0.09% of the time. With a minimum length of 100, only the real ORF remains.

## How to run

You need Python 3 (tested with Python 3.11). In a terminal, in this folder:

```
python sequence_tools.py
```

It runs 78 automated checks (`assert` statements) at the bottom of `sequence_tools.py`, prints two sample ORF reports, and ends with `All tests passed`. If any check fails, the program stops and shows which line failed. Each expected answer was worked out by hand before it was written down.

## How it was checked

- The GC content of the first 50 letters was counted by hand (28 of 50, 56%) before it was coded.
- The first 10 and last 6 codons were translated by hand (`MALWMRLLPL` and `LENYCN`) before the translation function was written.
- The reverse complement of the reverse complement of the whole gene equals the original sequence.
- `GAATTC`, the EcoRI cutting site, is its own reverse complement, and `ATGC` is not.
- The B chain has 2 cysteines and the A chain has 4, matching the known structure of insulin.
- Small DNA strings with known answers test the ORF functions, including an ORF that is only found on strand 2.
- On hemoglobin beta, the ORF the tool found matches NCBI's annotation: the same start, the same length, and the same first letters.
- When nested ORFs were merged, I wrote down which test outputs should change before running the program. Only those lines changed, and everything else stayed identical.
- The tests are 78 automated `assert` checks. I also changed one expected value on purpose to see the test fail, to confirm that a wrong answer cannot pass silently.
- The hemoglobin result is compared with the protein sequence printed in NCBI's own record, letter for letter.

## What I learned

- I changed one line in `reverse_complement` by mistake, writing `=` where I needed `+`. The function silently returned an empty string, and I only found out in a code review, because I had switched my tests off. Now I keep all my tests running.
- I write down what I expect a program to print before I run it. When the output did not match, for example when my reverse-complement tests were accidentally calling `complement`, the mismatch told me where to look.
- I can explain why the insulin protein is 110 amino acids and not 111: the last codon, TAG, is a stop signal and does not code for an amino acid. I also learned that the cell builds insulin as one long chain and then cuts it into pieces: signal peptide, B chain, C-peptide and A chain.

## Limits

- Only the letters A, C, G and T are accepted. Every function raises a clear error for anything else, including `N` (unknown base), so DNA files with unknown bases are not supported yet.
- It is tested on two genes only (insulin and hemoglobin beta).
- Only complete ORFs are reported (they must end at a stop codon). Only ORFs that share a stop codon are merged, so ORFs in different reading frames can still overlap.
- The piece positions of the insulin protein are written for insulin and are not general.
- The tests are plain `assert` checks that the file runs itself. They are not a test framework such as pytest, and they live in the same file as the code.

## Next

- Move the tests into their own file.
- Run it on more genes and compare each result with NCBI's annotation.

## Credits and honesty

- The gene sequences come from NCBI, a public database.
- The codon table is the standard genetic code. It was provided to me, not written by me.
- This is a learning project. I wrote the functions and the tests myself, with an AI tutor (Claude) that explained concepts, set exercises and reviewed my code.

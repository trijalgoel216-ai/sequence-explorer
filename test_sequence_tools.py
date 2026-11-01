"""Tests for sequence_tools. Run with:  python test_sequence_tools.py"""

from sequence_tools import (
    PAIRS,
    codon_table,
    gc_content,
    complement,
    reverse_complement,
    translate,
    count_bases,
    check_dna,
    orf_from,
    find_orfs,
    merge_nested,
    search_both_strands,
    list_orfs,
    print_report,
    read_fasta,
)

# ---- INSULIN DATA (used by the tests) ----
line1 = "AGCCCTCCAGGACAGGCTGCATCAGAAGAGGCCATCAAGCAGATCACTGTCCTTCTGCCATGGCCCTGTG"
line2 = "GATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCCTTTGTGAAC"
line3 = "CAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTCTTCTACA"
line4 = "CACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGC"
line5 = "AGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACC"
line6 = "AGCATCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTAGACGCAGCCCGCAGGCAGCCCCACACCCG"
line7 = "CCGCCTCCTGCACCGAGAGAGATGGAATAAAGCCCTTGAACCAGC"
full_sequence = line1 + line2 + line3 + line4 + line5 + line6 + line7

start_position = full_sequence.find("ATG")
gene = full_sequence[start_position:]
insulin_protein = translate(gene)

signal_peptide = insulin_protein[0:24]
b_chain = insulin_protein[24:54]
rr = insulin_protein[54:56]
c_peptide = insulin_protein[56: 87]
kr = insulin_protein[87:89]
a_chain = insulin_protein[89:110]

# ---- TESTS ----

def raises_value_error(function, argument):
    try:
        function(argument)
    except ValueError:
        return True
    return False

#GC CONTENT
assert gc_content("ATGC") == 50.0
assert gc_content("AAAA") == 0.0
assert gc_content("GGCC") == 100.0
assert gc_content("AGCCCTCCAGGACAGGCTGCATCAGAAGAGGCCATCAAGCAGATCACTGT") == 56.0
assert gc_content(full_sequence) == 63.9
assert len(full_sequence) == 465
assert gc_content("atgc") == 50.0

#COMPLEMENT
assert complement("ATGC") == "TACG"
assert complement("AAAA") == "TTTT"
assert complement("GGCC") == "CCGG"
assert complement("atgc") == "TACG"
assert raises_value_error(complement, "ATXC")

#REVERSE COMPLEMENT
assert reverse_complement("ATGC") == "GCAT"
assert reverse_complement("AAAC") == "GTTT"
assert reverse_complement("GGCC") == "GGCC"
assert reverse_complement(reverse_complement(full_sequence)) == full_sequence
assert reverse_complement("GAATTC") == "GAATTC"
assert reverse_complement("ATGC") != "ATGC"
assert reverse_complement("atgc") == "GCAT"
assert raises_value_error(reverse_complement, "ATXC")

#CODON TABLE
assert codon_table["GCC"] == "A"
assert codon_table["TAG"] == "*"
assert codon_table["CTG"] == "L"

#TRANSLATE
assert translate("atggccctg") == "MAL"
assert raises_value_error(translate, "ATGGC")

#INSULIN
assert len(insulin_protein) == 110
assert signal_peptide == "MALWMRLLPLLALLALWGPDPAAA"
assert b_chain == "FVNQHLCGSHLVEALYLVCGERGFFYTPKT"
assert rr == "RR"
assert c_peptide == "EAEDLQVGQVELGGGPGAGSLQPLALEGSLQ"
assert kr == "KR"
assert a_chain == "GIVEQCCTSICSLYQLENYCN"
assert b_chain.count("C") == 2
assert a_chain.count("C") == 4

#COUNT BASES
assert count_bases(full_sequence) == {"A": 91, "T": 77, "C": 156, "G": 141}
assert raises_value_error(count_bases, "ATXC")

#ORF_FROM
assert orf_from("CCATGGCCTAAGG", 2) == "MA"
assert orf_from("CCATGGCCGG", 2) is None
assert orf_from("CCATGGCC", 0) is None
assert orf_from(full_sequence, start_position) == insulin_protein
assert orf_from("ATGGCCGCC", 0) is None
assert orf_from("atggcctaa", 0) == "MA"
assert orf_from(full_sequence, 340) is None

#FIND_ORFS
assert find_orfs("CCATGGCCTAAGG") == {2: "MA"}
assert find_orfs("ATGATGTAA") == {0: "MM", 3: "M"}
assert find_orfs("CCCC") == {}
assert find_orfs(full_sequence)[start_position] == insulin_protein
assert list(find_orfs(full_sequence).keys()) == [59, 71, 441]
assert len(find_orfs(full_sequence)[71]) == 106
assert raises_value_error(find_orfs, "ATGNNNTAA")

#CHECK_DNA
assert check_dna("ATGC") is None
assert raises_value_error(check_dna, "ATGX")

#SEARCH BOTH STRANDS
assert search_both_strands("CCATGGCCTAAGG") == {"strand 1": {2: "MA"}, "strand 2": {}}
assert search_both_strands("TTAGGCCAT") == {"strand 1": {}, "strand 2": {0: "MA"}}
assert search_both_strands("CCCC") == {"strand 1": {}, "strand 2": {}}
result = search_both_strands(full_sequence)
assert list(result["strand 1"].keys()) == [59, 441]
assert list(result["strand 2"].keys()) == [110, 404]
assert len(result["strand 1"][59]) == 110
assert len(result["strand 2"][110]) == 49
assert raises_value_error(search_both_strands, "ATGAAANTAA")

#LIST_ORFS
assert len(list_orfs(full_sequence, 30)) == 2
assert len(list_orfs(full_sequence, 50)) == 1
assert len(list_orfs(full_sequence, 107)) == 1
assert len(list_orfs(full_sequence, 111)) == 0
assert len(list_orfs(full_sequence, 49)) == 2

#MERGE_NESTED
assert merge_nested({0: "MM", 3: "M"}) == {0: "MM"}
assert merge_nested({2: "MA"}) == {2: "MA"}
assert merge_nested({}) == {}

#READ_FASTA
assert len(read_fasta("insulin.fasta")) == 465
assert read_fasta("insulin.fasta") == full_sequence
hbb = read_fasta("hbb.fasta")
assert len(hbb) == 628

#HEMOGLOBIN BETA, checked against NCBI's annotation
ncbi_protein = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
assert list(merge_nested(find_orfs(hbb)).keys()) == [50, 114, 240, 468, 570]
best = list_orfs(hbb, 100)
assert len(best) == 1
assert best[0][0] == 147
assert best[0][1] == "strand 1"
assert best[0][2] == 50
assert best[0][3] == ncbi_protein
assert len(list_orfs(hbb, 30)) == 2

#STRAND SYMMETRY
flipped_hbb = reverse_complement(hbb)
assert search_both_strands(flipped_hbb)["strand 2"] == search_both_strands(hbb)["strand 1"]
assert search_both_strands(flipped_hbb)["strand 1"] == search_both_strands(hbb)["strand 2"]

flipped_insulin = reverse_complement(full_sequence)
assert search_both_strands(flipped_insulin)["strand 2"] == search_both_strands(full_sequence)["strand 1"]
assert search_both_strands(flipped_insulin)["strand 1"] == search_both_strands(full_sequence)["strand 2"]

#HEMOGLOBIN ALPHA 1, checked against NCBI's annotation (CDS 38..466)
hba1 = read_fasta("hba1.fasta")
ncbi_alpha = ("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYF"
              "PHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLL"
              "SHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR")
assert len(hba1) == 577
best_alpha = list_orfs(hba1, 100)
assert len(best_alpha) == 1
assert best_alpha[0][0] == 142
assert best_alpha[0][1] == "strand 1"
assert best_alpha[0][2] == 37
assert best_alpha[0][3] == ncbi_alpha
flipped_alpha = reverse_complement(hba1)
assert search_both_strands(flipped_alpha)["strand 2"] == search_both_strands(hba1)["strand 1"]

print_report(full_sequence, 30)
print_report(hbb, 30)

print("All tests passed")

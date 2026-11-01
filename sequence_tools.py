PAIRS = {"A":"T", "T":"A", "C":"G", "G":"C"}

def gc_content(seq):
    seq = seq.upper()
    seq_length = len(seq)
    g_count = seq.count('G')
    c_count = seq.count('C')
    gc_count = g_count + c_count
    gc_percent = gc_count/seq_length*100
    rounded_gc_percent = round(gc_percent, 1)
    return rounded_gc_percent


line1 = "AGCCCTCCAGGACAGGCTGCATCAGAAGAGGCCATCAAGCAGATCACTGTCCTTCTGCCATGGCCCTGTG"
line2 = "GATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCCTTTGTGAAC"
line3 = "CAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTCTTCTACA"
line4 = "CACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGC"
line5 = "AGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACC"
line6 = "AGCATCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTAGACGCAGCCCGCAGGCAGCCCCACACCCG"
line7 = "CCGCCTCCTGCACCGAGAGAGATGGAATAAAGCCCTTGAACCAGC"
full_sequence = line1 + line2 + line3 + line4 + line5 + line6 + line7


def complement(seq):
    seq = seq.upper()
    result = ""
    for letter in seq:
        if letter in PAIRS:
            result = result + PAIRS[letter]
        else:
            raise ValueError(f"Invalid DNA letter: {letter}")
    return result




def reverse_complement(seq):
    seq = seq.upper()
    result = ""
    for letter in seq:
        if letter in PAIRS:
            result = PAIRS[letter] + result
        else:
            raise ValueError(f"Invalid DNA letter: {letter}")
    return result


codon_table = {
    "TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",
    "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
}



 
def translate(dna):
    dna = dna.upper()
    protein = ""
    for i in range(0, len(dna), 3):
        codon = dna[i:i+3]
        if codon in codon_table:
            acid = codon_table[codon]
            if acid == "*":
                break
            else:
                protein = protein + acid
        else:
            raise ValueError(f"'{codon}' is NOT a Codon")
    return protein

start_position = full_sequence.find("ATG")
gene = full_sequence[start_position:]
insulin_protein = translate(gene)

signal_peptide = insulin_protein[0:24]
b_chain = insulin_protein[24:54]
rr = insulin_protein[54:56]
c_peptide = insulin_protein[56: 87]
kr = insulin_protein[87:89]
a_chain = insulin_protein[89:110]

def count_bases(seq):
    seq = seq.upper()
    bases = {"A":0, "T":0, "C":0, "G":0}
    for base in seq:
        if base in bases:
            bases[base] = bases[base] + 1
        else:
            raise ValueError(f"'{base}' is not a Base.")
    return bases


def orf_from(dna, start):
    dna = dna.upper()
    start_codon = dna[start:start+3]
    if start_codon != "ATG":
        return None
    protein = ""
    for j in range(start, len(dna), 3):
        codon_orf = dna[j:j+3]
        if len(codon_orf) < 3:
            return None
        codon_acid = codon_table[codon_orf]
        if codon_acid == "*":
            return protein
        else:
            protein = protein + codon_acid
    return None 

    
def find_orfs(dna):
    dna = dna.upper()
    check_dna(dna)
    orfs = {}
    for position in range(len(dna)):
        answer = orf_from(dna, position)
        if answer is not None:
            orfs[position] = answer
    return orfs

def search_both_strands(dna):
    strands = {"strand 1": merge_nested(find_orfs(dna)), "strand 2":merge_nested(find_orfs(reverse_complement(dna)))}
    return strands

def list_orfs(dna, min_length):
    data = search_both_strands(dna)
    found = []
    for strand_name in data:
        for position in data[strand_name]:
            protein = data[strand_name][position]
            if len(protein) >= min_length:
                found.append((len(protein), strand_name, position, protein))

    return sorted(found, reverse = True)



def print_report(dna, min_length):
    found_list = list_orfs(dna, min_length)
    print(f"ORF report: minimum length {min_length} amino acids")
    if len(found_list) == 0:
        print(f"No ORFs of at least {min_length} amino acids found.")
        return
    else:
        for len_protein, strand_name, position, protein in found_list:
            print(f"{strand_name} | start {position} | {len_protein} amino acids | {protein[:10]}...")
        print(f"Total: {len(found_list)}")



def read_fasta(path):
    sequence = ""
    with open(path) as file:
        for line in file:
            line = line.strip()
            if not line.startswith(">"):
                sequence = sequence + line
        return sequence


def merge_nested(orfs):
    kept = {}
    seen = {}
    for position in orfs:
        protein = orfs[position]
        end = position + 3 * len(protein) + 3
        if end not in seen:
            seen[end] = True
            kept[position] = protein
    return kept

def check_dna(seq):
    seq = seq.upper()
    for letter in seq:
        if letter not in PAIRS:
            raise ValueError(f"Invalid DNA letter: {letter}")



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

print_report(full_sequence, 30)
print_report(hbb, 30)

print("All tests passed")
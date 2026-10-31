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
    orfs = {}
    for position in range(len(dna)):
        answer = orf_from(dna, position)
        if answer is not None:
            orfs[position] = answer
    return orfs

def search_both_strands(dna):
    strands = {"strand 1": find_orfs(dna), "strand 2":find_orfs(reverse_complement(dna))}
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



# ---- TESTS ----

#GC CONTENT
print(gc_content("ATGC"))
print(gc_content("AAAA"))
print(gc_content("GGCC"))
print(gc_content("AGCCCTCCAGGACAGGCTGCATCAGAAGAGGCCATCAAGCAGATCACTGT"))
print(gc_content(full_sequence))
print(len(full_sequence))
print(gc_content("atgc"))

#COMPLEMENT
print(complement("ATGC"))
print(complement("AAAA"))
print(complement("GGCC"))
print(complement("atgc"))

try:
    print(complement("ATXC"))
except ValueError as error:
    print("Caught a problem:", error)

#REVERSE COMPLEMENT
print(reverse_complement("ATGC"))
print(reverse_complement("AAAC"))
print(reverse_complement("GGCC"))
print(reverse_complement(reverse_complement(full_sequence)) == full_sequence)
print(reverse_complement("GAATTC") == "GAATTC")
print(reverse_complement("ATGC") == "ATGC")
print(reverse_complement("atgc"))

try:
    print(reverse_complement("ATXC"))
except ValueError as error:
    print("Caught a problem:", error)

#CODON TABLE
print(codon_table["GCC"])
print(codon_table["TAG"])
print(codon_table["CTG"])

#TRANSLATE
print(translate("atggccctg"))
try:
    print(translate("ATGGC"))
except ValueError as error:
    print("Caught a problem:", error)

#INSULIN TESTS
print(len(insulin_protein))
print(f"Signal Peptide - {signal_peptide}: {len(signal_peptide)}")
print(f"B chain -{b_chain}: {len(b_chain)}")
print(f"RR Cut Site - {rr}: {len(rr)}")
print(f"C Peptide - {c_peptide}: {len(c_peptide)}")
print(f"KR Cut Site - {kr}: {len(kr)}")
print(f"A Chain - {a_chain}: {len(a_chain)}")
print(f"C's in Chain B: {b_chain.count('C')}")
print(f"C's in Chain A: {a_chain.count('C')}")

#COUNT BASES

print(count_bases(full_sequence))

try:
    print(count_bases("ATXC"))
except ValueError as error:
    print("Caught a problem:", error)

print(orf_from("CCATGGCCTAAGG", 2))
print(orf_from("CCATGGCCGG", 2))
print(orf_from("CCATGGCC", 0))
print(orf_from(full_sequence, start_position) == insulin_protein)
print(orf_from("ATGGCCGCC", 0))
print(orf_from("atggcctaa", 0))


print(find_orfs("CCATGGCCTAAGG"))
print(find_orfs("ATGATGTAA"))
print(find_orfs("CCCC"))
print(find_orfs(full_sequence)[start_position] == insulin_protein)

orfs = find_orfs(full_sequence)
for position in orfs:
    print(position, len(orfs[position]))

print(orf_from(full_sequence, 340))

print(search_both_strands("CCATGGCCTAAGG"))
print(search_both_strands("TTAGGCCAT"))
print(search_both_strands("CCCC"))

result = search_both_strands(full_sequence)
for strand in result:
    for position in result[strand]:
        print(strand, position, len(result[strand][position]))


print(len(list_orfs(full_sequence, 30)))
print(len(list_orfs(full_sequence, 50)))
print(len(list_orfs(full_sequence, 107)))
print(len(list_orfs(full_sequence, 111)))
print(len(list_orfs(full_sequence, 49)))

print_report(full_sequence, 30)
print_report(full_sequence, 111)

print("All tests finished")
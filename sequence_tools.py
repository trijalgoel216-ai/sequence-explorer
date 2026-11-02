import sys

"""sequence_tools: a small toolkit for reading DNA, finding ORFs and translating proteins."""

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
    

def main():
    if len(sys.argv) != 3:
        print("Usage: python sequence_tools.py FILE.fasta MIN_LENGTH")
        return
    else:
        file_name = sys.argv[1]
        try:
            min_length = int(sys.argv[2])
        except ValueError:
            print("Error: MIN_LENGTH must be a whole number, for example 100.")
            return
        
        try:
            print_report(read_fasta(file_name), min_length)
        except FileNotFoundError:
            print(f"Error: file '{file_name}' not found")
        except ValueError as error:
            print("Error:", error)


if __name__ == "__main__":
    main()
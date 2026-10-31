1 - What does .count("G") do?
Counts the number of G in the Sequence

2 - Why do I divide by the length and multiply by 100?
dividing by length says the proportion of gc content in the sequence and multiplying by 100 makes it out of 100.

3 - Why do I use round(..., 1)?
else, it would be a long, mesy decimal number with so many places after decimal.

4-will the whole gene’s GC percentage be higher or lower than the 56% of the first 50 letters?
I didn’t make a real prediction, I only said anything could happen. Next time I’ll pick a direction and write it before I run the code.

5 - Was your guess right?
No real guess, so nothing to check. The whole gene was 63.9%, higher than the first 50 letters (56.0%).

6 - What does 56% for the first 50 letters versus 63.9% for the whole gene tell you?
Some parts are more richer in GC content than others

7 - What is a codon?
A 3 letter code, that codes for amino acids.

8 - What does ATG do?
Its a Start Codon and it starts the protein recipe. ATG also codes for the amino acid methionine (M). It’s both the start signal and the first letter of the protein.

9 - What does a stop codon do?
Stops the protein recipe. A stop codon doesn’t code for any amino acid. It only ends the protein.

10 - Why would a computer be better at translating all 110 codons than you?
It can break the sequence in 3 letter code and map each to particualr amino acid easily and quickly.

11 - Predict this too: do the reverse complement of the reverse complement of full_sequence, and compare it with the original. Print:
reverse_complement(reverse_complement(full_sequence)) == full_sequence
Should it print True or False? Write your guess and a one-line reason.
Complementing twice returns each letter to what it was (A to T to A). Reversing twice returns the order to what it was. So the two steps cancel out, and the result is the original sequence.

12 - What does .find("ATG") do, and why did we need it before translating?
Finds the postion of the start codon ATG, we need it becuase translation starts from thsi codon and protein starts from this codon.find("ATG") gives the position of the first ATG. That works for insulin, but in other genes the first ATG isn’t always the real start.

13 - What does full_sequence[start_codon:] take?
It takes the part of sequence form the start codon till the end.

14 - What does break do in translate?
It breaks the loop when the it encounters the stop codon, cause translation stops.

15 - Why did my length maths give 110 and not 111?
That one last codon is stop codon which doesnt codes any amino acid

16 - Why is a silent wrong answer worse than a crash?
Because we will get wrong answer without knowing it was wrong.

17 - What does .upper() do, and why put it at the start of each function?
So every input is case insensitive. .upper() turns every lowercase letter into a capital. I put it first so the rest of the function only ever sees capitals.

18 - Why do the functions raise the error, while the code that calls them does the try and except?
Functions raise the error cause we are telling specifically to raise the error and also it is in functions that the specific problem occurs and it can't continue, while the code that calls them tries try and except: try literally tries to execute the fucntions successfully while except tells the code what to do when the functions encounters the problems, instead of crashing the terminal/program and keeps on continuing.The function knows what went wrong, but only the code that calls it knows what to do about it: show a message, skip it, or stop. So the function raises, and the caller decides.

19 - Why is it better to define PAIRS once?
Because if we want to do any change in the dictionary, we can do it at once. With one PAIRS, the two functions can never disagree.

20 - orf_from(full_sequence, start_position): will this give the same protein as insulin_protein?
Yes

21- why does the function return None instead of an empty string?
An empty string is a real value, and None means “no value”. The point is this: None says “there is no answer here”, while an empty string would be a real answer that happens to be empty. A function that returns None lets the caller ask “did I get an answer?” and not confuse “no ORF” with “an ORF with an empty protein”

22 - Why was your first condition always True?
Becuause one condtion was always True since a codon cant be all 3 types of stop codon at once and also other types of codon also returned True.with or, one true part makes the whole condition true.

23 - Why does find_orfs call orf_from at every position, and not only at the ATG’s?
orf_from already returns None when there’s no ATG at that position. So find_orfs doesn’t need to hunt for ATGs separately. It can ask every position and keep only the real answers.

24 - In the result dictionary, what does a key mean, and what does a value mean?
Key is position and value is protein

25 - why 340 is missing?
From position 340 the reading runs to the end of the DNA without meeting a stop codon, and ends with a 2-letter leftover, so orf_from returns None.

26 - why can’t I just search the top strand and skip the bottom?
Beacause oen strand makes diff protein, another strand produces diff protein.Genes can sit on either strand, so if I only search strand 1, I would miss every gene that is on strand 2.

27 - What is the difference between a gene and a genome?
Genome is the whole DNA of the humans, and gene is part of the genome which makes the proteins. this holds for any organism, not just humans.

28 - what does a position on strand 2 mean?
These are strand 2's own postions and diff form strand 1
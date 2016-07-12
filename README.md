# HaplotypeAssembly
This projets simulates the procedure to reassemble complementary haplotypes from read sequences.
Since these sequeneces can come from either haplotype, the reads must be processed to determine which haplotype they belong to.

The greedy algorithm sorts all the reads by the first occuring SNP and tries to match overlapping reads to determine which haplotype the reads originated from.
This is much faster than the baseline method which runs in exponential time as it compares all possible haplotypes to the reads.

Check out the [attached pdf](presentation/Haplotype Assembly_Neil Marion.pdf) for a more in-depth explanation.

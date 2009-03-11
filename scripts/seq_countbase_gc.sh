# The following will count the number of non N base pairs in a fasta file
#(
#    (grep -v '>' $1 | wc -m) | tr '\n' -;
#    grep -v '>' $1 | grep -o 'N' | wc -m;
#) | bc
# Now trying to remove the newline from the result
(
    (grep -v '>' $1 | wc -m) | tr '\n' -;
    grep -v '>' $1 | grep -o 'G|C' | wc -m;
) | bc; 
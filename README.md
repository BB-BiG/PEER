# PEER
Pattern Extraction through Entropy Retrieval

The source code 'PEER_DNA.c' is implementation of the algorithm PEER for DNA sequences.

To compile the source code 'PEER_DNA.c' following command can be used

$gcc PEER.c -DMAX=39000000 -o PEER_DNA

The option -DMAX indicates maximum length(in decimal) of the DNA sequences. We have used it up to 39000000 nt.

To run the executable file 'PEER_DNA' following command can be used

$PEER_DNA Input1.txt Input2.txt > output.txt

Files 'Input1.txt' and 'Input2.txt' contains the DNA sequences for which PEER similarity is to determined and  
file 'output.txt' will store the Wait vectors and and PEER similarity along with computation time. 

Note: The source code 'PEER_protein.c' is implementation of the algorithm PEER for protein sequences, compilation and execution 
are similar to that of verion for DNA sequence. 

The Hamming distance is one easy way of matching two bit streams. We use this assignment to understand the complexity of a quantum erasure channel, and with a QBER as high as 25%. Start by

Using arguments (inputs) to main()
Generating random numbers
Estimate the Hamming distance between arrays of random bits, without errors
For your submission, your program will have 3 inputs
N: the total number of bits ~ 1 million (more the better, but start small)
P: the fraction of bits that reach the receiver ~ 0.1
Q: the fraction of bits that have an error ~ 0.25
Imagine two parties, Alice and Bob.
A: transmits 1M bits
B: receives 10K bits, with some % errors
There is a time difference between A and B, since they are at different physical locations. This causes the bit arrays to be offset by an unknown value.

B chooses 2K bits and tells Alice which locations (times) he received the bits, and the corresponding bit values.
A must use the 2K bits to determine the offset.  The remaining 8K bits become a secret key between A and B.

The assignment has four parts
a) setup A and B arrays, create C as a subset of B
b) A will use C (with unknown offset and some % errors), plotting Hamming distance versus offset. Whenever the Hamming distance is a minimum, you have the true offset.
c) show that, in the presence of errors, your ability to determine the location of min(Hdist) improves as length(B') increases.
d) check the randomness of the secret key that you have simulated, and

Use dynamic memory allocation to create an array A[N] random bits of 0 and 1 (since N is an input to main, you will need to allocate/free memory)
Generate a random offset M (e.g. between N and N/10) to create B[N-M] such that B[0] = A[M]
Flip the bits (second row) in B[] array with a probability Q
Select P locations from B[] and store those bit values in a 2xP separate array (location, bit_value). This is C[P]
In quantum key distribution,  Bob tells Alice the locations of B, but not the bit_values, and hence complete secrecy of the key (password)
Now compute the Hamming distance of bit_values in C[] with the values in A[].  Vary the offset (from m=0 to m<N-M)
For better efficiency, use bit-wise logic operators
Compare A and B and calculate the bit error rate BER.
Check that your secret key in B was random (use the NIST tests)
https://github.com/kravietz/nist-sts
https://www.geeksforgeeks.org/hamming-distance-two-strings/

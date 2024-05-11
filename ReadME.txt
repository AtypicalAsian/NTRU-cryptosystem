20-April-2024
Nicolas Nguyen
CS401 - Prof. Duncan Buell
Denison University
_______________________________________________________________________________________________

		            	NTRU Encryption/Decryption
This project implements the NTRU cryptosystem, which is a a lattice-based encryption 
scheme. NTRU stands for "N-th degree TRUncated polynomial ring", referring to the 
mathematical structure it is based on. Unlike some other encryption methods, NTRU 
is resistant to attacks from both quantum and classical computers. Its security 
relies on the difficulty of finding short vectors in certain lattices, making it a 
promising candidate for post-quantum cryptography.

_______________________________________________________________________________________________

    Files                                Description
    ------------------    -------------------------------------------------------------------
    ntru.py   	           NTRU class implementation with encryption/decryption methods
    input.txt   	       File 1 with input paramters
    input2.txt             File 2 with input parameters
    main.py 	           Parses input data and instantiates the NTRU class to perform encryption
                           and decryption
    helper.py              Contains utility functions for ntru.py. Within this module,
                           some functions and methods are sourced and adapted from external 
                           sources (as detailed in the comments).

_______________________________________________________________________________________________
                        
                            USAGE

To perform NTRU encryption/decryption implementation using parameters
from an input file, run the main.py script by executing the command:

    - python main.py <filename>

The folder comes with 2 input files - input.txt and input2.txt - which contains 2 examples of 
the necessary parameters required for the NTRU encryption/decryption process.



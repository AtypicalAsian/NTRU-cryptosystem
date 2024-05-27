# NTRU Cryptosystem

The NTRU (N-th degree truncated polynomial ring) cryptosystem is a lattice-based public-key cryptographic algorithm designed to be secure against both classical and quantum attacks.
Unlike traditional cryptographic methods like RSA and ECC, which are vulnerable to quantum computers, NTRU leverages the hardness of lattice problems, specifically the Shortest Vector Problem (SVP),
which are believed to be resistant to quantum attacks because no efficient quantum algorithms have been discovered to solve them.

NTRU's efficiency in key generation, encryption, and decryption processes makes it particularly well-suited for applications requiring high-speed and
low-latency operations, such as IoT devices and real-time communication systems. As quantum computing technology continues to evolve, the need for robust post-quantum cryptographic solutions becomes increasingly critical.
NTRU stands out as a promising candidate due to its strong security foundations and operational efficiency, positioning it as a key player in the future landscape of secure digital communications.

This project's implementation is based on the research presented by Hoffstein, Pipher, and Silverman in their 1996 paper on <a href = "https://www.ntru.org/f/hps98.pdf"> NTRU cryptography<a/>.

## Project Organization

---

    ├── README.md          <- The top-level README for developers using this project.
    ├── ntru.py            <- NTRU class implementation with encryption/decryption methods
    ├── models             <- Contains utility functions for ntru.py. Within this module, some methods within this module are adapted from external sources (detailed in the comments).
    ├── main.py            <- Parses input file and instantiates the NTRU class
    ├── input1.txt         <- Input parameters
    ├── input2.txt         <- Input parameters

## Dependencies

```bash
$ pip3 install -r requirements.txt
```

## Clone this Repository

```bash
$ git clone https://github.com/AtypicalAsian/NTRU-cryptosystem.git
```

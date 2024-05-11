"""
Nicolas Nguyen
April 16, 2024
main.py
This program instantiates the NTRU class and performs encrytion and decryption on sample data
"""
import sys
from ntru import NTRU

def parse_input(file_name):
    """
    Parses input from a file and returns a dictionary containing parsed values.

    Params:
        file_path (str): The name of the file

    Returns:
        dict: A dictionary containing parsed values with keys 'f', 'g', 'm', and 'r'.
              'f', 'g', 'm', and 'r' are lists representing polynomial coefficients or
              integers representing specific parameters.
    """
    # Split input string by lines
    with open(file_name, 'r') as file:
        input_str = file.read()
    lines = input_str.strip().split('\n')

    # Initialize dictionaries to store parsed values
    parsed = {'f': [], 'g': [], 'm': [], 'r': []}

    for line in lines:
        key, value = line.split('::')
        key = key.strip()
        value = value.strip()

        # Parse comma-separated values into lists
        if key in ['f', 'g', 'm', 'r']:
            parsed[key] = list(map(int, value.split(',')))
        else:
            parsed[key] = int(value)
    return parsed


def main():
    """
    Main Function
    """
    filename = sys.argv[1]
    parsed = parse_input(filename)
    n_ring, p_mod = parsed['N'], parsed['p']
    q_mod, d_bound, f_poly = parsed['q'], parsed['d'], parsed['f']
    g_poly, message, r_poly_blinding = parsed['g'], parsed['m'], parsed['r']
    print("Parameters: N = ",n_ring, "p = ",p_mod, "q = ",q_mod)

    print("f(x)= ", f_poly)
    print("g(x)= ", g_poly)
    print("d   = ", d_bound, "\n\n")

    party_a = NTRU(n_ring,p_mod,q_mod)

    #Generate Public Key: [self.f_p, self.f_q, self.h]
    pub_key = party_a.create_pub_key(f_poly,g_poly, d_bound)

    print("f.p inverse = ", pub_key[0])
    print("f.q inverse = ", pub_key[1], "\n\n")
    print("Public Key: ",pub_key[2])

    party_b = NTRU(n_ring,p_mod,q_mod)
    party_b.set_h(pub_key[2])

    #Encryption - e
    print("\n ----Encryption----")
    print("Orginal Message: ", message)
    print("Blinding Polynomrial r: ", r_poly_blinding)
    encrypted = party_b.encrypt(r_poly_blinding, message)
    print("pr star h: ", encrypted[0])
    print("(pr star h) + m: ", encrypted[1])
    print("Encrypted Message: ", encrypted[2], "\n\n")


    #Decryption -m
    print("----Decryption----")
    decrypted = party_a.decrypt(encrypted[2])
    print("f star e: ", decrypted[0])
    print("f star e (modded down): ", decrypted[1])
    print("center-lift: ", decrypted[2])
    print("f.p star center-lift: ", decrypted[3])
    print("f.p star center-lift (modded down): ", decrypted[4])
    print("Decrypted Message: ", decrypted[5])


if __name__ == "__main__":
    main()

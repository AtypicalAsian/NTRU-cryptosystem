"""
Nicolas Nguyen
April 16, 2024
ntru.py
This program implements the NTRU class with encryption and decryption methods
"""
import sys
import math
import helper


class NTRU:
    """
    NTRU Class - secure cryptosystem post quantum computing
    """
    N = 7
    p = 3
    q = 41
    f = [-1,0,1,1,-1,0,1]
    g = [0,-1,-1,0,1,0,1]
    d = 2
    h = None
    r = None
    f_p = None
    f_q = None
    def __init__ (self, N_ring_order, p_modulus, q_modulus):
        """
        Constructor to initialize variables
        """
        self.N = N_ring_order
        self.p = p_modulus
        self.q = q_modulus

        #initialize a polynomial in N to use to mod down
        self.p_in_N = [-1] + [0] * (self.N - 1) + [1]

    def test_prime(self):
        """
        Check whether N is a prime or not
        """
        if self.N < 2:
            return False
        if self.N == 2:
            return True
        for i in range(2, int(self.N ** 0.5) + 1):
            if self.N % i == 0:
                return False
            return True

    def mod_residue(self, coeff_num, coeff_div, modulus):
        """
        Compute the residue of a polynomial division modulo a given polynomial.

        Paramters:
            coeff_num (list): Coefficients of the polynomial numerator.
            coeff_div (list): Coefficients of the polynomial divisor.
            modulus (list): Coefficients of the polynomial used for modulo operation.

        Returns:
            Coefficients of the residue obtained from the polynomial division modulo operation.
        """
        # Perform polynomial division and obtain the remainder
        _, residue = helper.pol_division(coeff_num, coeff_div)

        # Perform modulo operation on the remainder polynomial
        return helper.poly_mod(residue, modulus)

    def create_pub_key(self,f_poly, g_poly, d_coeff):
        """
        Generate the public key using the provided polynomials and coefficient.

        Parameters:
            f_poly (list): Coefficients of polynomial f.
            g_poly (list): Coefficients of polynomial g.
            d_coeff (list): integer bounds d

        Returns:
            None

        Raises:
            ValueError: If the conditions for generating the public key are not met.
        """
        self.f = f_poly
        self.g = g_poly
        self.d = d_coeff

		# Compute s_f using Extended Euclidean Algorithm
        _, s_f, _ = helper.poly_ext_euclidean(self.f, self.p_in_N)
        # Reduce s_f modulo p and q
        self.f_p = helper.poly_mod(s_f, self.p)
        self.f_q = helper.poly_mod(s_f, self.q)

        # Compute public key h
        self.h = self.mod_residue(helper.pol_multiplication(self.f_q, self.g), self.p_in_N, self.q)

        # Check if the generated public key passes all tests
        if not helper.tern_check(self.g, self.d, self.d):
            print("g must belong to T(d,d)")
            sys.exit()
        if not helper.tern_check(self.f, self.d+1, self.d):
            print("f must belong to T(d+1, d)")
            sys.exit()
        if math.gcd(self.N, self.p) != 1 or math.gcd(self.N, self.q)!=1:
            print("gcd(N,q) != 1 or gcd(N,p) != 1")
            sys.exit()
        if not self.test_prime():
            print("N must be prime!")
            sys.exit()
        if self.q <= (6 * self.d + 1) * self.p:
            print("q must be greater than (6*d+1)*p")
            sys.exit()
        return [self.f_p, self.f_q, self.h]


    def get_pub_key(self):
        """
        Return the public key
        """
        return self.h

    def set_h(self,new_pubkey_val):
        """
        Sets the public key to a given value
        """
        self.h = new_pubkey_val

    def encrypt(self,r_pol,mess):
        """
        Parameters:
            cipher_text (list): Coefficients of the ciphertext polynomial.

        Returns:
            list: Coefficients of the decrypted message polynomial.
        """
        prh = helper.pol_multiplication(helper.pol_multiplication([self.p], r_pol), self.h)
        prh_plus_m = helper.pol_addition(prh, mess)
        encrypted = self.mod_residue(prh_plus_m, self.p_in_N, self.q)
        return [prh, prh_plus_m, encrypted]


    def decrypt(self, cipher_text):
        """
        Decrypts the given ciphertext polynomial using the private key

        Parameters:
            cipher_text (list): Coefficients of the ciphertext polynomial.

        Return:
            list: Coefficients of the decrypted message polynomial.
        """
        # Compute the decrypted polynomial by multiplying with the private key polynomial f
        f_star_e = helper.pol_multiplication(self.f, cipher_text)

        # Reduce the decrypted polynomial modulo q
        f_star_e_mod = self.mod_residue(f_star_e, self.p_in_N, self.q)

        # Centerlift the decrypted polynomial
        center_lift = helper.centerlift(f_star_e_mod, self.q)

        # Multiply the centered polynomial with the polynomial f_p
        fp_star_centerlift = helper.pol_multiplication(self.f_p, center_lift)

        # Reduce the resulting polynomial modulo p_in_N and p
        fp_star_centerlift_mod = self.mod_residue(fp_star_centerlift, self.p_in_N, self.p)

        # Remove any trailing zeroes from the decrypted message
        truncated = helper.remove_zeroes(fp_star_centerlift_mod)

        ret = [-1 if x == 2 else x for x in truncated]
        return [f_star_e,f_star_e_mod,center_lift,fp_star_centerlift,fp_star_centerlift_mod,ret]
    
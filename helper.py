"""
Nicolas Nguyen
April 16, 2024
helper.py
This program contains some helper and utility functions for the NTRU implementation
"""
from fractions import Fraction as frac
from operator import add
from operator import neg

def extended_gcd(num_a,num_b):
    '''
    This function implements the Extended Euclidean Algorithm

    PARAMS:
        num_a: first integer
        num_b: second integer

    RETURN:
        A tuple containing three integers (gcd, x, y), where gcd is the greatest common divisor
        of num_a and num_b, and x and y are integers satisfying the equation ax + by = gcd(a, b)
    '''
    if num_a == 0:
        return num_b, 0, 1

    gcd, intermediate_x1, intermediate_y1 = extended_gcd(num_b % num_a, num_a)


    coeff_x = intermediate_y1 - (num_b//num_a) * intermediate_x1
    coeff_y = intermediate_x1

    return gcd, coeff_x, coeff_y



#======================================= POLYNOMIAL OPERATIONS ====================================
def remove_zeroes(pol):
    """
    This function removes the leading zeroes from a coefficient array

    PARAMETERS:
        pol (list): Coefficient array representing the polynomial.

    RETURNS:
        Coefficient array without leading zeroes
    """
    size = len(pol)
    if size == 0:
        return pol
    i = size - 1
    while i > 0:
        if pol[i] != 0:
            break
        i -= 1
    ret = pol[:i+1]
    return ret


def pol_addition(pol1, pol2):
    """
    This function performs addition operation on two polynomials.

    PARAMETERS:
        pol1 (list): Coefficient array representing the first polynomial.
        pol2 (list): Coefficient array representing the second polynomial.

    RETURNS:
        Coefficient array representing the polynomial obtained by adding pol1 and pol2.
    """
    # Ensure that the coefficient arrays have the same length
    len_pol1, len_pol2 = len(pol1), len(pol2)
    if len_pol1 < len_pol2:
        pol1 += [0] * (len_pol2 - len_pol1)
    elif len_pol1 > len_pol2:
        pol2 += [0] * (len_pol1 - len_pol2)

    # Perform addition of corresponding coefficients
    result = [a + b for a, b in zip(pol1, pol2)]

    # Remove any trailing zero coefficients
    while result and result[-1] == 0:
        result.pop()

    return remove_zeroes(result)

def pol_subtraction(pol1,pol2):
    """
    This function performs subtraction operation on two polynomials.

    PARAMETERS:
        pol1 (list): Coefficient array representing the first polynomial.
        pol2 (list): Coefficient array representing the second polynomial.

    RETURNS:
        list: Coefficient array representing the polynomial obtained by subtracting pol2 from pol1.
    """
    # Ensure that the coefficient arrays have the same length
    len_pol1, len_pol2 = len(pol1), len(pol2)
    if len_pol1 < len_pol2:
        pol1 += [0] * (len_pol2 - len_pol1)
    elif len_pol1 > len_pol2:
        pol2 += [0] * (len_pol1 - len_pol2)

    # Perform subtraction of corresponding coefficients
    pol2 = list(map(neg,pol2))
    res = list(map(add, pol1, pol2))
    return remove_zeroes(res)

def pol_multiplication(pol1,pol2):
    """
    This function performs multiplication operation on two polynomials

    PARAMETERS:
        pol1 (list): Coefficient array representing the first polynomial
        pol2 (list): Coefficient array representing the second polynomial

    RETURNS:
        Coefficient array representing the polynomial obtained by multiplying pol1 and pol2
    """
    # Determine the order of the resulting polynomial
    len_pol1, len_pol2 = len(pol1), len(pol2)
    res_ord = len_pol1 + len_pol2 - 2

    # Initialize the output coefficient array with zeros
    res = [0] * (res_ord + 1)

    # Perform multiplication of coefficients
    for i in range(len_pol1):
        for j in range(len_pol2):
            res[j + i] += pol1[i] * pol2[j]

    # Remove any trailing zero coefficients
    return remove_zeroes(res)



def pol_division(num,denom):
    """
    This function performs polynomial division.

    PARAMETERS:
        pol_num (list): Coefficient array representing the numerator polynomial.
        pol_denom (list): Coefficient array representing the denominator polynomial.

    RETURNS:
        list: Coefficient array representing the quotient polynomial.
        list: Coefficient array representing the remainder polynomial.
    """
    num, denom = list(map(frac,remove_zeroes(num))), list(map(frac,remove_zeroes(denom)))
    deg_num, deg_d = len(num)-1, len(denom)-1
    if deg_num>=deg_d:
        quotient=[0]*(deg_num-deg_d+1)
        while(deg_num>=deg_d and num!=[0]):
            divisor=list(denom)
            [divisor.insert(0,frac(0,1)) for i in range(deg_num-deg_d)]
            quotient[deg_num-deg_d]=num[deg_num]/divisor[len(divisor)-1]
            divisor=list(map(lambda x: x*quotient[deg_num-deg_d],divisor))
            num=pol_subtraction(num,divisor)
            deg_num=len(num)-1
        rem=num
    else:
        quotient=[0]
        rem=num
    return [remove_zeroes(quotient),remove_zeroes(rem)]



def poly_mod(coefficients,modulus):
    """
    Compute the coefficients of the polynomial modulo a given modulus.

    Parameters:
        coefficients (list): Coefficients of the polynomial.
        modulus (int): Modulus value.

    Returns:
        list: Coefficients of the polynomial modulo the given modulus.
    """
    if modulus == 0:
        raise ValueError("Error, diving by 0!")

    return [mod_fraction(coeff, modulus) for coeff in coefficients]


def centerlift(coefficients,modulus):
    """
    Centerlift a polynomial with respect to a modulus.

    Parameters:
        coefficients (list): Coefficients of the polynomial.
        modulus (int): Modulus value.

    Returns:
        list: Centerlifted coefficients of the polynomial.
    """
    upper_limit = modulus / 2
    lower_limit = -upper_limit
    centered_coefficients = []

    for coeff in coefficients:
        if coeff > upper_limit:
            centered_coefficients.append(coeff % (-modulus))
        elif coeff <= lower_limit:
            centered_coefficients.append(coeff % modulus)
        else:
            centered_coefficients.append(coeff)

    return centered_coefficients


def poly_ext_euclidean(coeff_a,coeff_b):
    """
    Computes the extended Euclidean algorithm for polynomials.

    Parameters:
        coeff_a (list): Coefficients of the first polynomial.
        coeff_b (list): Coefficients of the second polynomial.

    Returns:
        gcd_val: Coefficients of the greatest common divisor (GCD) of the two polynomials.
        s_out: Coefficients of the first polynomial after the algorithm.
        t_out: Coefficients of the second polynomial after the algorithm.
    Source: Adapted from https://www.youtube.com/watch?v=S88eItWC1QY
    """
    switch = False
    coeff_a=remove_zeroes(coeff_a)
    coeff_b=remove_zeroes(coeff_b)
    if len(coeff_a)>=len(coeff_b):
        coeff_a1, coeff_b1 = coeff_a, coeff_b
    else:
        coeff_a1, coeff_b1 = coeff_b, coeff_a
        switch = True
    Q,R=[],[]
    while coeff_b1 != [0]:
        [q,r]=pol_division(coeff_a1,coeff_b1)
        Q.append(q)
        R.append(r)
        coeff_a1=coeff_b1
        coeff_b1=r
    S=[0]*(len(Q)+2)
    T=[0]*(len(Q)+2)
    S[0],S[1],T[0],T[1] = [1],[0],[0],[1]
    for i in range(2, len(S)):
        S[i]=pol_subtraction(S[i-2],pol_multiplication(Q[i-2],S[i-1]))
        T[i]=pol_subtraction(T[i-2],pol_multiplication(Q[i-2],T[i-1]))

    gcd_val=R[-2]
    s_out=S[-2]
    t_out=T[-2]

   # Scale GCD and coefficients such that the leading term as coefficient 1
    scale_factor = gcd_val[-1]
    gcd_val = [x / scale_factor for x in gcd_val]
    s_out = [x / scale_factor for x in s_out]
    t_out = [x / scale_factor for x in t_out]

    # Return GCD and coefficients based on the switch flag
    return [gcd_val, t_out, s_out] if switch else [gcd_val, s_out, t_out]


def tern_check(coeff, no_1, no_m1):
    """
    Check if the polynomial is ternary given number of 1 and -1 coefficients

    Parameters:
        coeff (list): Coefficients of the polynomial.
        no_1 (int): Number of ones in the polynomial.
        no_m1 (int): Number of minus ones in the polynomial.

    Returns:
        bool: True if the polynomial is ternary with the given alpha and beta, False otherwise.
    """
    coeff_ones = coeff.count(-1) + coeff.count(1)
    return coeff_ones <= len(coeff) and coeff.count(1) == no_1 and coeff.count(-1) == no_m1

#==========================================================================================

def find_mod_inv(int_a,modulus):
    '''
    This function calculates the modular multiplicative inverse of an integer 'a' modulo 'm'.

    PARAMETERS:
        int_a (int): The integer for which the modular inverse is to be found.
        modulo (int): The modulus.

    RETURNS:
        int or None: The modular multiplicative inverse of 'int_a' modulo 'modulus', if it exists.
                     If no modular inverse exists (i.e., if 'int_a' and 'modulus' are not coprime),
                     None is returned.

    '''
    gcd_val, coeff_x, _ = extended_gcd(int_a, modulus)
    if gcd_val != 1:
        return None  # modular inverse does not exist
    return coeff_x % modulus

def mod_fraction(frac,modulus):
    '''
    This function calculates the modular inverse of a fraction 'f' modulo 'm'.

    PARAMETERS:
        f (Fraction): The fraction for which the modular inverse is to be found.
                      It should be provided as a Fraction object.
        m (int): The modulus.

    RETURNS:
        int: The modular inverse of the numerator of 'f' modulo 'm'.
    '''
    gcd_val,_,_ = extended_gcd(modulus, frac.denominator)
    if gcd_val != 1:
        return 0
    out = find_mod_inv(frac.denominator, modulus) * frac.numerator % modulus
    return out

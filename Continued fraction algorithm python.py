#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 20:22:15 2023

@author: Freja
"""

# This implementation uses the numpy and math modules.
import numpy as np
import math
import time

class GeneratorAE:
    '''
    An iterator class that generates a sequence of values A_i, E_i.

    Args:
        N (int): The input number used in the factorisation.
        k (int): A parameter used in the generation of values.

    Attributes:
        N (int): The input number.
        NK (int): The product of N and k.

    Methods:
        __init__(self, N, k): Initialises the GeneratorAE object with the given N and k values.
        
        __iter__(self): Resets the iterator and returns itself.
        
        __next__(self): Generates the next value in the sequence.
        
        AE(self): Returns the current value of A and its smallest remainder (E_i).
        
        smallestRemainder(self): Calculates and returns the smallest remainder of A**2 modulo N.
    '''

    def __init__(self, N, k):
        '''
        Initialises the GeneratorAE object with the given N and k values.
        Args:
            N (int): The input number used in the factorisation.
            k (int): A parameter used in the generation of values.
        '''
        self.N = N
        self.NK = N * k
        
    def __iter__(self):
        '''
        Resets the iterator and returns itself.
        Returns:
            GeneratorAE: The iterator object itself.
        '''
        self.a0 = math.isqrt(self.NK)
        self.a = 0
        self.m = 0
        self.d = 1
        self.Am1 = 1
        self.Am2 = 0
        self.A = 0
        return self
    
    def __next__(self):
        '''
        Generates the next value in the sequence using the algorithms described
            in my thesis.

        Returns:
            int: The next value in the sequence.
            
        Credit: Niels Lauritzen, who inspired me to combine the algorithm which calculates the
        terms of the continued fraction of the square root of k*N, and the algorithm which 
        calculates the A_i values from the terms in a single method. 
        '''
        self.a = (self.a0 + self.m) // self.d
        self.m = self.d * self.a - self.m
        self.d = (self.NK - self.m * self.m) // self.d
        self.A = (self.a * self.Am1 + self.Am2) % self.NK
        self.Am2 = self.Am1 % self.NK
        self.Am1 = self.A        
        return self
    
    def AE(self):
        '''
        Returns the current value of A and the smallest remainder of A**2 modulo N.

        Returns:
            tuple: A tuple containing the current value of A and the smallest remainder
                    of A**2 modulo N.
        '''
        return (self.A, self.smallestRemainder())

    def smallestRemainder(self):
        '''
        Calculates and returns the smallest remainder of A**2 modulo N.

        Returns:
            int: The smallest remainder of A**2 modulo N.
        '''
        E_pos = self.A * self.A % self.N 
        E_neg  = E_pos - self.N
        if abs(E_pos) < abs(E_neg):
            return E_pos
        else:
            return E_neg
        
        

class ContinuedFractionFactorisation:
    '''
    A class that implements the continued fraction factorisation algorithm
        to find factors of a given number N.

    Args:
        N (int): The number to be factored.
    Optional args:  
        max_runs (int) = The maximum number of times to add extra iterations. (Default = 20)
        max_time (int) = The maximum run time in seconds. (Default = 600)

    Attributes:
        N (int): The number to be factored.
        
        k (int): The value of k used in the factorisation process.
        
        upper_bound (int): The upper bound of the primes in the factor base. 
        
        iterations (int): The number of pairs (A_i,E_i) to be generated.
        
        factor_base (list): The factor base.
        
        len_primes (int): The length of the factor base.
        
        EAF_list (list): A list to store the values (E_i, A_i, vector of factors of E_i).
        
        max_time (int): The maximum allowed time for factorisation in seconds.
        
        start_time (float): The start time of the factorisation process.
        
        number_of_runs (int): The number of times the factorisation process has been restarted
                                with more iterations.
                                
        generator (GeneratorAE): An instance of the GeneratorAE class used for generating 
                                (A_i,E_i) values.
                                
        generator_iterations (iterator): An iterator object for the generator.

    Methods:
        __init__(self, N, max_runs = 20, max_time = 600): Constructor.
    
        primesLessThanP(self, p): Returns a list of primes less than or equal to p.
        
        legendre(self, k, p): Computes the Legendre symbol of k and p.
        
        kChoice(self): Chooses the appropriate value of k.
        
        upperBoundAndIterations(self): Determines the upper bound and number of iterations 
                                        based on the number of digits in N*k.
                                        
        factorBase(self): Generates the factor base.
        
        primeFactors(self, E): Computes a vector of the powers of the prime factors of a 
                                number e.
                                
        gaussEliminationMod2(self, M): Performs Gauss elimination modulo 2 on a matrix M.
        
        greatestCommonDivisor(self, n, m): Computes the greatest common divisor of n and m.
        
        checkSquares(self): Checks for non-trivial factors using the obtained square.
        
        generate(self): Generates A values and performs factorization.
        
        continueGeneration(self): Continues the generation process if the time limit 
                                        and max number of runs is not reached.
                                        
        timeUp(self): Checks if the elapsed time is more than max_time.
        
        printInfo(self): Prints N, k and initial number of iterations.
        
        smallFactor(self): Finds a prime factor of a small integer N < 10000.      
                                                              
        algorithm(self): Executes the continued fraction factorisation algorithm and returns the result.
    '''

    def __init__(self, N, max_runs = 20, max_time = 600):
        '''
        Initializes the ContinuedFractionFactorisation object with the given N value.

        Args:
            N (int): The number to be factored.
        '''
        self.N = N if N > 0 else -N
        self.k = self.kChoice()
        self.upper_bound, self.iterations = self.upperBoundAndIterations()
        self.factor_base = self.factorBase()
        self.EAF_list = []
        self.len_primes = len(self.factor_base)
        self.max_time = max_time
        self.max_runs = max_runs
        self.start_time = time.time()
        self.number_of_runs = 0
        self.generator = GeneratorAE(self.N, self.k)
        self.generator_iterations = iter(self.generator)
        
        
    def primesLessThanP(self, p):
        '''
        Returns a list of primes less than or equal to p.
        Args:
            p (int): The upper limit for generating primes.

        Returns:
            list: A list of primes less than or equal to p.
            
        Implemented with assistance from ChatGPT.
        '''
        sieve = [True] * (p + 1)
        sieve[0] = sieve[1] = False
        for i in range(2, int(p ** 0.5) + 1):
            if sieve[i]:
                sieve[i * i: p + 1: i] = [False] * len(sieve[i * i: p + 1: i])
        return [i for i in range(2, p + 1) if sieve[i]]

    def legendre(self, k, p):
        '''
        Computes the Legendre symbol of k and p.
        Args:
            k (int): The value of k.
            p (int): The prime.
        Returns:
            int: The Legendre symbol of k and p.
            
        Inspired by Stef on stackoverflow: 
        https://stackoverflow.com/questions/71020225/how-do-i-find-legendres-symbol
        '''
        if self.N * k % p == 0:
            return 0
        e = (p - 1) // 2
        r = pow(self.N * k, e, p) 
        return r - p if r > 1 else r

    def kChoice(self):
        '''
        Chooses the appropriate value of k.
        Returns:
            int: The chosen value of k.
        '''
        if self.N == 2 ** 128 + 1:
            return 257
        primes = self.primesLessThanP(31)    
        def k_len(k):
            quadr_primes = [i for i in primes if self.legendre(k, i) == 1 or self.legendre(k, i) == 0]
    
            if 2 in quadr_primes and 3 in quadr_primes:
                return len(quadr_primes)
            else: 
                return 0
        k_list = [k for k in self.primesLessThanP(300) if k_len(k) > 0]
        k_list = sorted(k_list, key=lambda x: k_len(x), reverse=True)
        return k_list[0]
    
    def upperBoundAndIterations(self):
        '''
        Determines the upper bound and number of iterations based on the number 
            of digits in N*k.
        Returns:
            tuple: A tuple containing the upper bound and the number of iterations.
        '''
        digits = len(str(self.N * self.k))
        primes_upper_bound = int(75 * math.exp(0.145 * digits)) 
        iterations = int(0.0011 * primes_upper_bound ** 2)
        
        return primes_upper_bound, iterations

    def factorBase(self):
        '''
        Generates the factor base by adding primes if N*k is a quadratic 
            residue modulo p, and adding the number -1.
        Returns:
            list: The factor base.
        '''
        primes = self.primesLessThanP(self.upper_bound)
        quadr_primes = [p for p in primes if self.legendre(self.k, p) == 1]
        return [-1] + quadr_primes
    
    def primeFactors(self, E):
        '''
        Computes a vector of the exponents of the the prime factors of 
            a given number E, only using the primes in the factor base.
        Args:
            E (int): The number to be factored.

        Returns:
            list: The exponents of the prime factors of E.

        '''
        factors_list = [0] * self.len_primes
        if E < 0:
            factors_list[0] += 1
            E = -E
        E_initial = E    
        E_sqrt = math.isqrt(E)
        for i, p in enumerate(self.factor_base):
            if i == 0:
                continue
            if p > E_sqrt:
                break
            if i >  self.len_primes // 2 and E == E_initial:
                break
            while E % p == 0:
                E = E // p
                factors_list[i] += 1
                if E == 1:
                    return factors_list
        return []
    
    def gaussEliminationMod2(self, M):
        '''
        Performs Gaussian elimination with modulo 2 on the given matrix.

        Args:
            M (numpy.ndarray): The matrix to perform Gaussian elimination on.

        Returns:
            list: The rows of the identity where the original matrix has a null 
                    row after Gaussian elimination.
            
        Credit: 
        ChatGPT has helped improve the efficiency of the function by introducing XOR and argmax.
        '''
        n, m = M.shape
        I = np.identity(n, dtype=M.dtype)
        M = np.concatenate((M, I), axis=1)
        for c in range(n):
            r = np.argmax(M[c:, c]) + c
            if M[r, c] == 0:
                continue
            M[[c, r]] = M[[r, c]]
            for r in range(c+1, n):
                if M[r, c] == 1:
                    M[r] ^= M[c]  # use XOR instead of addition
        null_rows = [row[m:] for row in M if not np.any(row[:m])]
        return null_rows

    def greatestCommonDivisor(self, n, m):
        '''
        Computes the greatest common divisor of two numbers.
        Args:
            n (int): The first number.
            m (int): The second number.

        Returns:
            int: The greatest common divisor of n and m.
        '''
        if n < 0:
            n = -n
        if m < 0:
            m = -m
        while m:
            n, m = m, n % m
        return n

    def checkSquares(self):
        '''
        Checks for non-trivial squares using the squares found by using the gaussEliminationMod2 
        method on the matrix of factorisations of E_n values.

        Returns:
            int: The non-trivial factor if found, otherwise 1.
        '''
        factors_list = [eaf[2] for eaf in self.EAF_list]
        A_n_list = np.array([eaf[1] for eaf in self.EAF_list], dtype=object)
        E_list = np.array([eaf[0] for eaf in self.EAF_list], dtype=object)
        A_mod2 = np.array(factors_list, dtype=np.int8) % 2
        zero_rows = list(self.gaussEliminationMod2(A_mod2))
        for row in zero_rows:
            prod_An = (np.prod(A_n_list[row == 1])) % self.N
            prod_E = (np.prod(E_list[row == 1]))
            sqrt_E = math.isqrt(prod_E)
            gcd_res = self.greatestCommonDivisor(self.N, (prod_An - sqrt_E) % self.N)
            if gcd_res != 1 and gcd_res != self.N:
                return gcd_res
        return 1

    def generate(self):
        '''
        Generates A_i, E_i values and tries to find a factor using the checkSquares method.

        Returns:
            int: The non-trivial factor if found, otherwise 1.
        '''
        for i in range(self.iterations):
            if self.timeUp():
                break
            next(self.generator_iterations)
            A, E = self.generator.AE()
            factors = self.primeFactors(E)
            if len(factors) > 0:
                self.EAF_list.append((E, A, factors))
                
        if len(self.EAF_list) > 0:
            res = self.checkSquares()
            if res > 1:
                return res
        return self.continueGeneration()

    def continueGeneration(self):
        '''
        Continues the generation process if the time limit and max number of runs is not 
        reached.

        Returns:
            int: The non-trivial factor if found, otherwise 1.
        '''
        if not self.timeUp() and self.number_of_runs < self.max_runs:
            self.number_of_runs += 1
            print("runs =", self.number_of_runs)
            self.iterations = 2*self.upper_bound
            return self.generate()
        else:
            if self.number_of_runs == self.max_runs:
                print("No non-trivial factor was found within a reasonable number of iterations.")
                print("N could be a prime number or a square.")
            else:
                print("No non-trivial factor was found within the time limit of %s seconds." % self.max_time)
            return 1

    def timeUp(self):
        '''
        Returns a boolean, true if the elapsed time is more than max, else false.
        '''
        elapsed_time = time.time() - self.start_time
        return elapsed_time > self.max_time 
    
    def printInfo(self):
        '''
        Prints N, k and initial number of iterations.
        '''
        print("N =", self.N)
        print("k =", self.k)
        print("Initial number of values of A_i and E_i =", self.iterations)
    
    def smallFactor(self):
        '''
        Returns a prime factor or None.
    
        Returns:
            p (int): prime factor of N.
        '''
        primes = self.primesLessThanP(math.isqrt(self.N))
        for p in primes:
            if N % p == 0:
                return p
        
    def algorithm(self):
        '''
        Executes the continued fraction factorisation algorithm.
        Checks if N is too small, then uses brute force to find a factor. 

        Returns:
            res (int): The factor which was found, either a non-trivial factor or 1.
        '''
        if self.N < 10000: 
            p = self.smallFactor()
            if p:
                return p
            print("N is a prime number.") 
            return 1
            
        self.printInfo()
    
        res = self.generate()
        return res



start_time = time.time()

# Try out the algorithm with different integers 

N = 2**128+1 # F_7
N = 2**64+1  # F_6
N = 12007001 

factor = ContinuedFractionFactorisation(N)  # Optional arguments are max_runs and max_time. 
res = factor.algorithm()

if res > 1:
    print("\nA factor of N is: %s \n" %res)
end_time = time.time()

print("Execution time:", round(end_time - start_time,4), "seconds")




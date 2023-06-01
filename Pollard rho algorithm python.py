#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 20:24:34 2023

@author: Freja
"""

import random
import math

def greatestCommonDivisor(n, m):
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

def primesLessThanP(p):
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


def smallNumber(N):
    '''
    Returns a prime factor or None.
    
    Returns:
        p (int): prime factor of N.
    '''
    for p in primesLessThanP(math.isqrt(N)):
        if (N % p) == 0:
            return p

def pollard(N):
    '''
    Returns a the factor of N.
            Parameters:
                N (int): The integer to factor
            Returns: 
                p (int): A factor of N
    '''
    if N < 1:
        N *= -1
    if N < 100:
        if smallNumber(N) != None:
            return smallNumber(N)
    
    # First, we define a recursive function with a random constant 
    a = random.randrange(1, N-1)  
    def f(v):
        return (v**2 - a) % N
    # Random initial value of x 
    x = random.randrange(2, N-1)
    y = x
    Q = 1
    # Choose how often we check the gcd
    loops = 10 if N < 10**20 else 100
    if N < 10**4:
        loops = 1
    # We try to find a factor of N within 10**12 loops 
    for i in range(10**12):
        # Update x and y such that y = x_2i
        x = f(x)
        y = f(y)
        y = f(y) 
        # Calculate the product of x_2i - x_i
        Q = (Q*(y-x)) % N 
        # For every 1 or 10 or 100 loops, we check if we have a factor of N
        if i % loops == 0:
            p = greatestCommonDivisor(Q,N)
            if p>1 and p<N:
                return p

import time

# Feel free to test the algorithm with different integers N
N = 2**64 + 1 
N = 2**128 + 1
N = 2**256 + 1
N = -1234567890987654321

start_time = time.time()
print("A factor of %s is %s. \n" % (N,pollard(N)))
end_time = time.time()

print("Execution time:", round(end_time - start_time,4), "seconds")
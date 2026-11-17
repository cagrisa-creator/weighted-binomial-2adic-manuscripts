"""Verificações diretas das novas fórmulas. Executar com Python 3."""
from fractions import Fraction
from math import comb
from pathlib import Path
import json

def valuation(n):
    return (n & -n).bit_length()-1 if n else float('inf')

def mahler(c, x, power, precision):
    modulus = 1 << precision
    last = (precision-1)//valuation(c)
    choose = 1
    total = 1
    for k in range(1,last+1):
        choose = choose*(x-k+1)//k
        total += pow(c,k,modulus)*pow(choose,power,modulus)
    return total%modulus

def reduced(c,z,n):
    # Evaluate at a congruent nonzero representative before dividing by z^3.
    z = z or (1<<n)
    e = valuation(z)
    top = mahler(c,2*z,3,n+7+3*e)-1
    denominator = 1<<(7+3*e)
    assert top%denominator == 0
    return (top//denominator)*pow(z>>e,-3,1<<n)%(1<<n)

def kernel(u,n):
    m=1<<n
    val=0
    for k in range(1,n+12):
        exponent=3*(k-1)-3*valuation(k)
        odd=k>>valuation(k)
        val+=(-1)**(k-1)*(1<<exponent)*pow(odd,-3,m)*pow(u,k-1,m)
    # omitted k>=n+12: 3(k-1)-3 floor(log2(k)) >= n
    return val%m

def lift(func,n):
    a=0
    for i in range(n):
        if func(a,i+1):
            a+=1<<i
        assert func(a,i+1)==0
    return a

def cubic():
    root=8*lift(kernel,48)
    assert root%64==8
    count=0
    for u in range(1,128,2):
        c=8*u
        a=lift(lambda z,n:reduced(c,z,n),32)
        assert valuation(a)==valuation(c-root)-4
        assert mahler(c,2*a,3,39)==1
        # Direct permutation test for the entire normalized map modulo 256.
        assert len({reduced(c,z,8) for z in range(256)})==256
        count+=1
    return {'unit_classes_tested':count,'root_precision_bits':32,
            'lambda_mod_2_51':root,'normalized_permutation_modulus':256}

def quadratic():
    count=0
    for d in [5,7,9,11]:
        j=(d-3)//2
        for unit in range(1,16,2):
            c=2+(1<<d)*unit
            def g(w):
                f=mahler(c,-1-4*(1<<j)*w,2,12+3*j)
                numerator=((1-c)*f-1)%(1<<(12+3*j))
                assert numerator%(1<<(6+3*j))==0
                return numerator>>(6+3*j)
            for a in [1,3]:
                assert ((g(a+8)-g(a))%64)==8*((1-unit)%8)
                count+=1
    return {'derivative_residues_checked':count}

if __name__=='__main__':
    report={'cubic':cubic(),'quadratic':quadratic()}
    print(report)
    Path(__file__).with_name('verificacao.json').write_text(json.dumps(report,indent=2))

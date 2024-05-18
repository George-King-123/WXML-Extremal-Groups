from math import factorial, log
from shared_code import abelian

def compare_bound():
  r = range (100, 1000, 100)
  for n in r:
    for k in r:
      for d in range(2, 6):
        print(f"n = {n}, k = {k}, d = {d}")
        bound = (factorial(d))**n * abelian(n, k)
        k_n = k**n 
        if k_n < 1000* bound:
          print("bad")


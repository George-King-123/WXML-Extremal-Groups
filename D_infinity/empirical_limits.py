# run with python -m D_infinity.empirical_limits

from math import factorial, sqrt
from shared_code import choose
from D_infinity.compute_s_n import FAST_FORMULA as formula

def gamma_k_n_growth(small_k, large_n, p): 
  N = large_n 
  k = small_k

  s_n = formula(N, k, p)
  empirical_limit = s_n / (N**(k-1))
  analytic_limit = 2**(k - 2 * p + 1) * choose(2 * p - 2, p- 1) / factorial(k - 1)
  
  print(f"{analytic_limit=}")
  print(f"{empirical_limit=}")

def gamma_n_n_growth(n, p): 
  val = formula(n, n, p=p)
  empirical_limit = val ** (1/n)
  print(f"{empirical_limit=}")

  analytic_limit = 3 + 2 * sqrt(2)
  print(f"{analytic_limit=}")

def main():
  # gamma_n_n_growth(n=200, p=2)
  # gamma_k_n_growth(small_k=4, large_n=1000, p = 1)
  pass

if __name__ == "__main__":
  main()
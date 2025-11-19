from math import factorial, sqrt
from shared_code import incl_range, format_large_num, choose
from compute_s_n import old_formula, new_formula

def gamma_D_inf(n, k, noisy = False):
  biggest = -1
  achieved_at = 0
  for t in incl_range(0, k):
    val = old_formula(n, k, t)
    if val > biggest:
      biggest = val
      achieved_at = t

  if noisy:
    print(f"γ({n=}, {k=}) = {format_large_num(biggest)}, achieved with {achieved_at} signs")

  return biggest 


def growth_rate_one_sign_assumption(small_k, large_n):
  val = old_formula(n = large_n, k = small_k, t = 1)
  theta_bound = large_n**(small_k - 1)
  print(f"Actual ratio: {val/theta_bound}")

  coefficient = (2**(small_k-1))/factorial(small_k-1)
  print(f"Conjectured ratio: {coefficient}")

def gamma_n_n_exponent_two_sign_assumption(n):
  val = old_formula(n, n, t=2)
  experimental = val ** (1/n)
  print("experimental: {}".format(experimental))

  conjecture = 3 + 2 * sqrt(2)
  print("conjecture: {}".format(conjecture))

def check_new_form_matches_old(N = 30, k = 30, noisy = True): 
  for n in range(1, N + 1):
    for k in range(1, n + 1):
      for t in range(1, k+1): 
        new = new_formula(n, k, t)
        old = old_formula(n, k, t) 
        if noisy: 
          print(f"{new=} {old=}")
        assert new == old 

def check_limits_with_dif_t_vals(): 
  t = 14
  k = 14
  N = 1000
  s_n = old_formula(N, k, t)
  quotient = s_n / (N**(k-1))
  conj_limit = 2**(k - 2 * t + 1) * choose(2 * t - 2, t- 1) / factorial(k - 1)
  
  print(f"{conj_limit = }")
  print(f"{quotient = }")

def gamma_n_n_limit():
  N = 400
  alpha = .7
  t = round(N * alpha)
  s_n = old_formula(N, N, t)
  conj_limit = 3 + 2 * sqrt(2)  
  emprirical_limit = s_n**(1/N) 
  print(f"{conj_limit = }")
  print(f"{emprirical_limit = }")

def main():
  gamma_D_inf(n=20, k = 7, noisy=True)

if __name__ == "__main__":
  main()
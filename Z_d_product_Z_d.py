# run with Python -m Z_d_product_Z_d

from shared_code import get_sparse_d_tuple, compute_Sn, incl_range
from group_operations import op_Z_d_Z_d
from math import comb, ceil, floor, factorial

# generates a set S \subset Z^d product Z_d with |S| = k,
# where the number of elements with i in the second component is distro_of_signs[i]
def make_S(k, d, distro_of_signs):
  if len(distro_of_signs) != d:
    raise ValueError("must say how many of each sign")
  
  if sum(distro_of_signs) != k:
    raise ValueError("total number of elements does not agree with the number of signs")
    
  S = set()
  for i in range(0, d):
    # add the right amount of elements of that sign
    for _ in range(0, distro_of_signs[i]):
      S.add((get_sparse_d_tuple(d), i))
  return S

# Returns |S^n| when S is constructed with the sign distrubtion as described
# in make_S
def compute_s_n_with_sign_distribution(n, d, distribution):
  S = make_S(d = d, k = sum(distribution), distro_of_signs= distribution)
  s_n = compute_Sn(S, n, op_Z_d_Z_d)
  return len(s_n)

def formula_when_all_zero_one_in_snd_comp(n, k, p, d):
  if p == 0:
    return comb(n + k -1, k - 1)
  
  if p == k:
    res = comb(ceil(n/d) + k - 1, k - 1) ** (n % d)
    res *= comb(floor(n/d) + k - 1, k - 1)**(d - (n % d))
    return res 

  res = comb(n + (k - p) - 1, k - p - 1)
  for i in incl_range(1, n):
    cur_term = comb(ceil(i/d) + p - 1, p-1)**(i % d) 
    cur_term *= comb(floor(i/d) + p - 1, p-1)**(d - (i % d))
    num_forms = min(d, i+1)
    cur_term *= comb(n - i + num_forms * (k - p) - 1, num_forms * (k - p) - 1)
    res += cur_term
  return res


# test that the formula for just 0s and 1s in the second component
# works for all values of d,n,k up to the given upper bounds
# these have to be small because we explicitly compute |S^n| in order 
# to check the formula 
# ex: test_formula_all_one_many_vals(5, 5, 5)
def test_formula_all_one_many_vals(d_upper, n_upper, k_upper):
  print("testing that the Z_d formula works")
  def formula_works(n, k, p, d):
    S = make_S(d = d, k = k, distro_of_signs= [k - p, p] + ([0] * (d-2)))
    s_n = compute_Sn(S, n, op_Z_d_Z_d)
    return len(s_n) == formula_when_all_zero_one_in_snd_comp(n, k, p, d)
  

  for d in incl_range(2, d_upper):
    for n in incl_range (1, n_upper):
      for k in incl_range(1, k_upper):
        for p in incl_range(0, k):
          if not formula_works(n, k, p, d):
            print(f"FAIL {n=} {k=} {p=} {d=}")
            return
          
  print("success")


# test \lim_{n \to \infty} formula(t, d, k, n) / n^{d(k-1)}. 
# This didn't end up being important in the paper
# This call doesn't take too long: test_predicted_limit(big_number= 10 ** 6, k=10, p=4, d=7)
def test_predicted_limit(big_number, k, p, d):
    def analytic_limit(k, p, d):
      numerator = factorial(d * (p - 1))
      denominator = ((factorial(p - 1)) ** d) * (d ** (d * (p - 1))) * factorial(d * (k - 1))
      return numerator / denominator

    empirical_limit = formula_when_all_zero_one_in_snd_comp(n=big_number, k=k, p=p, d=d) / big_number ** (d * (k - 1))
    print(f"{empirical_limit=}")
    
    analytic_limit_value = analytic_limit(k, p, d)
    print(f"{analytic_limit_value=}")

def find_maximizing_t_val(n, k, d): 
    maximizers = []
    maximimum = -1 

    for p in incl_range(0, k): 
        cur_val = formula_when_all_zero_one_in_snd_comp(n=n, k=k, p=p, d=d)
        if cur_val > maximimum:
            maximizers = [p]
            maximimum = cur_val 
        elif cur_val == maximimum:
            maximizers.append(p) 
    
    return maximizers

def main():
  # test_formula_all_one_many_vals(6, 6, 5)
  print("Not doing anything, look at main() for options of what this script can do")


if __name__ == "__main__":
  main()
from shared_code import get_sparse_d_tuple, compute_Sn, incl_range
from group_operations import op_Z_d_Z_d
from math import comb, ceil, floor, factorial

# generates a set S \subset Z^d product Z_d 
# with |S| = k,
# where the number of elements with i in the 
# second component is distro_of_signs[i]
# thus, k = |S| = sum(distro_of_signs)
#
# ex: make_S(3, 6, (1, 3, 2)) will make a set 
# S = {(x1, 0), (x2, 1), (x3, 1), (x4, 1), (x5, 2), (x6, 2)}
# where the x_j are 3-tuples of sparse integers 
# notice there is 1 element with a 0 sign, 3 elements with a 1 sign
# and 2 elements with a 2 sign
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

def show_we_only_need_one(n, k, d):
  pass

# formula for |S^n| when |S| = k, and there are t elements
# with a 1 in the 2nd component, k-t elements with a 0
# (i.e. no other signs)
# We can prove this formula and it matches with all simulations, 
# (see the test code)
def formula_when_all_one(n, k, t, d):
  # this is just the abelian case
  if t == 0:
    return comb(n + k -1, k - 1)
  
  if t == k:
    res = comb(ceil(n/d) + k - 1, k - 1) ** (n % d)
    res *= comb(floor(n/d) + k - 1, k - 1)**(d - (n % d))
    return res 

  res = comb(n + (k - t) - 1, k - t - 1)
  for i in range(1, n+1):
    cur_term = comb(ceil(i/d) + t - 1, t-1)**(i % d) 
    cur_term *= comb(floor(i/d) + t - 1, t-1)**(d - (i % d))
    num_forms = min(d, i+1)
    cur_term *= comb(n - i + num_forms * (k - t) - 1, num_forms * (k - t) - 1)
    res += cur_term
  return res

# We think the largest value of |S^n| is always achieved with just 0 and 1 
# in the second component (we DONT have a proof of this). 
# Thus we only need to check how many 1s to use, 
# i.e. try all possibly values of t from 0 to k in the above formula
def conjectured_gamma_k_n(k, n, d, noisy = False):
  best = -1
  achieved_at = -1
  for t in range(0, k+1):
    val = formula_when_all_one(n, k, t, d)
    if val > best:
      best = val
      achieved_at = t 
  
  if noisy:
    print("k = {}, n = {}, best at {} signs".format(k, n, achieved_at))
    
  return best


# Conjecture: when k is fixed, n gets large we get an exponent of n^{d(k-1)}
# with a coefficient of 1/(d(k-1)!)
# I think we had a start of a proof with this. It matches with the simulation below,
# n has to be quite large for it to converge. 
# e.x. try exponent(1000, 5, 3) and exponent(10000, 5, 3)
def exponent(n, k, d):
  print("\'Actual\'")
  size_s_n = conjectured_gamma_k_n(n=n, k=k, d=d)
  print(size_s_n/(n**(d * (k - 1))))
  
  print("Conjecture:")
  print(1/(factorial(d*(k - 1))))



# we are interested in \lim_{n \to \infty} formula(t, d, k, n) / n^{d(k-1)} 
def predicted_limit(k, t, d):
    numerator = factorial(d * (t - 1))
    denominator = ((factorial(t - 1)) ** d) * (d ** (d * (t - 1))) * factorial(d * (k - 1))
    return numerator / denominator

def test_predicted_limit(big_number, k, t, d):
    numerical_approx_of_limit = formula_when_all_one(n=big_number, k=k, t=t, d=d) / big_number ** (d * (k - 1))
    print(f"Numerical approximation of limit: {numerical_approx_of_limit}")
    predicted = predicted_limit(k, t, d)
    print(f"Predicted limit: {predicted}")

    print(f"Difference: {abs(numerical_approx_of_limit - predicted)}")

def find_best_t_val_limit(k, d):
    return max(range(1, k+1), key=lambda t: predicted_limit(k, t, d))

def main():
  # test_predicted_limit(big_number= 10 ** 6, k=10, t=4, d=7)
  print(find_best_t_val_limit(k=20, d=2))

def find_maximizing_t_val(n, k, d): 
    maximizers = []
    maximimum = -1 

    for t in incl_range(0, k): 
        cur_val = formula_when_all_one(n=n, k=k, t=t, d=d)
        if cur_val > maximimum:
            maximizers = [t]
            maximimum = cur_val 
        elif cur_val == maximimum:
            maximizers.append(t) 
    
    return maximizers

if __name__ == "__main__":
  main()
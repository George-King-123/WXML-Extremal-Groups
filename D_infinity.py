from math import floor, ceil, factorial, sqrt, comb as comb_original


# compute gamma(n, k) for D_infinity
# we have an explicit formula given n,k, and the number of signs,
# but we don't know what number of signs is best. So this function 
# has to check all possible number of signs between 0 and k to live up 
# to its name. Of course we use the formula for this
def gamma_D_inf(n, k, noisy = False):
  biggest = -1
  achieved_at = 0
  for t in range(0, k + 1):
    val = compute_s_n_with_formula_OLD_FORM(n, k, t)
    if val > biggest:
      biggest = val
      achieved_at = t

  if noisy:
    print("gamma(n = {}, k = {}) is {}, achieved with {} signs".format(n, k, biggest, achieved_at))

  return biggest 


def compute_s_n_new_formula(n, k, t):
  def comb(n, k):
    if (k > n):
      return 0
    if (k < 0 or n < 0):
      return 0
    if (k == 0):
      return 1
    return comb_original(n, k)
  
  def rangeIncl(start, end):
    return range(start, end + 1)

  if t == 0: 
    return comb(n + k - 1, k - 1)
  elif t == 1: 
    res = 0 
    for q in range(1, min(k-1, n-1) + 1):
      res += pow(2, q) * comb(k-1, q) * comb(n-1, q)
    return res + comb(n + k - 2, k - 2) + 1
  elif t < k:
    res = 0 
    for p in rangeIncl(1, min(t-1, n)):
      for ell in rangeIncl(p, n):
        for q in rangeIncl(0, min(k-t, n-ell)):
          for i in rangeIncl(q, n-ell):
            if (n - ell - i) % 2 != 0:
              continue 
            res += comb(t, p) * comb(ceil(ell/2)-1, p-1) * comb(floor(ell/2) + t - p - 1, t-p-1) * (2**q) * comb(k-t, q) * comb(i-1, q-1)
 
    # + .... -
    for q in rangeIncl(1, min(k-t, n-2)):
      for i in rangeIncl(q-1, n-2):
        if (n-i) % 2 != 0: 
          continue
        res += (2**q) * comb(k-t, q) * comb(i-1, q-1) 

    return res + comb(n + k - t - 1, k - t - 1)
  elif t == k: 
    res = 0 
    for p in rangeIncl(1, min(t-1, n)):
      for c in rangeIncl(0, floor(n-p/2)):
        res += comb(t, p) * comb(ceil(n/2 - c) - 1, p-1) * comb(floor(n/2 -c) + t-p-1, t-p-1)
    return res + (1 if n % 2 == 0 else 0) 

  else:
    raise Exception("t > k")


# computes |S^n| when |S| = k and there are t negative signs in S
# works for k >= 1, n >= 1
def compute_s_n_with_formula_OLD_FORM(n, k, t, noisy = False):
  def comb(n, k):
    if (k > n):
      return 0
    if (k < 0 or n < 0):
      return 0
    if (k == 0):
      return 1
    return comb_original(n, k)

  def I_0(t, k, n):
    return comb(n + (k - t) - 1, (k - t) - 1)

  # returns I_1(t, k, n)
  def I_1(t, k, n):
    res = 0
    for q in range (1, min(k-t, n-1) + 1):
      res += pow(2, q)*comb(k-t, q)*comb(n-2, q-1)
    res *= t
    return res

  # returns I_2(t, k, n)
  def I_2(t, k, n):
    res = 0
    for q in range (1, min(k-t, n-2) + 1):
      res += (pow(2, q) - 1) * comb(k-t, q) * comb(n-3, q-1) 
    # res *= t
    return res

  # returns R(t, k, n)
  def R(t, k, n): 
    if t in {0, 1}: 
      return 0

    res = 0
    for i in range(2, n + 1):
      omega_value = 0
      sigma_value = 0
      for q in range (1, min(k - t, n - i) + 1):
        omega_value += (pow(2, q)) * comb(k-t, q) * comb(n-i -1, q-1)

      for m in range (1, min(t-1, ceil(i / 2)) + 1):
        sigma_value += comb(t, m) * comb(ceil(i / 2) - 1, m - 1) * comb((t-m) + floor(i / 2) - 1, (t - m) - 1)
      
      if i == n: # note: when i = n, want the first summation to be 1
        omega_value = 1
      res += (omega_value * sigma_value)
    return res

  if t == 0:
    return comb(n + k - 1, k - 1)
  
  # set up for the even or odd case
  start = 1 if n % 2 == 0 else k
  cur = 2 if n % 2 == 0 else 3

  tot_I_0_contrib = 0 
  tot_I_1_contrib = 0
  tot_I_2_contrib = 0
  tot_R_contrib = 0

  while cur <= n: 
    tot_I_0_contrib += I_0(t, k, cur)
    tot_I_1_contrib += I_1(t, k, cur)
    tot_I_2_contrib += I_2(t, k, cur)
    tot_R_contrib += R(t, k, cur)
    cur += 2

  total = tot_I_0_contrib + tot_I_1_contrib + tot_I_2_contrib + tot_R_contrib + start
  if noisy:
    print(f"I_0: {round(tot_I_0_contrib/total, 4)}, I_1: {round(tot_I_1_contrib/total, 4)}, I_2: {round(tot_I_2_contrib/total, 4)}, R: {round(tot_R_contrib/total, 4)}")
  
  return total
  


# We think that for small k and large n that gamma(k, n) grows like n^{k-1}
# with a coefficient of 2^{k-1}/(k-1)!
# I believe we have the start of a proof of this? And the simulation below confirms
# the conjecture
# Note that when n >> k we are sure (but do not have a proof) that one sign is best,
# so I only check one sign
# ex: growth_rate_one_sign_assumption(50, 10000)
def growth_rate_one_sign_assumption(small_k, large_n):
  val = compute_s_n_with_formula_OLD_FORM(n = large_n, k = small_k, t = 1)
  theta_bound = large_n**(small_k - 1)
  print("Actual ratio: {}".format(val/theta_bound))

  coefficient = (2**(small_k-1))/factorial(small_k-1)
  print("Conjectured ratio: {}".format(coefficient))


# We think that gamma_n_n to the 1/n goes to 3 + 2 sqrt(3), check below
# Use two signs because that seems best for n = k, though we don't have a proof 
# and also aren't sure two signs always will win. But for the tractable values, up to 
# 700 or so, it does.
# ex: gamma_n_n_exponent_two_sign_assumption(200)
# see the old files for tests of two vs three signs and other nonsense functions
def gamma_n_n_exponent_two_sign_assumption(n):
  val = compute_s_n_with_formula_OLD_FORM(n, n, 2)
  experimental = val ** (1/n)
  print("experimental: {}".format(experimental))

  conjecture = 3 + 2 * sqrt(2)
  print("conjecture: {}".format(conjecture))

# compute partial sums of the generating function up to the given upper bounds, 
# evaluating at the given x and y values
def check_generating_function(nupper, kupper, xval, yval):
  def sum_terms():
    sum = 0
    for k in range(2, kupper + 1):
      for n in range(0, nupper + 1):
        sum += gamma_D_inf(n = n, k = k) * (xval)**n * (yval)**k
    return sum

  # the generating function for the abelian case is 1/(1 - x - y)
  def A_0(x, y):
    return y**2/(1-x-y)
  
  def A_1(x, y):
    return x * y - (y**2)/((1 - x - y - x * y) * (1 - x))
  
  def A_2(x, y):
    return x * A_1(x, y)
  
  def g(x, y):
    return A_0(x, y) + A_1(x, y) + A_2(x, y)
  
  term_sum = sum_terms()
  direct_eval = g(xval, yval)

  print("Summing terms gets {}".format(term_sum))
  print("Evaluating generating function gets {}".format(direct_eval))



# this is for a gut check that the function above is a reasonable way to check generating functions
# of course, xval and yval need to be small so that it converges
# ex: check_abelian_generating_function(30, 30, .1, .1)
def check_abelian_generating_function(nupper, kupper, xval, yval):
  def sum_terms():
    sum = 0
    for k in range(1, kupper + 1):
      for n in range(0, nupper + 1):
        sum += comb_original(n + k - 1, k - 1) * (xval)**n * (yval)**k
    return sum 
  
  # the generating function for the abelian case is 1/(1 - x - y)
  def g(x, y):
    return 1/(1-x-y)
  
  term_sum = sum_terms()
  direct_eval = g(xval, yval)

  print("Summing terms gets {}".format(term_sum))
  print("Evaluating 1/(1-x-y) gets {}".format(direct_eval))


def check_limits_with_dif_t_vals(): 
  t = 14
  k = 14
  N = 1000
  s_n = compute_s_n_with_formula_OLD_FORM(N, k, t)
  quotient = s_n / (N**(k-1))
  conj_limit = 2**(k - 2 * t + 1) * comb_original(2 * t - 2, t- 1) / factorial(k - 1)
  
  print(f"{conj_limit = }")
  print(f"{quotient = }")

def gamma_n_n_limit():
  N = 400
  alpha = .7
  t = round(N * alpha)
  s_n = compute_s_n_with_formula_OLD_FORM(N, N, t)
  conj_limit = 3 + 2 * sqrt(2)  
  emprirical_limit = s_n**(1/N) 
  print(f"{conj_limit = }")
  print(f"{emprirical_limit = }")

def main():
  compute_s_n_with_formula_OLD_FORM(2000, 10, 2, True)

if __name__ == "__main__":
  main()
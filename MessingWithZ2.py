from math import comb 

# change line 91 and mess around


# Let G = Z_2 \rtimes Z^2 

# compute the group operation in G given two ordered pairs g and h
# the second component is \pm 1
def op (g, h):
  if g[1] == 1:
    return ((g[0][0] + h[0][0], g[0][1] + h[0][1]), g[1]*h[1])
  elif g[1] == -1:
    return ((g[0][0] + h[0][1], g[0][1] + h[0][0]), g[1]*h[1])
  else:
    raise ValueError("bad 2nd component")


# given a set S of ordered pairs representing elements of D_infinity and an integer
# n >= 1, compute the set of all possible length n products of elements of S 
def compute_Sn(S, n):
  if n < 1:
    raise ValueError("should not give compute_Sn any value < 1")
  if n == 1:
    return S
  
  prod_with_one_less = compute_Sn(S, n - 1)
  s_n = set()
  for elt in S:
    for prod in prod_with_one_less:
      s_n.add(op(elt, prod))
  
  return s_n

# make a set S of size k, needs n to make integers sufficiently sparse
def make_S(n, k, num_neg):
  # give a method for assigning the sign of the ith element
  # e.g. return (-1)**i for alternating signs, return -1 for all negative
  def f(i):
    if i <= num_neg:
      return -1
    else:
      return 1

  S = {(((n+1)**i,(10*n+1)**i), f(i)) for i in range (1, k + 1)}
  return S

def run_sim(n, k):
  print("n = {} and k = {}".format(n, k))
  max_size = 0
  num_signs_max = 0
  for i in range(0, k+1):
    s_n = compute_Sn(make_S(n,k, i), n)
    cardinality_s_n = len(s_n)
    if cardinality_s_n > max_size:
      num_signs_max = i
      max_size = cardinality_s_n
    # Uncomment this to get actual sizes
    print("|S^n| = {} when we have {} negative signs".format(cardinality_s_n, i))
    if 0 < i and i < k:
      print("conjecture: {}".format(conj(n, k, i)))
  print("We achieved max size at {} signs".format(num_signs_max))

# formula that someone put in the overleaf. t is the number of negative signs
# this formula is correct, it matches the output of run_sim. So let's use
# get_values_from_conjecture
def conj(n, k, t):
  res = comb(n + (k - t) - 1, k - t- 1)
  for i in range(1, n+1):
    cur_term = comb(i - i//2 + t - 1, t-1) * comb(i//2 + t - 1, t-1)
    j_sum = 0
    for j in range(0, n-i+1):
      j_sum += comb(j + (k - t) - 1, (k - t) - 1) * comb(n -i - j + (k - t) -1, (k - t) - 1)
    cur_term *= j_sum
    res += cur_term
  return res

def get_values_from_conjecture(n, k):
  max_val = 0
  max_val_achieved_with_signs = 0
  for i in range(1, k):
    c = conj(n, k, i)
    if c > max_val:
      max_val = c
      max_val_achieved_with_signs = i
    #print("When we have {} signs we get |S^n| = {}".format(i, conj(n, k, i)))
  print("n = {}, k = {}, best with {} signs".format(n, k, max_val_achieved_with_signs))


#run_sim(n, k)
get_values_from_conjecture(2000, 10)

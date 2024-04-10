from math import comb 

# compute the group operation in D_infinity given two ordered pairs elt1 and elt2
def op (elt1, elt2):
  return (elt1[0] + elt1[1] * elt2[0], elt1[1] * elt2[1])

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
def make_S(n, k):
  #give a method for assigning the sign of the ith element
  # e.g. return (-1)**i for alternating signs, return -1 for all negative
  def f(i):
    if i <= 2:
      return -1
    else:
      return 1
  
  S = {((3*n+1)**i, f(i)) for i in range(1, k+1)}
  return S


def make_full_table(n_upper_bound, k_upper_bound):
  # results[k][n] will give will give |S^n| when |S| = k 
  results = [[0]*(n_upper_bound+1) for i in range(0, k_upper_bound + 1)]

  for n in range (1, n_upper_bound + 1):
    for k in range (1, k_upper_bound + 1):
      S = make_S(n, k)
      l = len(compute_Sn(S, n))
      results[k][n] = l

  # print everything out, just copy to Latex
  def print_results(results):
    print("\\begin{array}{" + "c|" * n_upper_bound + "}")

    for j in range (1, n_upper_bound):
      print(j, end = " & ")

    print(n_upper_bound, end = " \\\\")
    print("\n\hline")


    for j in range(1, k_upper_bound + 1):
      arr = results[j]
      l = len(arr)
      for i in range(1, l - 1):
        print(arr[i], end = " & ")
      print(arr[l - 1], end = " \\\\")
      print()
    print("\\end{array}")
  print_results(results)

# closed form for |S^n| when |S| = k and S has one element with a minus sign
# note: not correct
def conj_one_sign(k, n):
  def g(n, k, j):
    return comb(j - 1 + k - 1, k-1) * comb(n - j - 1 + k- 1, k - 1)
  res = comb(n + k - 2, k - 2)
  for j in range (1, n+1):
    res += g(n, k, j)

  if (k % 2 == 1):
    res -= g(n, k, (k+1)//2)
    c = (k+1)//2 - 1
    res += comb(c + k - 1, k-1) * (comb(c + k - 1, k-1) - 1) + 1
  return res

def make_nk_table(max):
  results = [0] * (max + 1)
  for i in range (1, max + 1):
    S = make_S(i, i)
    l = len(compute_Sn(S, i))
    results[i] = l


  z_vals = [0] + [comb(2 * n - 1, n - 1) for n in range (1, max + 1)]
  ratios = [0] + [results[i]/z_vals[i] for i in range(1, max + 1)]

  def print_results():
    print("\\begin{array}{c|c|c}")
    print("n & |S^n| & |S^n|/\\gamma_\mathbb{Z}(n) \\\\")
    print("\\hline")
    
    for i in range (1, max + 1):
      print("{} & {} & {} \\\\".format(i, results[i], round(ratios[i], 2)))

    print("\\end{array}")
  print_results()


make_full_table(10, 5)




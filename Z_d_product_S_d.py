# run with python -m Z_d_wreath_product_S_d

from math import factorial
from group_operations import op_wreath_prod
from shared_code import compute_Sn, get_sparse_d_tuple, WC, format_large_num, incl_range
import Z_d_product_Z_d
from tqdm import tqdm as loading_bar
import itertools

# returns a set of all tuples of length d having elements in
# {0, 1, ..., d-1} representing bijections
def make_s_d(d: int):
  res = set(itertools.permutations(range(d), d))
  assert len(res) == factorial(d)
  return res 

def perm_to_cycle_name(perm_tuple) -> str:
    d = len(perm_tuple)
    seen = [False] * d
    cycles = []

    for i in range(d):
        if seen[i]:
           continue 
        
        j = i
        cycle = []
        while not seen[j]:
          seen[j] = True
          cycle.append(j)
          j = perm_tuple[j]
        if len(cycle) > 1:
          cycles.append(cycle)

    if not cycles:
        return "e"

    cycle_strs = ["(" + " ".join(str(x) for x in cycle) + ")" for cycle in cycles]

    return "".join(cycle_strs)

def perms_to_str(perms):
   return "; ".join(perm_to_cycle_name(perm) for perm in perms)

# very, very slow. We try all possible second components, so this is the true value of gamma
def compute_gamma_k_n(k, n, d, noisy=False): 
  print(f"Computing γ({n=}, {k=}) for Z_d \\rtimes S_d, {d=}")
  Sd = make_s_d(d)
  all_possible_snd_component_lists = list(itertools.combinations_with_replacement(Sd, k))
  assert len(all_possible_snd_component_lists) == len(WC(k, factorial(d)))

  maximum = -1
  achieving_max = []        

  print(f"Iterating over all {len(all_possible_snd_component_lists)} sets to find max")
  for snd_component_list in loading_bar(all_possible_snd_component_lists):
    S = set(
      [(get_sparse_d_tuple(d), snd_comp) for snd_comp in snd_component_list]
    )

    val_of_sn = len(compute_Sn(S=S, n=n, group_op=op_wreath_prod))
    if val_of_sn > maximum: 
      maximum = val_of_sn 
      achieving_max = [snd_component_list]
    elif val_of_sn == maximum:
      achieving_max.append(snd_component_list)
    
    if noisy:
      print(f"{perms_to_str(snd_component_list)} gets value {format_large_num(val_of_sn)}")

  print(f"gamma({n=}, {k=}) = {format_large_num(maximum)}")
  print("Achieved with the following second components")
  for perm_list in achieving_max:
     print(perms_to_str(perm_list))
     

  return maximum

# return the set mentioned in the paper, with identity in the second component except for (12) and (123...d)
def get_paper_set(k, d):
    perm_id = tuple([i for i in range(0, d)])
    perm_12 = tuple([1, 0] + list(range(2, d))) 
    perm_cyc = tuple([i + 1 for i in range(0, d - 1)] + [0])
    
    return set(
       [(get_sparse_d_tuple(d), perm_id) for _ in range(k-2)] + [(get_sparse_d_tuple(d), perm_12), (get_sparse_d_tuple(d), perm_cyc)]
    )

def compute_Sn_paper_set(n, k, d):
    S = get_paper_set(k, d)
    s_n = compute_Sn(S, n, op_wreath_prod)
    return len(s_n)

# h_n(r) is the number of n \times n magic squares with sum r, as defined in the overleaf
# d = 3 throughout this file because we actually have a formula for h_3(r)
def h_3(r):
    return (r+1) * (r+2) * (r**2 + 3*r + 4) // 8

# f_d(k, n) = \sum_{x \in \WC(n, k)} \prod_{i=1}^{k} h_d(x_i) 
def f_3(k, n):
    WC_n_k = WC(n, k)
    tot = 0
    for x in WC_n_k:
        prod = 1 
        for i in range(0, k):
            prod *= h_3(x[i])
        tot += prod
    
    return tot

def lower_bound_3(k, n):
  tot = 0
  C_3 = 6
  for i in range(C_3, n):
     tot += (f_3(k-2, n-i) * (i - C_3 + 1))
  return tot

def see_how_tight_bounds_are(n, k):
  d = 3
  print(f"Upper bound: f_3(k, n) = {format_large_num(f_3(k, n))}")
  print(f"Lower bound = {format_large_num(lower_bound_3(k,n))}")
  print(f"S^n paper set: {format_large_num(compute_Sn_paper_set(n=n, k=k, d=d))}")

def see_how_tight_bounds_are_check_all_sets(n, k):
  d=3
  see_how_tight_bounds_are(n, k)
  compute_gamma_k_n(k = k, n = n, d = d)

# Since Z_d is a special case of S_d, we can test the S_d code by making sure a set 
# that only has "shift by one" permutations agrees with the formula we came up 
# with for Z_d
# using the given n, k, d, and t many components with 1 in the 2nd component
# with k - t elements having 0 in the 2nd component
def restrict_to_Z_d_works(n, k, d, p):
  def get_s_d_set():
    # shift_by_one, identity \in Z_d \subseteq S_d.
    shift_by_one = tuple(list(range(1, d)) + [0])
    identity = tuple(range(0, d))
    
    return set(
     [(get_sparse_d_tuple(d), shift_by_one) for _ in range(p)] + 
     [(get_sparse_d_tuple(d), identity) for _ in range(k-p)]  
    )
    
  val_from_s_d_sim = len(compute_Sn(S = get_s_d_set(), n = n, group_op = op_wreath_prod))
  val_from_z_d_formula = Z_d_product_Z_d.formula_when_all_zero_one_in_snd_comp(n=n, k=k, p=p, d=d)

  return val_from_s_d_sim == val_from_z_d_formula

def test_a_bunch_of_restrictions(n_upper, k_upper, d_upper):
  print("Testing S_d code by seeing that it lines up with Z_d")
  for n in incl_range(1, n_upper):
    for k in incl_range(1, k_upper):
      for d in incl_range(1, d_upper):
        for p in incl_range(0, k):
          if not restrict_to_Z_d_works(n=n, k=k, d=d, p=p):
            print("fail")

  print("success")

def main():
  #  see_how_tight_bounds_are(n=15, k=3)
  # see_how_tight_bounds_are_check_all_sets(n=15, k=3)
  # test_a_bunch_of_restrictions(6, 6, 3)
  pass

if __name__ == "__main__": 
    main()

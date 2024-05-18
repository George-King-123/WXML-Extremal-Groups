from shared_code import get_sparse_d_tuple, compute_Sn
from group_operations import op_wreath_prod
import Z_d_product_Z_d
from Z_d_wreath_product_S_d import compute_s_n_and_how_we_achieved_each_elt, make_tagged_set, generate_all_transpositions_in_s_d 
import D_infinity
from math import comb

def test_transpositions():
  def test_d_val(d):
    M = generate_all_transpositions_in_s_d(d)
    assert len(M) == comb(d, 2)

  for d in range(2, 20):
    test_d_val(d)

test_transpositions()

# Since Z_d is a special case of S_d, we can test the S_d code by making sure a set 
# that only has "shift by one" permutations agrees with the formula we came up 
# with for Z_d
# using the given n, k, d, and t many components with 1 in the 2nd component
# with k - t elements having 0 in the 2nd component
def restrict_to_Z_d(n, k, d, t):
  def get_s_d_set():
    # shift_by_one in Z_d
    shift_by_one = tuple([i for i in range(1, d)] + [0])
    identity = tuple([i for i in range(0, d)])
    
    S = set()
    
    # add something to the set k times 
    for i in range(0, k):
      if i < t:
        S.add((get_sparse_d_tuple(d), shift_by_one))
      else:
        S.add((get_sparse_d_tuple(d), identity))
    return S
  
  val_from_s_d_sim = len(compute_Sn(S = get_s_d_set(), n = n, group_op = op_wreath_prod))
  val_from_z_d_formula = Z_d_product_Z_d.formula_when_all_one(n = n, k =k, t=t, d=d)

  return val_from_s_d_sim == val_from_z_d_formula

# ex: test_a_bunch_of_restrictions(6, 6, 3, restrict_to_Z_d)
def test_a_bunch_of_restrictions(n_upper, k_upper, d_upper, restrict_fun):
  for n in range(1, n_upper + 1):
    for k in range(1, k_upper + 1):
      for d in range(3, d_upper + 1):
        for t in range(0, k + 1):
          if not restrict_fun(n=n, k = k, d= d, t = t):
            print("fail")

  print("success")

#test_a_bunch_of_restrictions(6, 6, 3, restrict_to_Z_d)

def test_compute_s_n_and_how_we_achieved_each_elt(n, k, d):
  tagged_set = make_tagged_set(k, d)
  elt_to_LL = compute_s_n_and_how_we_achieved_each_elt(tagged_set, n)
  total_paths = 0
  for elt in elt_to_LL.keys():
    total_paths += len(elt_to_LL[elt])
  return total_paths == k**n

#test_compute_s_n_and_how_we_achieved_each_elt_lots_of_vals(5, 5, 4)
def test_compute_s_n_and_how_we_achieved_each_elt_lots_of_vals(n_upper, k_upper, d_upper):
  for n in range (1, n_upper+1):
    for k in range (1, k_upper + 1):
      for d in range(1, d_upper + 1):
        if not test_compute_s_n_and_how_we_achieved_each_elt(n, k, d):
          print("failure")
          return
  print("success")


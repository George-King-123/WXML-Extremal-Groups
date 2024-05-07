from Z_d_product_Z_d import make_S, formula_when_all_one
from shared_code import compute_Sn
from group_operations import op_Z_d_Z_d

# test that the formula for just 0s and 1s in the second component
# works for all values of d,n,k up to the given upper bounds
# these have to be small because we explicitly compute |S^n| in order 
# to check the formula 
# ex: test_formula_all_one_many_vals(5, 5, 5)
def test_formula_all_one_many_vals(d_upper, n_upper, k_upper):
  def conj_works(n, k, t, d):
    S = make_S(d = d, k = k, distro_of_signs= [k - t, t] + ([0] * (d-2)))
    s_n = compute_Sn(S, n, op_Z_d_Z_d)
    return len(s_n) == formula_when_all_one(n, k, t, d)
  

  for d in range(2, d_upper + 1):
    for n in range (1, n_upper + 1):
      for k in range(1, k_upper + 1):
        for t in range(0, k+1):
          if not conj_works(n, k, t, d):
            print("we failed, t= {}".format(t))
            return
          
  print("success")

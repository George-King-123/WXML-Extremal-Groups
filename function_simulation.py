# mathematical idea (applies to this and function_simulation_fast):
# Assume S = {(x_1, e_1), ..., (x_k, e_k)} with x_1, ..., x_k a very 
# sparse set of integers. With this we can assume that there will no collisions,
# i.e. c_1x_1 + ... + c_k x_k = d_1 x_1 + ... + d_k x_k implies that c_i = d_i for all i
# thus we can represent these sums just as k tuples of the coefficients c_i
# so we can look at all the functions from [n] to [k] and compute these tuples. We then 
# put those sum tuples inside the larger tuple that contains the sign


# everything in this file is 1 indexed
# so we have to append a [0] to the start of some things, which we just ignore
# remember range(i, j) = [i, i+1, ..., j-1]

# define a function f: [n] -> [k] via a list of length n where each element is an integer 1-k
# we have an ignored [0] at the start, so f(i) is exactly f[i] 

# given an f: [n] -> [k], the values of n and k, and the number of negative signs in S, 
# compute the product ((c_1, ..., c_k), e)
def compute_prod(f, n, k, num_neg):
  # this is correct: e.g. if num_neg = k the second range is empty and adds nothing
  S = [0] + [-1 for i in range(0, num_neg)] + [1 for i in range (num_neg, k)]

  # uses the formula from the latex docuemnt
  # (it does take some thinking to line them up)
  components_of_tuple = [0] + [0] * k
  cur_sign = 1
  for i in range(1, n+1):
    components_of_tuple[f[i]] += cur_sign
    cur_sign *= S[f[i]]

  # the tuple won't have the extra [0] at the start
  k_tuple_sum = tuple(components_of_tuple[1:])
  return (k_tuple_sum, cur_sign)

# -this method is very simple: compute all possible functions from [n] to [k] and compute the
#   product that we get from that function, return the set of all possible products
# -runtime is at least O(k^n)
def compute_all_fcns_and_set(n, k, num_neg):
  all_tuples = set()
  def pick(lst):
    # our function is done
    if len(lst) == n+1:
      all_tuples.add(compute_prod(lst, n, k, num_neg))
      return 
    for i in range(1, k+1):
      pick(lst + [i])
  pick([0])
  return all_tuples

# computes |S^n| for |S| = k, choosing S
# to have all possible amonts of negative signs 
# between 0 and k
def run_sim(n, k):
  print("n = {} and k = {}".format(n, k))
  for negs in range(0, k+1):
    cardinality_s_n = len(compute_all_fcns_and_set(n, k, negs))
    print("|S^n| = {} when we have {} negative signs".format(cardinality_s_n, negs))
  

run_sim(7, 7)

  






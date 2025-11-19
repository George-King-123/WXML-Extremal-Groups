# run with python -m test.D_infinity_test_code

from D_infinity import compute_s_n_old_form
# import openpyxl

# verifies the formula for D_infinity given in the main file, checking
# all values of n, k up to the given upper bounds and all values of t 
# up to each k. 
# Since this actually runs the simulation, numbers need to be small e.g.
# run_check_on_formula(7, 7)
def run_check_on_formula(n_upper, k_upper, formula_fcn, noisy = False):
  def formula_works(n, k, t):
    formula = formula_fcn(n, k, t)
    empirical = compute_size_s_n_simulation(n, k, t)
    
    if noisy:
      if formula != empirical:
        print(f"formula = {formula:_}, empirical = {empirical:_}")
        print("formula doesn't work for n = {}, k = {}, t = {}".format(n, k, t))
        print()
      else: 
        print(f"formula works for n = {n}, k = {k}, t = {t}")
   
    return formula == empirical 

  works_on_all = True
  for n in range(1, n_upper + 1):
    for k in range(1, k_upper + 1):
      for t in range(0, k + 1):
        if not formula_works(n, k, t):
          works_on_all = False 
  
  if works_on_all:
    print("Formula worked on all values")
  else:
    print("formula failed, re-run with noisy if no other output")


# used to make the table I uploaded to the latex document that shows the 
# best number of signs for each n,k up to ~80. I banished it to the test code 
# file because no one should look at this code
def maximize_signs_excel_file(n_upper, k_upper, filename):
  wb = openpyxl.Workbook()
  ws = wb.active

  def write(r, c, data):
    cur_cell = ws.cell(row = r, column = c)
    cur_cell.value = data

  #ws.write(0, 1, "n")
  #ws.write(1, 0, "k")
  #ws.write(0, 0, "k \\ n")
  write(1, 1, "k \\ n")

  for i in range(2, n_upper + 1):
    write (1, i, str(i))


  for i in range(2, k_upper + 1):
    write(i, 1, str(i))

  def best_number_of_signs(n, k):
    biggest = -1
    achieved_at = 0
    tie_list = []
    for t in range(0, k + 1):
      val = compute_s_n_old_form(n, k, t)
      if val > biggest:
        biggest = val
        achieved_at = t
        tie_list = [str(achieved_at)]
      elif val == biggest:
        tie_list.append(str(t))
    write(k, n, ",".join(tie_list))
  

  for n in range(2, n_upper + 1):
    for k in range(2, k_upper + 1):
      best_number_of_signs(n = n, k = k)
    wb.save(filename)

def test_formula():
  run_check_on_formula(n_upper = 8, k_upper = 8, formula_fcn=compute_s_n_old_form, noisy=True)

def main():
  test_formula()

if __name__ == "__main__":
  main()

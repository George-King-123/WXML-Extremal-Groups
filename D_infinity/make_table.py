# run with python -m D_infinity.make_table

import D_infinity.compute_s_n as compute_s_n
from shared_code import incl_range
import pandas as pd 
import matplotlib.pyplot as plt



def find_maximizing_t_values(n, k):
    maximizers = []
    maximimum = -1 

    for t in incl_range(0, k): 
        cur_val = compute_s_n.formula(n=n, k=k, t=t)
        if cur_val > maximimum:
            maximizers = [t]
            maximimum = cur_val 
        elif cur_val == maximimum:
            maximizers.append(t) 
    
    return maximizers

def make_table(N_max, K_max) -> pd.DataFrame:
    table = pd.DataFrame(
        index=incl_range(2, N_max),
        columns=incl_range(2, K_max),
        dtype=object
    )

    # n = 1 or k = 1 is not interesting
    for n in incl_range(2, N_max):
        for k in incl_range(2, K_max):
            table.loc[k, n] = find_maximizing_t_values(n=n, k=k)

    return table 

# TODO: make more readable
def render_latex(table: pd.DataFrame): 
    table_str = table.map(lambda lst: ", ".join(map(str, lst)))

    # Set index name for LaTeX top-left cell
    table_str = table_str.reset_index().rename(columns={'index': r"$k \backslash n$"})

    num_cols = table_str.shape[1]
    col_format = "c|" + "c" * (num_cols - 1)
    latex_code = table_str.to_latex(index=False, escape=False, column_format=col_format)
    print(latex_code)

def main():
    render_latex(make_table(N_max=30, K_max=30))

if __name__ == "__main__":
    main()
# run with python -m D_infinity.make_table

import D_infinity.compute_s_n as compute_s_n
from shared_code import incl_range
import pandas as pd 

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
            table.loc[n, k] = find_maximizing_t_values(n=n, k=k)

    return table 

# TODO: make more readable
def render_latex(table: pd.DataFrame) -> str: 
    def get_column_format(): 
        CENTERED_COL = 'c'
        DIVIDER = '|'
        
        num_cols = table.shape[1]
        return CENTERED_COL + DIVIDER + (CENTERED_COL * (num_cols))

    def comma_separate_list(list_of_maximizers): 
        return ", ".join([str(t) for t in list_of_maximizers])
    
    def hacky_set_labels(df):
        return df.reset_index().rename(columns={'index': r"$n \backslash k$"})

    formatted_table = table.map(comma_separate_list)
    formatted_table = hacky_set_labels(formatted_table)

    latex_code = formatted_table.to_latex(index=False, escape=False, column_format=get_column_format())
    return latex_code

def main():
    print(render_latex(make_table(N_max=20, K_max=10)))

if __name__ == "__main__":
    main()
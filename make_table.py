# run with python -m make_table

import D_infinity.compute_s_n
import Z_d_product_Z_d
from shared_code import incl_range
from tqdm import tqdm
import pandas as pd 

def make_table(N_max, K_max, f_nk, n_min=2, k_min=2) -> pd.DataFrame:
    table = pd.DataFrame(
        index=incl_range(n_min, N_max),
        columns=incl_range(k_min, K_max),
        dtype=object
    )

    for n in incl_range(n_min, N_max):
        for k in tqdm(incl_range(k_min, K_max)):
            table.loc[n, k] = f_nk(n=n, k=k)

    return table 

# TODO: This function isn't pretty to look at, but does get the latex table looking right 
def render_latex(table: pd.DataFrame, highlight_diagonal=True) -> str: 
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

    if highlight_diagonal:
        for k_val in formatted_table:
            for n_val in formatted_table[k_val].index:
                if k_val == n_val:
                    formatted_table.loc[k_val, n_val] = "\\cellcolor{lightgray}" + str(formatted_table.loc[k_val, n_val])

    formatted_table = hacky_set_labels(formatted_table)

    latex_code = formatted_table.to_latex(index=False, escape=False, column_format=get_column_format())
    return latex_code

def d_infty():
    return render_latex(make_table(N_max=20, K_max=10, f_nk=D_infinity.compute_s_n.find_maximizing_p_values))

def s_d(d, n_min=2, k_min=2, N_max=20, K_max=10):
    def f(n, k):
        return Z_d_product_Z_d.find_maximizing_t_val(n=n, k=k, d=d)
    return render_latex(make_table(N_max=N_max, K_max=K_max, f_nk=f, n_min=n_min, k_min=k_min))

def main():
    print(d_infty())

if __name__ == "__main__":
    main()
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

def render_latex(table): 
    table_str = table.map(lambda lst: ", ".join(map(str, lst)))
    table_str.index.name = r"$k \backslash n$"
    print(table_str.to_latex(index=True, escape=False))

def render_matplot_lib(table): 
    # pivot from above

    fig, ax = plt.subplots(figsize=(max(6, 0.6*len(table.columns)), max(4, 0.5*len(table.index))))
    ax.axis('off')

    # cellText expects a 2D list of strings
    cell_text = table.astype(str).values.tolist()
    table = ax.table(
        cellText=cell_text,
        rowLabels=table.index.astype(str),
        colLabels=table.columns.astype(str),
        cellLoc='center',
        loc='center'
    )

    # table.auto_set_font_size(False)
    # table.set_fontsize(10)
    # table.auto_set_column_width(col=list(range(len(table.columns))))

    plt.tight_layout()
    plt.show()


def main():
    render_latex(make_table(N_max=10, K_max=10))
    # render_matplot_lib(make_table(N_max=10, K_max=10))

if __name__ == "__main__":
    main()
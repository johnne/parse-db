#!/usr/bin/env python

import matplotlib.pyplot as plt, pandas as pd
from matplotlib_venn import venn3
from argparse import ArgumentParser

def main():
    parser = ArgumentParser(description='''Create a Venn diagram plot of the Metacyc database for categories Biosynthesis, 
    Generation of Precursor Metabolites and Energy, and Degradation/Utilization/Assimilation''')
    
    parser.add_argument("-m", "--metacyctab", required=True,
            help="Tabular file of Metacyc database. See process_metacyc_db.py.")
    parser.add_argument("-o", "--outfile", default="metacyc_venn.png",
            help="Name of outfile plot. Defaults to metacyc_venn.png")

    args = parser.parse_args()

    metacyc_df = pd.read_csv(args.metacyctab, header=0, sep="\t", index_col=0)
    core = ["Generation of Precursor Metabolites and Energy"]
    syn = ["Biosynthesis"]
    deg = ["Degradation/Utilization/Assimilation"]
    ## Make Venn diagram
    c = set(metacyc_df.loc[(metacyc_df.Category1.isin(core)),"Pathway"])
    b = set(metacyc_df.loc[(metacyc_df.Category1.isin(syn)),"Pathway"])
    d = set(metacyc_df.loc[(metacyc_df.Category1.isin(deg)),"Pathway"])
    BC = b.intersection(c).difference(d)
    BD = b.intersection(d).difference(c)
    CD = c.intersection(d).difference(b)
    BCD = CD.intersection(b)
    B = b.difference(c.union(d))
    C = c.difference(d.union(b))
    D = d.difference(b.union(c))
    #make subsets as: (B,C,BC,D,BD,CD,BCD) 
    venn3(subsets=(len(B),len(C),len(BC),len(D),len(BD),len(CD),len(BCD)), set_labels=["Biosynthesis","Core","Degradation"])
    plt.title("#Pathway")
    plt.savefig(args.outfile, dpi=300, width=150, height=150, bbox_inches="tight")


if __name__ == '__main__':
    main()

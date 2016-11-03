#!/usr/bin/env python

import sys,pandas as pd, numpy as np
from argparse import ArgumentParser

def file_to_df(f):
    hin = open(f)
    r = []
    for line in hin:
        if line[0]=="#": continue
        line = line.rstrip()
        line = line.rsplit()
        items = [line[0],line[3],line[4],line[6],line[7],line[11],float(line[13]),int(line[17]),int(line[18])]
        r.append(items)
    hin.close()
    df = pd.DataFrame(r)
    df.columns = ["target","query_name","accession","seq_eval","seq_score","dom_eval","dom_score","from","to"]
    return df

def checkoverlap(r):
    alis = []
    store = []
    for i in range(0,len(r)): 
        f = int(r.iloc[i,6])
        t = int(r.iloc[i,7])
        this_ali = range(f,t+1)
        this_ali_len = len(this_ali)
        ## Check overlap
        if len(alis)==0: 
            store.append(i)
            alis.append(this_ali)
            continue
        overlapping = False
        for ali in alis:
            if len(ali)<this_ali_len: shortest = len(ali)
            else: shortest = this_ali_len
            o = set(this_ali).intersection(set(ali))
            if (float(len(o))/shortest)>0.1:
                overlapping = True
                break
                
        if not overlapping:
            alis.append(this_ali)
            store.append(i)
    return store

def parse_trusted(df,t):
    trusted = pd.read_csv(args.trusted, sep="\t", index_col=0, header=None, names=["Family","Sequence","Domain"], dtype={"Sequence":np.float64,"Domain":np.float64})
    ## Intersect with trusted
    ta = list(set(trusted.index).intersection(set(df.accession)))
    tn = list(set(trusted.index).intersection(set(df.query_name)))
    df1 = pd.merge(df,trusted.loc[ta],left_on="accession",right_index=True)
    df2 = pd.merge(df,trusted.loc[tn],left_on="query_name",right_index=True)

    dft = pd.concat([df1,df2])
    dft.index = dft.target
    dft.drop("target",axis=1, inplace=True)
    parsed = dft[dft.dom_score>=dft.Domain]
    return parsed

def main():
    parser = ArgumentParser(descrption='''Parses HMMER output (--tblout and --domtblout output) and reports 1 or 
    several non-overlapping hits matching thresholds per query''')

    parser.add_argument("-i", "--infile", required=True, type=str,
            help="HMMER results file")
    parser.add_argument("-e", "--evalue", type=float, default=1e-5,
            help="E-value to use as threshold. Defaults to 1e-5.")
    parser.add_argument("-t", "--trusted", type=str,
            help="Provide trusted score cutoffs for HMMs to parse with")
    
    df = file_to_df(f)
    
    if args.trusted: df = parse_trusted(df,args.trusted)
    else: df = df.loc[df.seq_eval<args.evalue]

    ## Count hits
    c = pd.DataFrame(df.groupby(level=0).count().ix[:,0])
    
    ## Get queries with single hits
    single = list(c[c.query_name==1].index)

    ## Get queries with multiple hits
    multi = list(set(c[c.query_name>1].index))

    d = df.loc[single]
    for gene in multi:
        r = df.loc[gene]
        rs = r.sort_values("dom_score",ascending=False)
        store = checkoverlap(rs)
        d = pd.concat([d,rs.iloc[store]])
    
    ## Join multiple non-overlapping hits for each query
    df_join = d.ix[:,0:2].groupby(level=0).agg(lambda x: '|'.join(set(x)))
    
    df_join.to_csv(sys.stdout, sep="\t")
    
if __name__ == '__main__':
    main()

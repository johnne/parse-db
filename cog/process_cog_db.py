#!/usr/bin/env python

import sys, pandas as pd, logging, urllib

def make_coginfo(fun="ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/fun2003-2014.tab",names="ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cognames2003-2014.tab"): 
    try: 
        logging.info("Reading info from "+fun)    
        cogfun = pd.read_csv(fun, sep="\t",header=0,index_col=0,names=["Code","Name"])
        logging.info("Reading info from "+names)
        cognames = pd.read_csv(names, sep="\t", header=0,index_col=0, usecols=[0,1], names=["COG","Code"])
    except urllib.error.URLError: sys.exit("Connection failed")

    coginfo = pd.DataFrame(columns=["Code","Name"])
    for cog in cognames.index:
        codes = list(cognames.loc[cog,"Code"])
        funcs = list(cogfun.loc[codes,"Name"])
        tmp = pd.DataFrame(columns=coginfo.columns,index=[cog]*len(codes),data={"Code":codes,"Name":funcs})
        coginfo = pd.concat([coginfo,tmp])
    return coginfo

def main():
    logging.basicConfig(format='%(levelname)s:%(message)s',level=logging.INFO)
    coginfo = make_coginfo()
    coginfo.to_csv(sys.stdout, sep="\t")

if __name__ == '__main__':
    main()

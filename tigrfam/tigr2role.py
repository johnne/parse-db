#!/usr/bin/env python

from argparse import ArgumentParser
import pandas as pd, sys, logging

def main():
    parser = ArgumentParser()
    parser.add_argument("-l", "--link", type=str,required=True,
            help="TIGRFAM role link table")
    parser.add_argument("-n", "--names", type=str, required=True,
            help="TIGRFAM role name table")
    args = parser.parse_args()

    rl = pd.read_csv(args.link, sep="\t", header=None, names=['tigrfam','role_id'])
    rn = pd.read_csv(args.names, sep="\t", header=None, usecols=[1,2,3], names=['role_id','role_type','role_name'])
    tr = pd.DataFrame(columns=['tigrfam','mainrole','sub1role'])
    
    for i,fam in enumerate(list(set(rl.tigrfam))):
        role_id = list(rl.loc[rl.tigrfam==fam,'role_id'])[0]
        try: mainrole = list(rn.loc[(rn.role_id==role_id) & (rn.role_type=='mainrole:'),'role_name'])[0]
        except IndexError: mainrole="Unknown"
        try: sub1role = list(rn.loc[(rn.role_id==role_id) & (rn.role_type=='sub1role:'),'role_name'])[0]
        except IndexError: sub1role="Unknown"
        tmp = pd.DataFrame(columns=['tigrfam','mainrole','sub1role'],index=[i],data={'tigrfam':fam,'mainrole':mainrole,'sub1role':sub1role})
        tr = pd.concat([tr,tmp])
    tr.index = tr.tigrfam
    tr.drop('tigrfam',axis=1,inplace=True)
    tr.to_csv(sys.stdout,sep="\t")
    

if __name__ == '__main__':
    main()

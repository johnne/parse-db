#!/usr/bin/env python

import pandas as pd, sys, logging
from argparse import ArgumentParser
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

def add_term_data(terms,term_ids,term_name,parents):
    parents = list(set(parents))
    index = term_ids*len(parents)
    tmp = pd.DataFrame(columns=["name","parent"],index=index, data={"name": [term_name]*len(index),
        "parent":parents*len(term_ids)})
    terms = pd.concat([terms,tmp])
    return terms

def parse_terms(infile):
    lines = i = 0
    terms = pd.DataFrame(columns=["name","parent"])
    term_ids = []
    term_name = ""
    parents = []
    top_parents = []
    with open(infile, 'r', encoding="utf-8", errors="ignore") as f:
        for line in f:
            lines+=1
            i+=1
            if i>=10000:
                logging.info("Processed "+str(lines)+" lines")
                i=0
            if line[:6]=="[Term]":
                if term_ids!=[]:
                    if parents ==[]: parents = ['root']
                    terms = add_term_data(terms,term_ids,term_name,parents)
                term_ids = []
                name = "" 
                parents = []
            elif line[:4]=="id: " or line[:8]=="alt_id: ": term_ids.append(line.rstrip().split(" ")[-1])
            elif line[:6]=="name: ": term_name = line.rstrip()[6:]
            elif line[:6]=="is_a: ": parents.append(line.rstrip().split(" ")[1])
    logging.info("Processed "+str(lines)+" lines")
    return terms

def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", type=str, required=True,
            help="Gene ontology basic file ('go-basic.obo')")
    args = parser.parse_args()

    logging.info("Reading go file")
    terms = parse_terms(args.infile)
    terms.to_csv(sys.stdout,sep="\t")

if __name__ == "__main__":
    main()

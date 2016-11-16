#!/usr/bin/env python

import pandas as pd, sys, logging
from argparse import ArgumentParser
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

def find_all_parents(graph, start, end, path=[]):
        path = path + [start]
        if start == end:
            return [path]
        if not start in graph.keys():
            return []
        paths = []
        for node in graph[start]:
            if node not in path:
                newpaths = find_all_parents(graph, node, end, path)
                for newpath in newpaths:
                    paths.append(newpath)
        return paths

def make_graph(df):
    logging.info("Creating graph")    
    grouped = df.loc[:,["term","parent"]].groupby('term')
    dfg = grouped.aggregate(lambda x: list(x))
    graph = dfg.loc[:,"parent"].to_dict()
    return graph

def get_name(term,df):
    try: return list(set(df.loc[df.term==term,"name"]))[0]
    except IndexError: return term

def write_hierarchy(df,graph,goids,root):
    r = df.groupby("term").first().loc[:,"name"]
    logging.info("Creating hierarchy")
    sys.stdout.write("term\tname\tparent_names\tparent_ids\n")
    if goids: terms = goids
    else: terms = list(graph.keys())
    total = str(len(terms))
    j = 0
    processed = []
    for i,term in enumerate(terms,start=1):
        j+=1
        if j==500:
            logging.info("Processed "+str(i)+"/"+total)
            j = 0
        name = get_name(term,df)
        paths = find_all_parents(graph,term,root)
        for path in paths:
            path = path[1:-1]
            path.reverse()
            names = "|".join(list(r.loc[path].values))
            sys.stdout.write(term+"\t"+name+"\t"+names+"\t"+"|".join(path)+"\n")

def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", required=True,
            help="Infile table version of Gene Ontology database with columns 'term_id,name,parent'")
    parser.add_argument("-g", "--goids", nargs="*",
            help="List of Gene ontology ids to create hierarchy table for")
    parser.add_argument("-r", "--root", default='root', type=str,
            help="End point of hierarchy (GO term id). Defaults to 'root' which generates the full hierarchy")

    args = parser.parse_args()

    df = pd.read_csv(args.infile, header=0, sep="\t",names=["term","name","parent"])
    logging.info("Read table with "+str(len(df))+" rows")
    if args.goids: logging.info("Read "+str(len(args.goids))+" GO term ids")
    graph = make_graph(df)
    
    write_hierarchy(df,graph,args.goids,args.root)

if __name__ == '__main__':
    main()

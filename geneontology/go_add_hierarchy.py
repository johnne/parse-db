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
    graph = {}
    for i in df.index:
        term = df.loc[i,"term"]
        parent = df.loc[i,"parent"]
        try: graph[term].append(parent)
        except KeyError: graph[term] = [parent]
    return graph

def get_name(term,df):
    try: return list(set(df.loc[df.term==term,"name"]))[0]
    except IndexError: return term

def write_hierarchy(df,graph):
    sys.stdout.write("term"+"\t"+"name"+"\t"+"parent_names"+"\t"+"parent_ids\n")
    for term in graph:
        name = get_name(term,df)
        paths = find_all_parents(graph,term,'root')
        for path in paths:
            path = path[1:-1]
            path.reverse()
            names = []
            for item in path: names.append(get_name(item,df))
            sys.stdout.write(term+"\t"+name+"\t"+"|".join(names)+"\t"+"|".join(path)+"\n")

def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", required=True,
            help="Infile table version of Gene Ontology database with columns 'term_id,name,parent'")
    parser.add_argument("-g", "--goids", nargs="*",
            help="List of Gene ontology ids to create hierarchy table for")

    args = parser.parse_args()

    df = pd.read_csv(args.infile, header=0, sep="\t",names=["term","name","parent"])
    logging.info("Read table with "+str(len(df))+" rows")
    if args.goids: 
        df = df.loc[(df.term.isin(args.goids)) | (df.parent.isin(args.goids))]
        logging.info("Trimmed table to "+str(len(df))+" rows")

    graph = make_graph(df)
    
    write_hierarchy(df,graph)

if __name__ == '__main__':
    main()

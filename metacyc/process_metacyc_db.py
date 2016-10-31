#!/usr/bin/env python

import pandas as pd
from argparse import ArgumentParser

def ec_parse(line):
    line = line.rstrip()
    ec = line.split("-")[-1]
    ec+=".-"*(3-ec.count("."))
    return ec

def add_ec_data(ec2path,ecs,pwys):
    ecs = list(set(ecs))
    pwys = list(set(pwys))
    df = pd.DataFrame(columns=["Pathway"],index=ecs*len(pwys), data={"Pathway":pwys*len(ecs)})
    return pd.concat([ec2path,df])

def parse_reactions(reactions):
    ec2path = pd.DataFrame(columns=["Pathway"])
    ecs = []
    pwys = []
    with open(reactions,'r', encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line[:2]=="//": 
                if len(pwys)>0: ec2path = add_ec_data(ec2path,ecs,pwys)
                ecs = []
                pwys = []
            elif line[:12]=="EC-NUMBER - ": ecs.append(ec_parse(line))
            elif line[:13]=="IN-PATHWAY - ": pwys.append(line.rstrip().split(" ")[-1])
    return ec2path

def add_cpd_data(cpds,cpd,name,ty):
    tmp = pd.DataFrame(columns=["Name","Type"],index=[cpd],data={"Name":name,"Type":ty})
    cpds = pd.concat([cpds,tmp])
    return cpds

def parse_compounds(compounds):
    cpds = pd.DataFrame(columns=["Name","Type"])
    cpd = name = ty = ""
    with open(compounds, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line[:2]=="//": 
                cpds = add_cpd_data(cpds,cpd,name,ty)
                cpd = name = ty = ""
            elif line[:12]=="UNIQUE-ID - ": cpd = line.rstrip()[12:]
            elif line[:8]=="TYPES - ": ty = line.rstrip()[8:]
            elif line[:14]=="COMMON-NAME - ": name = line.rstrip()[14:].replace("<i>","").replace("</i>","")
            
    return cpds

def add_pwy_data(pwys,pathway,name,parents):
    parents = list(set(parents))
    tmp = pd.DataFrame(columns=["Name","Parent"],index=[pathway]*len(parents), data={"Name": [name]*len(parents),
                                                                                     "Parent":parents})
    pwys = pd.concat([pwys,tmp])
    return pwys

def add_links(links,pathway,line,p):
    items = line.rstrip()[17:].replace("(","").replace(")","").split(" ")
    cpd = items[0]
    try: links[pathway]
    except KeyError: links[pathway] = {}
        
    for pwy_link in items[1:]:
        if not pwy_link in p: continue
        
        try: links[pwy_link]
        except KeyError: links[pwy_link] = {}
            
        try: 
            links[pathway][pwy_link].append(cpd)
            links[pathway][pwy_link] = list(set(links[pathway][pwy_link]))
        except KeyError: links[pathway][pwy_link] = [cpd]
            
        try: 
            links[pwy_link][pathway].append(cpd)
            links[pwy_link][pathway] = list(set(links[pwy_link][pathway]))
        except KeyError: links[pwy_link][pathway] = [cpd]
    return links        

def parse_pathways(pathways,p):
    pwys = pd.DataFrame(columns=["Parent"])
    links = {}
    pathway = name = ""
    parents = []
    with open(pathways, 'r', encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line[:2]=="//": 
                pwys = add_pwy_data(pwys,pathway,name,parents)
                pathway = name = "" 
                parents = []
            elif line[:12]=="UNIQUE-ID - ": pathway = line.rstrip()[12:]
            elif line[:8]=="TYPES - ": parents.append(line.rstrip()[8:])
            elif line[:14]=="COMMON-NAME - ": name = line.rstrip()[14:].replace("<i>","").replace("</i>","")
            elif line[:17]=="PATHWAY-LINKS - (": links = add_links(links,pathway,line,p)
    return (pwys,links)

def parse_links(links,cpds):
    link_count = pd.DataFrame(columns=["Count"])
    graph_df = pd.DataFrame(columns=["Node1","Node2","CPD"])
    for key1 in links.keys():
        nodes = []
        edges = []
        for key2 in links[key1]:
            nodes.append(key2)
            cpd_names = []
            for c in links[key1][key2]:
                try: cpd_names.append(cpds.loc[c,"Name"])
                except KeyError: cpd_names.append(c)
            edges.append("|".join(cpd_names))
        tmp = pd.DataFrame(columns=["Node1","Node2","CPD"], data={"Node1": [key1]*len(nodes), "Node2": nodes, "CPD": edges})
        graph_df = pd.concat([graph_df,tmp])

    graph_df.index = list(range(0,len(graph_df)))
    return graph_df

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

def create_parent_hier(cl_graph, cl_names, pwy_parents, top):
    all_parents = pd.DataFrame(columns=["Class","Hierarchy"])
    max_levels = 0
    for p in pwy_parents:
        parent_lists = find_all_parents(cl_graph,p,top)
        if parent_lists == []: continue
        for parents in parent_lists:
            parents.reverse()
            for parent in parents: 
                parents[parents.index(parent)] = cl_names[parent]
                if len(parents)>max_levels: max_levels = len(parents)
                
            tmp = pd.DataFrame(index=[0],columns=["Class","Hierarchy"],data={"Class": p, "Hierarchy": "|".join(parents).lstrip("|")})
            all_parents = pd.concat([all_parents,tmp])
    all_parents.index = list(range(0,len(all_parents))) 
    
    cls = pd.DataFrame(columns=["Category"+str(x) for x in list(range(1,max_levels+1))])
    for i in all_parents.index:
        c = all_parents.loc[i,"Class"]
        h = all_parents.loc[i,"Hierarchy"]
        items = h.split("|")
        items+=[items[-1]]*(max_levels-len(items))
        data = {}
        for i,item in enumerate(items,start=1):data["Category"+str(i)] = item
        tmp = pd.DataFrame(columns=cls.columns,index=[c],data=data)    
        cls = pd.concat([cls,tmp])
    return cls
    
def parse_classes(classes, pwy_parents):
    cl_names = {}
    cl_graph = {} 
    cl = name = ""
    parents = []
    with open(classes, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line[:2]=="//": 
                cl_names[cl] = name
                cl_graph[cl] = parents
                cl = name = ""
                parents = []             
            elif line[:12]=="UNIQUE-ID - ": cl = line.rstrip()[12:]
            elif line[:8]=="TYPES - ": parents.append(line.rstrip()[8:])
            elif line[:14]=="COMMON-NAME - ": name = line.rstrip()[14:].replace("<i>","").replace("</i>","")
    top = "Pathways"
    cls = create_parent_hier(cl_graph, cl_names, pwy_parents, top)
    return cls

def main():
    parser = ArgumentParser()
    parser.add_argument("-d", "--datadir", type=str, default="./", 
            help="Metacyc data dir, containing classes.dat, reactions.dat, pathways.dat and compounds.dat. Defaults to current directory.")

    args = parser.parse_args()

    reactions = args.datadir+"/reactions.dat"
    classes = args.datadir+"/classes.dat"
    pathways = args.datadir+"/pathways.dat"
    compounds = args.datadir+"/compounds.dat"

    ec2path = parse_reactions(reactions)

    cpds = parse_compounds(compounds)

    p = list(set(ec2path.Pathway))
    (pwys,links) = parse_pathways(pathways,p)

    graph_df = parse_links(links, cpds)

    pwy_parents = list(set(pwys.Parent))
    cls = parse_classes(classes,pwy_parents)

    df = pd.merge(pwys,cls,left_on="Parent",right_index=True)
    df = pd.merge(ec2path,df,left_on="Pathway",right_index=True)
    cols = [0,1]+list(range(3,len(df.columns)))
    df.ix[:,cols]

    df.to_csv("metacyc.ec2pathcats.tab", sep="\t")
    graph_df.to_csv("metacyc.pathway_graph.tab", sep="\t")

if __name__ == "__main__":
    main()

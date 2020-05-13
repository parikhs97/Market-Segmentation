from igraph import Graph
import pandas as pd 
import sys
from scipy import spatial
import collections

#creating graphs 
def createGraph():
    attributes = pd.read_csv('data/fb_caltech_small_attrlist.csv').as_matrix()
    graph = Graph.Read_Edgelist('data/fb_caltech_small_edgelist.txt', directed=False)
    return attributes,graph

#list of members of the grpah 
def members(graph):
    member = []
    for g in range(graph.vcount()):
        member.append(g)
    return member
    
#phase 1 
def phase1(graph,members,attributes,alpha):
    l = len(members)
    for i in range(l):
        maxPositiveGain = 0.0
        maxMember = -1
        for j in range(l):
            if members[i] == members[j]:
                continue
            
            #cosine similarity
            indices = []
            for x, y in enumerate(members):
                if y== members[j]:
                    indices.append(x)
            sim = 0 
            for ele in indices:
                sim = sim + spatial.distance.cosine(attributes[i],attributes[ele])
            sim = sim / len(indices)
            
            temp = members.copy()
            Q_NewMan_i = graph.modularity(temp)
            temp[i] = members[j]
            Q_NewMan_j = graph.modularity(temp)
            Delta_Q_NewMan = Q_NewMan_j - Q_NewMan_i
            
            Delta_Q_NewMan_Attr = sim 
            Delta_Q1 = (alpha * Delta_Q_NewMan)
            Delta_Q2 = (float(1 - alpha) * Delta_Q_NewMan_Attr)
            Delta_Q = ( Delta_Q1 + Delta_Q2 )
            
            if Delta_Q > maxPositiveGain:
                    maxPositiveGain = Delta_Q
                    maxMember = j
            
        if maxPositiveGain > 0 :
            members[i] = members[maxMember]

    return members

#phase 2
def phase2(graph,members,attributes,alpha):
    graph.contract_vertices(members,combine_attrs="mean")
    graph.simplify(combine_edges=sum,multiple = True, loops = True)
    phase2Members = phase1(graph,members,attributes,alpha)
    return phase2Members

#creating clusters(dictionary)
def createCluster(phase2Members):
    clusters = collections.defaultdict(list)
    for i in range(len(phase2Members)):
        clusters[phase2Members[i]].append(i)
    return clusters
    
#result into test file 
def outputFile(clusters,alpha):
    if alpha == 0.5: 
        fileName = "communities_"+str(5)+".txt"
    else:
        fileName = "communities_"+str(alpha)+".txt"
    file = open(fileName,"w+")
    for c in clusters.values():
        for n in range(len(c)):
            if n != len(c)-1:
                file.write(str(c[n]) + ",")
            else:
                file.write(str(c[n]))    
        file.write("\n")
    file.close()
    
 

alpha = float(sys.argv[1])     
attributes,graph = createGraph()
member= members(graph)
phase1Members = phase1(graph,member,attributes,alpha)
#print(len(set(phase1Members)))
phase2Members = phase2(graph,member,attributes,alpha)
#print(len(set(phase2Members)))
clusters = createCluster(phase2Members)
#print(clusters)
outputFile(clusters,alpha)
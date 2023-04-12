import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import sys
import numpy as np
import networkx as nx
from tkinter import *
root=Tk()
INT_MAX = sys.maxsize
######################################################## Filing ################################################
def returngraphtooriginal(edges,G):
    Temp = list(G.edges.keys())
    G.remove_edges_from(Temp)
    G.add_weighted_edges_from(edges)



global FILE_NAME,dif_x,dif_src,dif_directed,dif_undirected
FILE_NAME=""
dif_x=0
dif_src=0
dif_directed=[[]]
dif_undirected=[[]]

class FileExtract:

    def getgraphing(self,fname,G):

        global FILE_NAME,dif_x,dif_directed,dif_undirected,dif_src
        if(True):
            f = open(fname, 'r')
            for _ in range(2):
                next(f)
            x = -1000
            i = 0
            directed = [[]]
            undirected = [[]]
            for line in f.readlines():
                if (int(x) == -1000):
                    x = int(line)
                    # print(x)
                    directed = np.zeros((x, x))
                    undirected = np.zeros((x, x))
                    indices_zero = directed == 0
                    indices_zero = undirected == 0
                    undirected[indices_zero] = INT_MAX
                    directed[indices_zero] = INT_MAX
                    # print(k)


                elif (int(i) < int(x)):
                    line = line.replace('\n', '')
                    if line != '':
                        fields = line.split('\t')
                        # print(i)
                        # print(fields)
                        G.add_node(int(fields[0]), pos=(float(fields[1]), float(fields[2])))
                        #positions.append([float(fields[i]),float(fields[2])])
                        i += 1


                elif (int(i) >= int(x)):
                    line = line.replace('\n', '')
                    if line != '':
                        fields1 = line.split('\t')

                        if len(fields1) != 1:
                            fields1.pop()
                            # print(fields1)
                            z = int(fields1[0])

                            j = 1
                            while (j < len(fields1)):

                                if (j != 0):
                                    # print(z, int(float(fields1[j])), int(float(fields1[j + 2])))
                                    if z != int(fields1[j]):
                                        directed[z][int(fields1[j])] = float(fields1[j + 2]) / 10 ** 7
                                        G.add_edge(z, int(fields1[j]))
                                        if directed[int(fields1[j])][z] == 0:
                                            directed[int(fields1[j])][z] = float(fields1[j + 2]) / 10 ** 7

                                    if undirected[z][int(fields1[j])] > float(fields1[j + 2]) / 10 ** 7:
                                        undirected[z][int(fields1[j])] = float(fields1[j + 2]) / 10 ** 7
                                        undirected[int(fields1[j])][z] = float(fields1[j + 2]) / 10 ** 7

                                    j += 4
                        else:

                            b = int(fields1[0])
            f.close()

            #pos = nx.get_node_attributes(G, 'pos')

            #nx.draw(G, pos, with_labels=1, arrows=False, node_color='green')
            FILE_NAME=fname
            dif_x = x
            dif_directed = directed
            dif_undirected = undirected
            dif_src = b
            return [x, b, directed, undirected]
        else:

            return [dif_x, dif_src, dif_directed, dif_undirected]

###################################################################################################################
################################## Dijkstra #######################################################################
class Graph():
 
    def __init__(self, vertices):
        self.V = vertices
        self.graph = [[0 for column in range(vertices)]
                    for row in range(vertices)]
 
    def printSolution(self, dist, source):
        text1.insert(INSERT,("Vertex \tDistance from Source: %d\n" %source))
        #print ("Vertex \tDistance from Source")
        for node in range(self.V):
            if(dist[node]==INT_MAX):
                text1.insert(INSERT,("%d \t Infinity\n" %(node)))
                #print (node, "\t", "Infinity")
            else:
                text1.insert(INSERT,("%d \t %d\n" %(node, dist[node])))
                #print (node, "\t", dist[node])
 
    def minDistance(self, dist, sptSet):
 
        min = INT_MAX
        min_index=0
        for u in range(self.V):
            if dist[u] < min and sptSet[u] == False:
                min = dist[u]
                min_index = u
 
        return min_index
 
    def dijkstra(self, src):
 
        dist = [INT_MAX] * self.V
        dist[src] = 0
        sptSet = [False] * self.V
 
        for cout in range(self.V):
 
            x = self.minDistance(dist, sptSet)
 
            sptSet[x] = True
 
            for y in range(self.V):
                if self.graph[x][y] > 0 and sptSet[y] == False and \
                dist[y] > dist[x] + self.graph[x][y]:
                        dist[y] = dist[x] + self.graph[x][y]
 
        self.printSolution(dist,src)
#################################################################################################################
################################floyd warshal################################################
INF = INT_MAX
def floyd(G,nV):
    dist = G
    for r in range(nV):
        for p in range(nV):
            for q in range(nV):
                dist[p][q] = min(dist[p][q], dist[p][r] + dist[r][q])
    sol(dist,nV)

def sol(dist,nV):
    for p in range(nV):
        text1.insert(INSERT,("Source: %d\n" %p))
        #print("Source: ",end=" ")
        #print(p)
        for q in range(nV):
            if(dist[p][q] == INF):
                text1.insert(INSERT,("INF "))
                #print("INF", end=" ")
            else:
                text1.insert(INSERT,("%d " %dist[p][q]))
                #print(dist[p][q], end="  ")
        text1.insert(INSERT,"\n")
        print("\n")
###############################################################################################################
########################################Prims##################################################################
def isValidEdge(u, v, inMST):
    if u == v:
        return False
    if inMST[u] == False and inMST[v] == False:
        return False
    elif inMST[u] == True and inMST[v] == True:
        return False
    return True
 
def primMST(cost,source,V,G):
    inMST = [False] * V
    graphedges=[]
    inMST[0] = True
 
    edge_count = 0
    mincost = 0
    while edge_count < V - 1:
 
        minn = INT_MAX
        a = -1
        b = -1
        for i in range(V):
            for j in range(V):
                if cost[i][j] < minn:
                    if isValidEdge(i, j, inMST):
                        minn = cost[i][j]
                        a = i
                        b = j
 
        if a != -1 and b != -1:
            text1.insert(INSERT,("Edge %d: (%d, %d) cost: %d\n" %(edge_count,a,b,minn)))
            #print("Edge %d: (%d, %d) cost: %d" %
            #     (edge_count, a, b, minn))
            graphedges.append([a,b,minn])
            edge_count += 1
            mincost += minn
            inMST[b] = inMST[a] = True
    text1.insert(INSERT,("Minimum cost = %d\n" %mincost))
    returngraphtooriginal(graphedges,G)
    plt.title('Prims MST')
    #print("Minimum cost = %d" % mincost)
###############################################################################################################
################################################ Kruskal's Algo ###############################################
def find(i,parent):
    while parent[i] != i:
        i = parent[i]
    return i
 
def union(i, j, parent):
    a = find(i,parent)
    b = find(j,parent)
    parent[a] = b
 
def kruskalMST(cost,V,G):
    graphedges=[]
    parent = [i for i in range(V)]
    mincost = 0 # Cost of min MST
 
    for i in range(V):
        parent[i] = i
 
    edge_count = 0
    while edge_count < V - 1:
        min = INF
        a = -1
        b = -1
        for i in range(V):
            for j in range(V):
                if find(i,parent) != find(j,parent) and cost[i][j] < min:
                    min = cost[i][j]
                    a = i
                    b = j
        union(a, b, parent)
        text1.insert(INSERT,('Edge {}:({}, {}) cost:{}\n'.format(edge_count, a, b, min)))
        #print('Edge {}:({}, {}) cost:{}'.format(edge_count, a, b, min))
        graphedges.append([a,b,min])
        edge_count += 1
        mincost += min
    text1.insert(INSERT,("Minimum cost= {}\n".format(mincost)))
    returngraphtooriginal(graphedges,G)
    plt.title('Kruskal MST')
    #print("Minimum cost= {}".format(mincost))


def getglobaldata(filename,G):
    F1 = FileExtract()
    GlobalData = F1.getgraphing(filename,G)
    return GlobalData

###################################################################################################################
########################################## BellmanFord ############################################################
def BellmanFord(graph, V, E, src):
 
    dis = [sys.maxsize] * V
 
    dis[src] = 0
 
    for i in range(V - 1):
        for j in range(E):
            if dis[graph[j][0]] + \
                   graph[j][2] < dis[graph[j][1]]:
                dis[graph[j][1]] = dis[graph[j][0]] + \
                                       graph[j][2]
 
    for i in range(E):
        x = graph[i][0]
        y = graph[i][1]
        weight = graph[i][2]
        if dis[x] != sys.maxsize and dis[x] + \
                        weight < dis[y]:
            text1.insert(INSERT,("Graph contains negetive weight cycle\n"))
            #print("Graph contains negative weight cycle") 
    text1.insert(INSERT,("Vertex Distance from Source\n"))
    #print("Vertex Distance from Source")
    for i in range(V):
        if(dis[i]==sys.maxsize):
            text1.insert(INSERT,("%d\t\tInfinity\n" % (i)))
        else:
            text1.insert(INSERT,("%d\t\t%d\n" % (i, dis[i])))
        #print("%d\t\t%d" % (i, dis[i]))
def edges(G,V):
    edges=[]
    count=0
    for i in range(V):
        for j in range(V):
            if(G[i][j]!=INT_MAX):
                list=[i,j,G[i][j]]
                edges.append(list)
                count=count+1
    edgesandcount=[edges,count]
    return edgesandcount

#########################################################################################################################
################################################## BoruvkaMST ###########################################################
def findboruvka(parentboruvka, i):
    if parentboruvka[i] == i:
        return i
    return findboruvka(parentboruvka, parentboruvka[i])

def unionboruvka(parentboruvka, rank, x, y):
        xroot = findboruvka(parentboruvka, x)
        yroot = findboruvka(parentboruvka, y)

        if rank[xroot] < rank[yroot]:
            parentboruvka[xroot] = yroot
        elif rank[xroot] > rank[yroot]:
            parentboruvka[yroot] = xroot
        else :
            parentboruvka[yroot] = xroot
            rank[xroot] += 1


def boruvkaMST(graph,V,G):
        graphedges=[]
        edgecount=0
        parentboruvka = []; rank = []; 
        cheapest =[]
        numTrees = V
        MSTweight = 0
        for node in range(V):
            parentboruvka.append(node)
            rank.append(0)
            cheapest =[-1] * V
        while numTrees > 1:
            for i in range(len(graph)):
                u,v,w =  graph[i]
                set1 = findboruvka(parentboruvka, u)
                set2 = findboruvka(parentboruvka ,v)
                if set1 != set2:     
                    if cheapest[set1] == -1 or cheapest[set1][2] > w :
                        cheapest[set1] = [u,v,w] 
                    if cheapest[set2] == -1 or cheapest[set2][2] > w :
                        cheapest[set2] = [u,v,w]  
            for node in range(V):
                if cheapest[node] != -1:
                    u,v,w = cheapest[node]
                    set1 = findboruvka(parentboruvka, u)
                    set2 = findboruvka(parentboruvka ,v)
                    if set1 != set2 :
                        MSTweight += w
                        unionboruvka(parentboruvka, rank, set1, set2)
                        text1.insert(INSERT,("Edge %d: (%d, %d) cost: %d\n" %(edgecount,u,v,w)))
                        edgecount=edgecount+1;
                        graphedges.append([u,v,w])
                        #print ("Edge %d-%d with weight %d included in MST" % (u,v,w))
                        numTrees = numTrees - 1
            cheapest =[-1] * V            
        text1.insert(INSERT,("Minimum cost = %d\n" %MSTweight))
        returngraphtooriginal(graphedges,G)
        plt.title('Boruvka MST')
        #print ("Weight of MST is %d" % MSTweight)
###################################################################################################################


def setGraph():
    text1.delete(1.0,END)
    #print(preVInput)
    #if(preVInput != int(clicked1.get())):
    G=nx.Graph()
    filename="input"+clicked1.get()+".txt"
    list=getglobaldata(filename,G)
    undirected=list[2]
    plt.title('Input Graph')
    if (clicked2.get()=="Dijkstra"):
        #print("Inserted")
        #list=getglobaldata(filename)
        graph1=Graph(list[0])
        graph1.graph=list[2]
        graph1.dijkstra(list[1])
    elif (clicked2.get()=="FloydWarshal"):
        #list=getglobaldata(filename)
        graph2=Graph(list[0])
        graph2.graph=list[2]
        indices_zero=graph2.graph == INT_MAX
        graph2.graph[indices_zero]=INF
        floyd(graph2.graph,list[0])
    elif (clicked2.get()=="Prims"):
        #list=getglobaldata(filename,G)
        graph3=Graph(list[0])
        graph3.graph=list[2]
        primMST(graph3.graph,list[1],list[0],G)
    elif (clicked2.get()=="Kruskal"):
        #list=getglobaldata(filename)
        graph4=Graph(list[0])
        graph4.graph=list[2]
        kruskalMST(graph4.graph,list[0],G)
    elif (clicked2.get()=="BellmanFord"):
        #list=getglobaldata(filename)
        graph5=Graph(list[0])
        graph5.graph=list[2]
        edgescount=edges(graph5.graph,list[0])
        BellmanFord(edgescount[0],list[0],edgescount[1],list[1])
    elif( clicked2.get()=="Boruvka"):
        #list=getglobaldata(filename)
        graph6=Graph(list[0])
        graph6.graph=list[2]
        edgescount=edges(graph6.graph,list[0])
        boruvkaMST(edgescount[0],list[0],G)
    elif (clicked2.get()=="Clustering_Coefficient"):
        text1.insert(INSERT,("Clustering Coefficient: {}".format(nx.average_clustering(G))))
    elif (clicked2.get()=="None"):
        text1.insert(INSERT,"Select Algorithm!")
        #plt.title('Original Graph')
    
    pos = nx.get_node_attributes(G, 'pos')
    nx.draw(G, pos, with_labels=1, arrows=False, node_color='brown')
    plt.show()
    plt.close()
root.title("Design & Analysis of Algorithms Project")
root.geometry("1024x750")
root.minsize(1024,750)
#root.maxsize(1024,750)
root.configure(bg="maroon")
TITLE= Label(text="Design & Analysis of Algorithms Project",bg="green",fg="white",font=("Bel MT",32,"bold"),borderwidth=3, relief=GROOVE)
TITLE.pack()
f1=Frame(root,bg="grey",borderwidth=3,width=100)
f1.pack()
clicked1=StringVar()
clicked1.set("10")
l1=Label(root,text="No Of Nodes",bg="black",fg="red",font="Ariel")
clicked2=StringVar()
clicked2.set("None")
l2=Label(root,text="Algorithm?",bg="black",fg="red",font="Ariel")
l1.pack()
selectNodes=OptionMenu(f1,clicked1,"10","20","30","40","50","60","70","80","90","100")
selectNodes.pack()
f2=Frame(root,bg="grey",borderwidth=3,width=100)
f2.pack()
selectAlgo=OptionMenu(f2,clicked2,"None","Dijkstra","Prims","Kruskal","FloydWarshal","BellmanFord","Boruvka","Clustering_Coefficient")
selectAlgo.pack()
l2.pack()
bframe=Frame(root,bg="Maroon")
bframe.pack()
preVInput=0
b1=Button(bframe,text="Select",command=setGraph,width=30,bg="Green",fg="white")
b1.pack()
scroll=Scrollbar(root)
scroll.pack(side=RIGHT,fill=Y)
text1=Text(root,width=500,height=26,bd=5,relief="groove",wrap="word",font=('BelMT',12),yscrollcommand=scroll.set)
text1.pack(fill=Y)
scroll.config(command=text1.yview)

root.mainloop()
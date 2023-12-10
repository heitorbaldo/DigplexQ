'''
This code plots the connectivity digraph such that the node positions are  
in the 10-20 system configuration.
'''

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

#fixed positions :
pos = {0: (0,15), 1: (2,16), 2: (4,16), 3: (6, 15), 
       4: (-2,13), 5: (0,12), 6: (2,12), 7: (4,12), 8: (6,12), 9: (8, 13),
       10: (-2,8), 11: (0, 8), 12: (2, 8), 13: (4, 8), 14: (6, 8), 15: (8, 8),
       16: (-2,3), 17: (0, 4), 18: (2, 4), 19: (4, 4), 20: (6, 4), 21: (8, 3),
       22: (2, 0), 23: (4, 0)}


#fixed positions (26 channels):
pos_labels_1020 = {"Fp1": (2,18), "Fp2": (4,18), 
                   "F9": (-2,16.4), "F7": (0,15.4), "F3": (2,15.4), "F4": (4,15.4), "F8": (6,15.4), "F10": (8, 16.4),  
                   "Fc5": (0,12.4), "Fc1": (2,12.4), "Fc2": (4,12.4), "Fc6": (6,12.4),
                   "T3": (-0.7,9.4), "C3": (2, 9.4), "C4": (4, 9.4), "T4": (6.7, 9.4),
                   "Cp5": (0,7), "Cp1": (2, 7), "Cp2": (4, 7), "Cp6": (6, 7),
                   "T5": (-0.7,3), "P3": (2, 4), "P4": (4, 4), "T6": (6.7, 3),
                   "O1": (2, 1.8), "O2": (4, 1.8)}


#fixed positions (24 channels):
pos_labels_1020_extended = {"F9": (0,15), "Fp1": (2,16), "Fp2": (4,16), "F10": (6, 15), 
       "T1": (-2,13), "F7": (0,12), "F3": (2,12), "F4": (4,12), "F8": (6,12), "T2": (8, 13),
       "T9": (-2,8), "T3": (0, 8), "C3": (2, 8), "C4": (4, 8), "T4": (6, 8), "T10": (8, 8),
       "P9": (-2,3), "T5": (0, 4), "P3": (2, 4), "P4": (4, 4), "T6": (6, 4), "P10": (8, 3),
       "O1": (2, 0), "O2": (4, 0)}


#Setting a nxn matrix (n channels):
def set_matrix(M, n):
    A = np.zeros((n,n))
    arr = np.array([i for i in range(len(M))])
    A[arr[:, None], arr] += M
    return A


#Nullifies the principal diagonal
def diag_null(M):
    for i in range(len(M)):
        M[i,i]=0
    return M

def to_binary(M):
    for i in range(len(M)):
        for j in range(len(M)):
            if(M[i, j] > 0.001):
                M[i, j] = 1
            else:
                M[i, j] = 0
    return M

def to_pdc_gct(M, T):
    L = to_binary(T)
    for i in range(len(M)):
        for j in range(len(M)):
            if(L[i, j] > 0):
                L[i, j] = M[i, j]
            else:
                L[i, j] = 0
    return L


#Plotting the weighted digraph (24 nodes):
def plot_digraph_24(M, n, k, color='black', save=False):
        S = np.fliplr(np.flipud(M))
        A = set_matrix(S.transpose(), n)
        B = diag_null(A)
        G = nx.from_numpy_matrix(B, create_using=nx.MultiDiGraph())
        weights = [12*G[u][v][0].get("weight") for u,v in G.edges()]
    
        labels_1020 = {0: "F9", 1: "Fp1", 2: "Fp2", 3: "F10", 4: "T1", 5: "F7", 6: "F3",
                       7: "F4", 8: "F8", 9: "T2", 10: "T9", 11: "T3", 12: "C3", 13: "C4",
                       14: "T4", 15: "T10", 16: "P9", 17: "T5", 18: "P3", 19: "P4", 20: "T6",
                       21: "P10", 22: "O1", 23: "O2"}
        H = nx.relabel_nodes(G, labels_1020)

        plt.figure(3,figsize=(7,7),dpi=300) 
        nx.draw_networkx(H, pos=pos_labels_1020_extended, with_labels=True, node_color="#ffffff", node_size=500,
                connectionstyle='arc3, rad = 0.06', edge_color = color, font_weight='bold', width=weights)
        
        if save == True:
            name = "Digraph_"+str(k)+".png"
            plt.savefig(name, format="PNG")
            plt.show()
        else:
            plt.show()
        
        
#Plotting the weighted digraph (26 nodes):
def plot_digraph_26(M, n, k, color='black', save=False):
        S = np.fliplr(np.flipud(M))
        A = set_matrix(S.transpose(), n)
        B = diag_null(A)
        G = nx.from_numpy_matrix(B, create_using=nx.MultiDiGraph())
        weights = [10*G[u][v][0].get("weight") for u,v in G.edges()]
       
        labels_1020 = {0: "Fp1", 1: "Fp2", 2: "F9", 3: "F7", 4: "F3", 5: "F4", 6: "F8",
                       7: "F10", 8: "Fc5", 9: "Fc1", 10: "Fc2", 11: "Fc6", 12: "T3", 13: "C3",
                       14: "C4", 15: "T4", 16: "Cp5", 17: "Cp1", 18: "Cp2", 19: "Cp6", 20: "T5",
                       21: "P3", 22: "P4", 23: "T6", 24: "O1", 25: "O2"} 
        
        H = nx.relabel_nodes(G, labels_1020)
        plt.figure(3,figsize=(7,7),dpi=300) 
        nx.draw_networkx(H, pos=pos_labels_1020, with_labels=True, node_color="black", node_size=480,
                connectionstyle='arc3, rad = 0.06', edge_color=color, font_color='white', font_size='10', width=weights) 
        
        if save == True:
            name = "Digraph_"+str(k)+".png"
            plt.savefig(name, format="PNG")
            plt.show()
        else:
            plt.show()
        

#Plotting the weighted connectivity digraph with significant GCT:
def plot_digraph_gct(M, L, n, color):
        Q = to_pdc_gct(M, L)
        S = np.fliplr(np.flipud(Q))
        A = set_matrix(S.transpose(), n)
        B = diag_null(A)
        G = nx.from_numpy_matrix(B, create_using=nx.MultiDiGraph())
        weights = [12*G[u][v][0].get("weight") for u,v in G.edges()]
    
        labels_1020 = {0: "F9", 1: "Fp1", 2: "Fp2", 3: "F10", 4: "T1", 5: "F7", 6: "F3",
               7: "F4", 8: "F8", 9: "T2", 10: "T9", 11: "T3", 12: "C3", 13: "C4",
               14: "T4", 15: "T10", 16: "P9", 17: "T5", 18: "P3", 19: "P4", 20: "T6",
               21: "P10", 22: "O1", 23: "O2"}
        H = nx.relabel_nodes(G, labels_1020)
        
        plt.figure(3,figsize=(7,7),dpi=300) 
        nx.draw_networkx(H, pos=pos_labels_1020_extended, with_labels=True, node_color="#ffffff", node_size=500,
            connectionstyle='arc3, rad = 0.06', edge_color = color, width=weights)
        #plt.savefig('plotgraph.png', dpi=300, bbox_inches='tight')
        plt.show()
        
    
    
#Plotting the dynamic weighted connectivity digraph:
#def plot_dynamic_digraph(DnD, time):
from time import sleep
from tqdm import trange, tqdm
from multiprocessing import Pool, RLock, freeze_support
from tqdm import tqdm 
from networkx.algorithms import bipartite
from scipy.stats import binned_statistic as bin1d

import os
import numpy as np 
import networkx as nx 
import pandas as pd 
import matplotlib.pyplot as plt 
import json 
import itertools
import random


def plot_dd(List,label,axes=False):
    
    if type(axes) == bool:

        fig = plt.figure(figsize=(4,4))
        ax = fig.add_subplot(111)

    else:
        ax = axes    

    if min(List) != 0 :

        bins = np.logspace(np.log10(min(List)),np.log10(max(List)),50)
    else:
        bins = np.logspace(np.log10(1),np.log10(max(List)),50)

    pk,k,_ = bin1d(List,List,statistic='count',bins=bins)

    N = k[1:] - k[:-1]

    x,y = k[1:],pk/(sum(pk)*N)

    ax.loglog(x,y,'o',color='lightgray',alpha=0.3)
#         ax.hist(degrees,bins=bins,density=True)
    ax.set_xlabel('k',fontsize=14)
    ax.set_ylabel('P(k)',fontsize=14)
#         ax.set_xscale('log')
#         ax.set_yscale('log')
    ax.set_title(f'{label} Degree Distribution')
    
    return ax 
    

#### Defining class for the network growth model

class ramasco_model:
    
    
    def __init__(self,m,n,n_steps,n_top=0,n_btm=0,init_edges:int=0,init=False):
        
        self.m = m ## number of new actors 
        self.n = n ## total actors added for a every new movie
        self.initialize = init
        if self.n < self.m:
            raise ValueError('n should be larger than m')
        
        self.totalTime = n_steps ## Time Steps for evolution
        self.topNodes = []
        self.btmNodes = []
        self.edges = []
        self.n_edges = np.zeros(self.totalTime)
        
        if self.initialize:
            
            self.init_topNodes = range(n_top)
            self.init_btmNodes = range(n_steps,n_steps+n_btm)
            
            if not init_edges <= self.n*self.m:
                raise ValueError('# of initial edges cannot be more than the product of n*m in unweighted graph')
            
            
            self.topNodes = self.init_topNodes
            self.btmNodes = self.init_btmNodes
            self.edges = random.sample(set(itertools.product(self.init_topNodes,self.init_btmNodes)),init_edges)
        
        return None
    
    def create_Graph(self):
        
        G = nx.Graph()
        G.add_nodes_from(self.topNodes,bipartite=0)
        G.add_nodes_from(self.btmNodes,bipartite=1)
        G.add_edges_from(self.edges)
        return G
    
    def simulate(self):
        
        G = self.create_Graph()
        
        if not self.initialize:
            
            self.topNodes.append(0)
            self.btmNodes.extend(list(range(self.totalTime,self.totalTime+self.n)))
            self.edges.extend(list(zip([0]*self.n,self.btmNodes)))
            
            self.n_edges[0] = len(self.edges)
            potential_old_targets = self.btmNodes
            
            for t in range(1,self.totalTime):
        
                new_btmNodes = list(range(max(self.btmNodes)+1,max(self.btmNodes)+1+self.m))
                new_topNode = t
                
                newActors_edges = list(zip([new_topNode]*self.m,new_btmNodes))
                
                if self.n-self.m > len(potential_old_targets):
                    old_actors = random.sample(potential_old_targets,len(potential_old_targets))
                else:
                    old_actors = random.sample(potential_old_targets,self.n-self.m)
                
                oldActors_edges = list(zip([new_topNode]*len(old_actors),old_actors))
                
                potential_old_targets.extend(old_actors)
                potential_old_targets.extend(new_btmNodes)
                
                self.edges.extend(newActors_edges)
                self.edges.extend(oldActors_edges)
                self.topNodes.append(new_topNode)
                self.btmNodes.extend(new_btmNodes)
                
                self.n_edges[t] = len(self.edges)
            
            G.add_nodes_from(self.topNodes,bipartite=0)
            G.add_nodes_from(self.btmNodes,bipartite=1)
            G.add_edges_from(self.edges)
            
            return G
        
        else:
            print('Not programmed yet')
            return None
            
    
    def plot_dd(self,Graph,set_indx,label,axes=False):
        
        if type(axes) == bool:
        
            fig = plt.figure(figsize=(4,4))
            ax = fig.add_subplot(111)
        
        else:
            ax = axes
        
        if type(set_indx) == int:
        
            nodeSet = [i for i,j in Graph.nodes(data=True) if j['bipartite']==set_indx]
            degrees = [Graph.degree(i) for i in nodeSet]
        
        else:
            
            degrees = list(dict(Graph.degree()).values())
        
        degrees.sort()
        x,y = np.unique(degrees,return_counts=True)
        
        if min(degrees) != 0 :
            
            bins = np.logspace(np.log10(min(degrees)),np.log10(max(degrees)),20)
        else:
            bins = np.logspace(np.log10(1),np.log10(max(degrees)),20)
        
        pk,k,_ = bin1d(degrees,degrees,statistic='count',bins=bins)
        
        
        N = k[1:] - k[:-1]
        
        x,y = k[1:],pk/(sum(pk)*N)

#         y = y/sum(y)
        
        ax.plot(x,y,'o',color='darkred')
#         ax.hist(degrees,bins=bins,density=True)
        ax.set_xlabel('k',fontsize=14)
        ax.set_ylabel('P(k)',fontsize=14)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_title(f'{label} Degree Distribution')
        
        return None
    
    
    def plot_degcorr(self,Graph,label,axes=False):
        
        if type(axes) == bool:
        
            fig = plt.figure(figsize=(4,4))
            ax = fig.add_subplot(111)
        
        else:
            ax = axes
        
        x,y = [],[]
        
        
        for (i,j) in Graph.edges():
            
            x.append(Graph.degree(i))
            y.append(Graph.degree(j))
        
        k_nn,bins,_, = bin1d(x,y,statistic='mean',bins=30) 
        
        ax.loglog(bins[1:],k_nn,'--.',color='darkred')
        ax.set_xlabel('k',fontsize=14)
        ax.set_ylabel(r'$k_{nn}(k)$',fontsize=18)
        ax.set_title('Near Nbr Degree dist.')

        return None


class user_object_model(ramasco_model):
    
    
    def __init__(self,w,m,n,n_steps,gamma_e,gamma_i,init=False):
        
        self.w = w
        self.m = m
        self.n = n
        
#         if (self.w < self.m) or (self.w < self.n):
            
#             raise ValueError('w should be greater than m and n') # here we assume that user's attention is finite 
        
        self.gamma_e = gamma_e
        self.gamma_i = gamma_i
        self.totalTime = n_steps
        
        self.topNodes = []   ### users
        self.btmNodes = []   ### objects
        self.edges = []
        self.initialize = init
        self.n_edges = np.zeros(self.totalTime)
        
        return None
    
    def get_prob(self,nodeList,gamma):
        
        '''
        this gives the probability for selection of every node based on the attachment kernel. 
                    A(k) = (k+gamma)/sum(k+gamma)
                    
        Note that if gamma --> infty then it is a purely random attachment whereas when gamma --> 0 it
        is purely preferential attachment 
        
        '''
        
        x,px = np.unique(nodeList,return_counts=True)
        px = (px+gamma)/np.sum(px+gamma)
        
        return x,px
    
    
    def simulate(self):
        
        G = self.create_Graph()
        
        if not self.initialize:
            
            self.topNodes.append(0)
            self.btmNodes.extend(list(range(self.totalTime,self.totalTime+self.w)))
            
            selected_targets = []
            
            potential_targets = selected_targets + self.btmNodes
            
            if self.w < self.m: ### if number of new objects is less than new edges by a user we set m=w
                
                self.edges.extend(list(zip([0]*self.w),potential_targets))
                selected_targets.extend(self.btmNodes)
            else:
                
                NewUsr_targets = np.random.choice(potential_targets,self.m,replace=False)
                
                self.edges.extend(list(zip([0]*self.m,NewUsr_targets)))
                
                selected_targets.extend(NewUsr_targets)
            
            self.n_edges[0] = len(self.edges)
            
#             for t in tqdm(range(1,self.totalTime)):
            for t in range(1,self.totalTime):
                self.topNodes.append(t)
                new_btmNodes = list(range(max(self.btmNodes)+1,max(self.btmNodes)+1+self.w))
                self.btmNodes.extend(new_btmNodes)
                
                potential_targets = selected_targets + self.btmNodes
                
                if (len(self.btmNodes) < self.m) or (len(self.btmNodes) < self.n): 
                    ### Change this for new edges and old edges
                    self.edges.extend(list(zip([t]*self.m),self.btmNodes))
                
                else:
                    
                    potential_NewUsr_targets,prob_e = self.get_prob(potential_targets,self.gamma_e)
                    
                    NewUsr_targets = np.random.choice(potential_NewUsr_targets,self.m,replace=False,p=prob_e)
                    
                    potential_OldUsr_targets,prob_i = self.get_prob(potential_targets,self.gamma_i)
                    
                    OldUsr_targets = np.random.choice(potential_OldUsr_targets,self.n,replace=False,p=prob_i)
                    
                    randomOldusr = np.random.choice(self.topNodes[:-1],self.n)
                    
                    self.edges.extend(list(zip(randomOldusr,OldUsr_targets)))
                    self.edges.extend(list(zip([t]*self.m,NewUsr_targets)))
                    
                    selected_targets.extend(NewUsr_targets)
                    selected_targets.extend(OldUsr_targets)
                
                    
                self.n_edges[t] = len(self.edges)
            
            G.add_nodes_from(self.topNodes,bipartite=0)
            G.add_nodes_from(self.btmNodes,bipartite=1)
            G.add_edges_from(self.edges)
            
            return G
        else:
            print("not programmed yet")
            return None 
    
    
    def analytical_prob(k,r,s,k0):
        
        return (r/s)*((k0)**(r/s))*(1/(k+k0)**(1+(r/s)))
    
    
    
def worker(procnum):
    """Function run by worker processes"""
#     print(f"Worker {procnum} running in process {os.getpid()}")
    
    return procnum*procnum

def multi_sim(i):
    
    G = ramasco_model(1,3,1000,n_top=2,n_btm=3,init_edges=2,init=False)
    G1 = G.simulate()
    
    btmNodes = [n for n,j in G1.nodes(data=True) if j['bipartite']==1]
    btmNodes.sort()
    degrees = [G1.degree(i) for i in btmNodes]
    return degrees

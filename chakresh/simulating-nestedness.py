from essentials import *

def get_nestedness(BiMat):
    
    '''
    takes the biAdjacency matrix of the bipartite graph and
    gives the nestedness calculated from the R script 
    '''
    
    np.savetxt('graph_matrices/test.csv',BiMat,delimiter=',')
    
    import subprocess

    cmd = r'''cd /Users/chakreshsingh/Documents/bipartite_growth_models/ | Rscript mat_temp.R'''

    x = float(subprocess.check_output(cmd, shell=True))
    
    return x


gamma_e = np.random.randint(1,50,30)
gamma_i = np.random.randint(1,50,30)

gamma_e = list(map(int,sorted(gamma_e)))
gamma_i = list(map(int,sorted(gamma_i)))



w,m,n = 5,3,2

nt = 514

Z = []

for i in tqdm(gamma_e):
    
    for j in gamma_i:
        
        model = user_object_model(w,m,n,nt,i,j)
        G = model.simulate()
        B = [n for n,j in G.nodes(data=True) if j['bipartite']==1]
        T = [n for n,j in G.nodes(data=True) if j['bipartite']==0]

        BiMat = bipartite.biadjacency_matrix(G,row_order=B,column_order=T).toarray()

        Z.append(get_nestedness(BiMat))

        
N = {'gamma_e':list(gamma_e),'gamma_i':list(gamma_i),'Z':Z}
        
with open('simulation_nestedness_v2.json','w') as f:
    
    json.dump(N,f,indent=1)
    
print('-----------------------------DONE-------------------------------')
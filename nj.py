import numpy as np

'''
Source:
https://en.wikipedia.org/wiki/Neighbor_joining
https://www.tenderisthebyte.com/blog/2022/08/31/neighbor-joining-trees/
https://moodle.ufrgs.br/pluginfile.php/5635441/mod_resource/content/0/Phylogenetic_inference.pdf
https://medium.com/geekculture/phylogenetic-trees-implement-in-python-3f9df96c0c32
'''

def neighbor_joining(distance_matrix) -> np.matrix:
    d = np.copy(distance_matrix)
    n = d.shape[0]

    for _ in range(n-2):
        n = d.shape[0]
        q_matrix = np.zeros((n,n))
        
        for i in range(n):
            for j in range(n):
                if i != j:
                    q_matrix[i,j] = (n-2)*d[i,j] - d[i, 0:].sum() - d[0:, j].sum()
        
        i, j = np.unravel_index(np.argmin(q_matrix), q_matrix.shape)

        print(f"{q_matrix=}")
        print(f"{i=} {j=}")
        
        delta_ij = (d[i, 0:].sum() + d[0:, j].sum()) / (n-2)
        branch_length_i = 1/2 * (d[i,j] + delta_ij)
        branch_length_j = 1/2 * (d[i,j] - delta_ij)
        
        '''
        Nova matriz removendo as linha/coluna i e j de D 
        e acrescentando uma linha/coluna m tal que para qualquer
        k -> Dk,m = (Dk,i + Dk,j - Di,j)/2
        '''

        #dprime = np.copy(d) 
        #dprime = np.delete(dprime, i, 0)
        #dprime = np.delete(dprime, i, 1)
        #dprime = np.delete(dprime, j, 0)
        #dprime = np.delete(dprime, j, 1)

        dprime = np.zeros((n-1, n-1))
        for k in range(n-1):
            dprime[i, k] = (d[k+1,i] + d[k+1,j] - d[i,j])/2
            dprime[k, i] = (d[k+1,i] + d[k+1,j] - d[i,j])/2

        d = np.copy(dprime)
        
        print(f"{d=}")

'''
 	L. braziliensis 	T. rangeli 	T. cruzi 	T. gambiae
L. braziliensis 	0.000 	0.010 	0.300 	0.280
T. rangeli 	0.010 	0.000 	0.280 	0.270
T. cruzi 	0.300 	0.280 	0.000 	0.015
T. gambiae 	0.280 	0.270 	0.015 	0.000
'''

distance_matrix = np.matrix([[0.000, 0.010, 0.300, 0.280], 
               [0.010, 0.000, 0.280, 0.270],
               [0.300, 0.280, 0.000, 0.015], 
               [0.280, 0.270, 0.015, 0.000]])

CONVERSION_DICT = {
    0: 'L. braziliensis',
    1: 'T. rangeli',
    2: 'T. cruzi',
    3: 'T. gambiae'
}

d2 = np.matrix([[0, 5, 9, 9, 8],
                [5, 0, 10, 10, 9],
                [9, 10, 0, 8, 7],
                [9, 10, 8, 0, 3],
                [8, 9, 7, 3, 0]])

d3 = np.matrix([[0, 13, 21, 22],
                [13, 0, 12, 13],
                [21, 12, 0, 13],
                [22, 13, 13, 0]])

neighbor_joining(d3)
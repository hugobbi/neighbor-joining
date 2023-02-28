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
    q_matrix = np.zeros((n,n))

    for k in range(n-2):
        for i in range(n):
            for j in range(n):
                if i != j:
                    q_matrix[i,j] = (n-2)*d[i,j] - d[i, 0:].sum() - d[0:, j].sum()
        
        i, j = np.unravel_index(np.argmin(q_matrix), q_matrix.shape)

        #dij = (D[i][j] + (np.sum(D[i]) - np.sum(D[j])) / (n-2)) / 2
        
        #dist_i_new_node = 1/2*d[i,j] + (1/(2*(n-2)))*(d.sum[i, 1:]-d.sum[1:, j])
        
        print(d)
        print(q_matrix)

        print(f"{i=} {j=}")

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

neighbor_joining(d2)
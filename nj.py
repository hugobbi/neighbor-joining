import numpy as np

'''
Source:
https://en.wikipedia.org/wiki/Neighbor_joining
https://www.tenderisthebyte.com/blog/2022/08/31/neighbor-joining-trees/
https://moodle.ufrgs.br/pluginfile.php/5635441/mod_resource/content/0/Phylogenetic_inference.pdf
https://medium.com/geekculture/phylogenetic-trees-implement-in-python-3f9df96c0c32
'''

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

def sum_array(matrix, n, i):
    sum = 0
    for k in range(1, n):
        sum += matrix(i, k)

def neighbor_joining(d) -> np.matrix:
    n = d.shape[0]
    q_matrix = np.zeros((n,n))

    for i in range(n):
        for j in range(n):
            q_matrix[i,j] = (n-2)*d[i,j] - d[i, 1:].sum() - d[1:, j].sum()

    print(d)
    print(q_matrix)


neighbor_joining(distance_matrix)
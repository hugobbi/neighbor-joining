import numpy as np
from Bio import Phylo

'''
Source:
https://en.wikipedia.org/wiki/Neighbor_joining
https://www.tenderisthebyte.com/blog/2022/08/31/neighbor-joining-trees/
https://moodle.ufrgs.br/pluginfile.php/5635441/mod_resource/content/0/Phylogenetic_inference.pdf
https://medium.com/geekculture/phylogenetic-trees-implement-in-python-3f9df96c0c32
'''

def neighbor_joining(distance_matrix, conversion_dict) -> np.matrix:
    d = np.copy(distance_matrix)
    n = d.shape[0]

    tree_labels = [conversion_dict[str(i)] for i in range(n)]
    print(tree_labels)

    while n > 2:
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

        dprime = np.copy(d) 
        dprime = np.delete(dprime, j, 0)
        dprime = np.delete(dprime, j, 1)
        remaining_idx = [r_idx for r_idx in range(n) if r_idx != i]
        for k, r_idx in enumerate(remaining_idx): # i row/column will become new node row/column
            new_distance = (d[r_idx,i] + d[r_idx,j] - d[i,j])/2
            dprime[i, k] = new_distance
            dprime[k, i] = new_distance
            if k == i:
                dprime[i, k] = 0
                dprime[k, i] = 0

        d = np.copy(dprime)
        n = n-1
        print(f"{d=}")

        new_node = f"({tree_labels[i]}, {tree_labels[j]})"
        new_tree_labels = [tree_label for tree_label in tree_labels if tree_label not in {tree_labels[i], tree_labels[j]}]
        new_tree_labels.append(new_node)
        tree_labels = new_tree_labels.copy()

    final_label = f"({tree_labels[0]}, {tree_labels[1]});"
    return final_label

d1 = np.matrix([[0.0000000, 1.4230670, 	1.2464736, 	2.4829678, 	1.6023069, 	0.7448415],
                [1.4230670,	0.0000000,  1.8985344, 	3.1350286, 	2.2543677, 	1.9430337],
                [1.2464736,	1.8985344, 	0.0000000, 	1.4820114, 	0.6013505, 	1.7664403],
                [2.4829678,	3.1350286, 	1.4820114, 	0.0000000, 	1.2678046, 	3.0029345],
                [1.6023069,	2.2543677, 	0.6013505, 	1.2678046, 	0.0000000, 	2.1222736],
                [0.7448415,	1.9430337, 	1.7664403, 	3.0029345, 	2.1222736, 	0.0000000]])

d2 = np.matrix([[0, 5, 9, 9, 8],
                [5, 0, 10, 10, 9],
                [9, 10, 0, 8, 7],
                [9, 10, 8, 0, 3],
                [8, 9, 7, 3, 0]])

d3 = np.matrix([[0, 13, 21, 22],
                [13, 0, 12, 13],
                [21, 12, 0, 13],
                [22, 13, 13, 0]])

d4 = np.matrix([[0.000, 0.010, 0.300, 0.280], 
                [0.010, 0.000, 0.280, 0.270],
                [0.300, 0.280, 0.000, 0.015], 
                [0.280, 0.270, 0.015, 0.000]])

CONVERSION_DICT = {
    "0": "A",
    "1": "B",
    "2": "C",
    "3": "D",
    "4": "E",
    "5": "F"
}

tree = neighbor_joining(d2, CONVERSION_DICT)
print(tree)
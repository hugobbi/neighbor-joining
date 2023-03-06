import numpy as np
from newick import loads

def neighbor_joining(d, ids):
    while d.shape[0] > 2:
        q = np.zeros(d.shape)
        n = d.shape[0]
        # calculate q_matrix
        for i in range(n):
            for j in range(n):
                if i != j:
                    q[i,j] = (n-2)*d[i,j] - d[i].sum() - d[j].sum()

        i, j = np.unravel_index(np.argmin(q), q.shape)

        # calculate distances from new node and its decandants
        dist_i_to_u = (1/2 * d[i, j]) + ((d[i].sum() - d[j].sum()) / (2*(n-2)))
        dist_j_to_u = d[i, j] - dist_i_to_u
        new_node = f"({ids[i]}: {round(dist_i_to_u, 4)}, {ids[j]}: {round(dist_j_to_u, 4)})"
        new_ids = [new_node]
        for id in ids:
            if id not in {ids[i], ids[j]}:
                new_ids.append(id)

        # calculates new d matrix
        dprime = np.zeros((n-1, n-1))
        ij_indexes = [i, j]
        dprime[1:, 1:] = np.delete(np.delete(d, ij_indexes, axis=1), 
                            ij_indexes, axis=0)
        
        # calculate distances for new node
        dist_k_to_u = 1/2 * (d[i] + d[j] - d[i, j]) 
        dist_k_to_u = np.delete(dist_k_to_u, ij_indexes) # remove joined taxons
        dist_k_to_u = np.concatenate([[0], dist_k_to_u])
        
        for k in range(n-1): # assign new distances for the new node
           dprime[0, k] = dprime[k, 0] = dist_k_to_u[k]

        d = np.copy(dprime)
        ids = new_ids.copy()

        print(f"{i} {j}")

    final_ids = f"{ids[0], ids[1]}"
    return final_ids

d1 = np.array([[0.0000000, 1.4230670, 	1.2464736, 	2.4829678, 	1.6023069, 	0.7448415],
                [1.4230670,	0.0000000,  1.8985344, 	3.1350286, 	2.2543677, 	1.9430337],
                [1.2464736,	1.8985344, 	0.0000000, 	1.4820114, 	0.6013505, 	1.7664403],
                [2.4829678,	3.1350286, 	1.4820114, 	0.0000000, 	1.2678046, 	3.0029345],
                [1.6023069,	2.2543677, 	0.6013505, 	1.2678046, 	0.0000000, 	2.1222736],
                [0.7448415,	1.9430337, 	1.7664403, 	3.0029345, 	2.1222736, 	0.0000000]])

CONVERSION_DICT = {
    "0": "A",
    "1": "B",
    "2": "C",
    "3": "D",
    "4": "E",
    "5": "F",
    "6": "G",
    "7": "H"
}

ids = ["A", "B", "C",  "D",  "E",  "F"]

tree = neighbor_joining(d1, ids)
print(tree)
print(loads(tree)[0].ascii_art())
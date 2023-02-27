import numpy as np

'''
 	L. braziliensis 	T. rangeli 	T. cruzi 	T. gambiae
L. braziliensis 	0.000 	0.010 	0.300 	0.280
T. rangeli 	0.010 	0.000 	0.280 	0.270
T. cruzi 	0.300 	0.280 	0.000 	0.015
T. gambiae 	0.280 	0.270 	0.015 	0.000
'''

D = np.matrix([[0.000, 0.010, 0.300, 0.280], 
               [0.010, 0.000, 0.280, 0.270],
               [0.300, 0.280, 0.000, 0.015], 
               [0.280, 0.270, 0.015, 0.000]])

CONVERSION_DICT = {
    0: 'L. braziliensis',
    1: 'T. rangeli',
    2: 'T. cruzi',
    3: 'T. gambiae'
}

def neighbor_joining(distance_matrix) -> np.matrix:
    print(distance_matrix)


neighbor_joining(D)
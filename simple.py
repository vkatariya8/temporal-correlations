import numpy as np
import math

zero_matrix = [[1, 0, 0, 0], [0, 0, 0.5, 0], [0, 0.5, 0, 0], [0, 0, 0, 0]]
zero_matrix = np.asarray(zero_matrix)
one_matrix = [[0, 0, 0, 0], [0, 0, 0.5, 0], [0, 0.5, 0, 0], [0, 0, 0, 1]]
one_matrix = np.asarray(one_matrix)
alpha = math.sqrt(0.8)
beta = math.sqrt(0.2)

result_matrix = alpha**2 * zero_matrix + beta**2 * one_matrix

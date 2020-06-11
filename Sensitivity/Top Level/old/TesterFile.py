import itertools as it
import numpy as np

#a  =  list(set(list(it.permutations([1,1,1,1,1,-1,-1,-1,-1,-1,0,0,0,0,0],5))))
a  =  (list(it.permutations([1,1,1,1,1,-1,-1,-1,-1,-1,0,0,0,0,0,1,1,1,1,1,-1,-1,-1,-1,-1,0,0,0,0,0],5)))
print(a)

# for i,order in enumerate(list(set(list(it.permutations([1,1,1,1,1,-1,-1,-1,-1,-1,0,0,0,0,0],5))))):
#     S = np.zeros([5, 3])
#     for l in range(3):
#         for q in range(5):
#             S[q,l] = S[q,l]+np.array(order)[q]
#     print(S)
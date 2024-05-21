#%%
from sklearn.datasets import load_iris
import numpy as np
import archetypes as arch


wine=load_iris()
X = wine.data
target=wine.target
aa_kwargs = {
    "n_archetypes": 4,
    "n_init": 5,
    "max_iter": 200,
    "verbose":False,
    "tol":1e-4,
    
}
model = arch.AA(**aa_kwargs,algorithm_init="furthest_sum")
model.fit(X)

model.archetypes_

#%% 
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 10))
arch.simplex(model.alphas_, c=target, alpha=0.5, show_circle=False, show_direction=True)

plt.show()
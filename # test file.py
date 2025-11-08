# test file

import numpy as np
import matplotlib.pyplot as plt
import random

def randomwalk3D(n):
    x, y, z = np.zeros(n), np.zeros(n), np.zeros(n) # de

    directions = ["up", "down", "left", "right", "in", "out"] # defining the direction which can be taken in 3d space
    for i in range(1, n):
        step = random.choice(directions)  # randomizing steps taken 
        
        if step == "right":
            x[i] = x[i - 1] + 1
            y[i] = y[i - 1] 
            z[i] = z[i - 1] 
        elif step == "left":
            x[i] = x[i - 1] - 1
            y[i] = y[i - 1] 
            z[i] = z[i - 1]
        elif step == "up":
            x[i] = x[i - 1] 
            y[i] = y[i - 1] + 1
            z[i] = z[i - 1]  
        elif step == "down":
            x[i] = x[i - 1] 
            y[i] = y[i - 1] - 1 
            z[i] = z[i - 1] 
        elif step == "in":
            x[i] = x[i - 1] 
            y[i] = y[i - 1] 
            z[i] = z[i - 1] - 1 
        elif step == "out":
            x[i] = x[i - 1] 
            y[i] = y[i - 1] 
            z[i] = z[i - 1] + 1

    return x, y, z

# parameters
n_walkers = 5

x_data, y_data, z_data = randomwalk3D(100)

ax = plt.subplot(1, 1, 1, projection='3d')
ax.plot(x_data, y_data, z_data, alpha=0.9)
ax.scatter(x_data[-1], y_data[-1], z_data[-1])

plt.show
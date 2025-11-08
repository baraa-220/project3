#Task 2 topic 2
'To make a fast function of random walkers in 3D space we can implement the same method as above but optimize it for speed.'

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random

def randomwalk3D(n_walkers, n_steps):
    x, y, z = np.zeros(n_steps), np.zeros(n_steps), np.zeros(n_steps) # defining 3D arrays

    directions = ["up", "down", "left", "right", "in", "out"] # defining the direction which can be taken in 3d space
    for i in range(1, n_steps):
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

n_walkers = 3  # amount of walkers
n_steps = 100  # amount of steps

fig = plt.figure()  
ax = fig.add_subplot(111, projection='3d')  # Creating a 3D space and a figure in 3D

for _ in range(n_walkers):   # Making a loop for the walkers
    x, y, z = randomwalk3D(1, n_steps) 
    ax.plot(x,y,z)
ax.set_title('random walks in 3D')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.show()

#%%

def randomwalk3D(n_walkers, n_steps):
    x, y, z = np.zeros(n_steps), np.zeros(n_steps), np.zeros(n_steps) # defining 3D arrays

    directions = ["up", "down", "left", "right", "in", "out"] # defining the direction which can be taken in 3d space
    for i in range(1, n_steps):
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

n_walkers = 3  # amount of walkers
n_steps = 100  # amount of steps

fig = plt.figure()  
ax = fig.add_subplot(111, projection='3d')  # Creating a 3D space and a figure in 3D

for _ in range(n_walkers):   # Making a loop for the walkers
    ax.plot(x,y,z)
ax.set_title('random walks in 3D')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.show()

#%%


def DNArandomwalk3D(n_steps):
    hydrogen = np.zeros(n_steps)
    oxygen = np.zeros(n_steps)
    carbon = np.zeros(n_steps)

    directions = ["up", "down", "left", "right", "in", "out"] # defining the direction which can be taken in 3d space
    for i in range(1, n_steps):
        step = random.choice(directions)  # randomizing steps taken 
        
        if step == "right":
            hydrogen[i] = hydrogen[i - 1] + 1
            oxygen[i] = oxygen[i - 1] 
            carbon[i] = carbon[i - 1] 
        elif step == "left":
            hydrogen[i] = hydrogen[i - 1] - 1
            oxygen[i] = oxygen[i - 1] 
            carbon[i] = carbon[i - 1]
        elif step == "up":
            hydrogen[i] = hydrogen[i - 1] 
            oxygen[i] = oxygen[i - 1] + 1
            carbon[i] = carbon[i - 1]  
        elif step == "down":
            hydrogen[i] = hydrogen[i - 1] 
            oxygen[i] = oxygen[i - 1] - 1 
            carbon[i] = carbon[i - 1] 
        elif step == "in":
            hydrogen[i] = hydrogen[i - 1] 
            oxygen[i] = oxygen[i - 1] 
            carbon[i] = carbon[i - 1] - 1 
        elif step == "out":
            hydrogen[i] = hydrogen[i - 1] 
            oxygen[i] = oxygen[i - 1] 
            carbon[i] = carbon[i - 1] + 1

    return hydrogen, oxygen, carbon

def accessible_volume(hydrogen, oxygen, carbon, grid_size=50): # Creating a 3d grid for the volume
    grid = np.zeros((grid_size, grid_size, grid_size))
    coordinates = np.vstack((hydrogen, oxygen, carbon)).T + grid_size / 2 # keep coordinates as positive values
    coordinates = np.clip(coordinates, 0, grid_size - 1) # keep the coordinates inside the grid

    coordinates_int = coordinates.astype(int)
    grid[coordinates_int[:,0], coordinates_int[:,1], coordinates_int[:,2]] = 1
    filled = np.sum(grid)
    total = grid_size ** 3
    return total - filled

    
n_simulations = 5
n_steps = 50


for i in range(n_simulations):
    hydrogen, oxygen, carbon = DNArandomwalk3D(n_steps)
    vol = accessible_volume(hydrogen, oxygen, carbon)
    

print("Accesible volume: ", int(vol))
#%%

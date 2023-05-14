import math
import random

N = 30      # number of UE (10, 20, or 30)
R = 500     # cell radius
X = []      # legal x coordinates
Y = []      # legal x coordinates
i = 0       # looping index

while i < N:
    # generate a point in the square
    x = random.uniform(-R, R)
    y = random.uniform(-R, R)
    # check if out of circle
    r = math.sqrt(x**2 + y**2)
    if r < R and x != 0:
        theta = math.atan2(y, x)        # range in [-pi, pi]
        # fit theta into the 1st sector of the 6
        while theta < 0:
            theta += math.pi / 3
        while theta > math.pi / 3:
            theta -= math.pi / 3
        xp = r * math.cos(theta)
        yp = r * math.sin(theta)
        # check if out of hexagon
        if yp < math.sqrt(3) * (R - xp):
            X.append(x)
            Y.append(y)
            i += 1
print("X =", X, ";")
print("Y =", Y, ";")

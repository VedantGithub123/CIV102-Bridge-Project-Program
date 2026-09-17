# Libraries
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from copy import deepcopy
import numpy as np
from math import pi
from functools import lru_cache

# Constants for material properties
matboardTensileStrength = 30
matboardCompressiveStrength = 6
matboardShearStrength = 4
matboardYoungsModulus = 4000
matboardPoissons = 0.2
cementShearStrength = 2

# Constants for bridge geometry
bridgeLength = 1200
matboardThickness = 1.27

# Load case
case1 = False

# Parameters to describe the geometry of the bridgs
# Rectangles are described with two tuples: [(topLeftX, topLeftY), (bottomRightX, bottomRightY)]
# Cross section is described as a list of rectangles: [rect1, rect2, ...]
crossSectionKeyPoints = [
    [[[(0, 0), (100, -1.27)], [(10, -1.27), (10+1.27+5, -2*1.27)], [(90-1.27-5, -1.27), (90, -2*1.27)], [(10, -1.27*2), (10+1.27, -1.27-75)], [(90-1.27, -1.27*2), (90, -1.27-75)]], 0],
    [[[(0, 0), (100, -1.27)], [(10, -1.27), (10+1.27+5, -2*1.27)], [(90-1.27-5, -1.27), (90, -2*1.27)], [(10, -1.27*2), (10+1.27, -1.27-75)], [(90-1.27, -1.27*2), (90, -1.27-75)]], 1200]
] # [[cross section, position]] ...]

glueKeyPoints = [
    [[-1.27, 6.27*2, 0], [-1.27, 6.27*2, 1200]], # First Glue Tab
    # 2nd Glue Tab
    # ...
] # [[[y-value, width, position]...] ...]

diaphragmPositions = [
    0,
    400,
    800,
    1200
    # ..
]

type1FailureKeyPoints = [
    [[0, 1.27, 80, 0], [0, 1.27, 80, 1200]], # First Place
    # ...
] # [[[y-value, thickness, width, position]...] ...]

type2FailureKeyPoints = [
    [[0, 1.27, 10, 0], [0, 1.27, 10, 1200]], # First Place
    [[0, 1.27, 10, 0], [0, 1.27, 10, 1200]], # 2nd Place
    # ...
] # [[[y-value, thickness, width, position]...] ...]

type3FailureKeyPoints = [
    [[-1.27, 1.27, 0], [-1.27, 1.27, 1200]], # First Place
    # 2nd Place
    # ...
] # [[[y-value, thickness, position]...] ...]

type4FailureKeyPoints = [
    [[75, 1.27, 0], [75, 1.27, 1200]], # First Place
    [[75, 1.27, 0], [75, 1.27, 1200]], # 2nd Place
    # ...
] # [[[height, thickness, position]...] ...]
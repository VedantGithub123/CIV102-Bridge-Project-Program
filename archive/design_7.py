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
    [[[(0, 0), (100, -1.27)], [(12.5+1.27, -1.27), (100-12.5-1.27, -1.27*2)], [(7.5, -1.27), (12.5+1.27, -2*1.27)], [(100-12.5-1.27, -1.27), (100-7.5, -2*1.27)], [(12.5, -1.27*2), (12.5+1.27, -1.27-125)], [(100-12.5-1.27, -1.27*2), (100-12.5, -1.27-125)]], 0],
    [[[(0, 0), (100, -1.27)], [(12.5+1.27, -1.27), (100-12.5-1.27, -1.27*2)], [(7.5, -1.27), (12.5+1.27, -2*1.27)], [(100-12.5-1.27, -1.27), (100-7.5, -2*1.27)], [(12.5, -1.27*2), (12.5+1.27, -1.27-125)], [(100-12.5-1.27, -1.27*2), (100-12.5, -1.27-125)]], 250],
    [[[(0, 0), (100, -1.27)], [(12.5+1.27, -1.27), (100-12.5-1.27, -1.27*2)], [(7.5, -1.27), (12.5+1.27, -2*1.27)], [(100-12.5-1.27, -1.27), (100-7.5, -2*1.27)], [(12.5, -1.27*2), (12.5+1.27, -1.27-150)], [(100-12.5-1.27, -1.27*2), (100-12.5, -1.27-150)]], 420],
    [[[(0, 0), (100, -1.27)], [(12.5+1.27, -1.27), (100-12.5-1.27, -1.27*2)], [(7.5, -1.27), (12.5+1.27, -2*1.27)], [(100-12.5-1.27, -1.27), (100-7.5, -2*1.27)], [(12.5, -1.27*2), (12.5+1.27, -1.27-150)], [(100-12.5-1.27, -1.27*2), (100-12.5, -1.27-150)]], 780],
    [[[(0, 0), (100, -1.27)], [(12.5+1.27, -1.27), (100-12.5-1.27, -1.27*2)], [(7.5, -1.27), (12.5+1.27, -2*1.27)], [(100-12.5-1.27, -1.27), (100-7.5, -2*1.27)], [(12.5, -1.27*2), (12.5+1.27, -1.27-125)], [(100-12.5-1.27, -1.27*2), (100-12.5, -1.27-125)]], 950],
    [[[(0, 0), (100, -1.27)], [(12.5+1.27, -1.27), (100-12.5-1.27, -1.27*2)], [(7.5, -1.27), (12.5+1.27, -2*1.27)], [(100-12.5-1.27, -1.27), (100-7.5, -2*1.27)], [(12.5, -1.27*2), (12.5+1.27, -1.27-125)], [(100-12.5-1.27, -1.27*2), (100-12.5, -1.27-125)]], 1200],
] # [[cross section, position]] ...]

glueKeyPoints = [
    [[-1.27, 85, 0], [-1.27, 85, 1200]], # First Glue Tab
    # 2nd Glue Tab
    # ...
] # [[[y-value, width, position]...] ...]   

diaphragmPositions = [
    0,
    120,
    260,
    420,
    780,
    940,
    1080,
    1200
    # ..
]

type1FailureKeyPoints = [
    [[0, 1.27*2, 75, 0], [0, 1.27*2, 75, 1200]], # First Place
    # ...
] # [[[y-value, thickness, width, position]...] ...]

type2FailureKeyPoints = [
    [[0, 1.27, 12.5, 0], [0, 1.27, 12.5, 1200]], # First Place
    [[0, 1.27, 12.5, 0], [0, 1.27, 12.5, 1200]], # 2nd Place
    # ...
] # [[[y-value, thickness, width, position]...] ...]

type3FailureKeyPoints = [
    [[-1.27, 1.27, 0], [-1.27, 1.27, 1200]], # First Place
    # 2nd Place
    # ...
] # [[[y-value, thickness, position]...] ...]

type4FailureKeyPoints = [
    [[125, 1.27, 0], [125, 1.27, 250], [150, 1.27, 420], [150, 1.27, 780], [125, 1.27, 950], [125, 1.27, 1200]], # First Place
    [[125, 1.27, 0], [125, 1.27, 250], [150, 1.27, 420], [150, 1.27, 780], [125, 1.27, 950], [125, 1.27, 1200]], # First Place
    # ...
] # [[[height, thickness, position]...] ...]
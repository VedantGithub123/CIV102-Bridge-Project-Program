# Program to analyze stresses in a given geometry, will determine modes of failure and the maximum load
# This program is designed to account for varying cross-sections, increasing its complexity
# Team 507
# Bridge Iteration 1
# All units in newtons and millimeters
# Assumptions: Diaphragms are infinitely thin
#              Glue tabs are negligible for local buckling
#              For shear buckling, side plates are assumed to be rectangles 
#              Neglect self-weight of the bridge
# Notes: This program is designed to calculate safety factors by comparing the maximum allowable Shear Force
#        and/or Bending Moment with the Shear Force and/or Bending Moment the beam faces. As such, it does not
#        compare the stresses to find the factors of safety
# 
#        Creating an envelope is done by simulating the train at every position on the bridge and creating an SFD/BMD
#        from that. Then the maximum at each point is taken.
# 
#        The maximum allowable forces/moments is done by calculating the allowable forces/moments
#        at every point along the bridge.
# 
#        The safety factors are calculated by taking the ratio of the encountered Shear Force and/or Bending Moment
#        to the calculated safe Shear Force and/or Bending Moment at every point on the bridge and finding the minimum

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


############################################################################################################

# GEOMETRIC PROPERTIES CALCULATIONS

# Getting the cross section based on the position
def getCrossSection(position):
    crossSection = []
    
    # Interpolates between the key points
    for i in range(len(crossSectionKeyPoints)):
        if position == crossSectionKeyPoints[i][1]:
            crossSection = deepcopy(crossSectionKeyPoints[i][0])
            break

        if crossSectionKeyPoints[i][1] > position:
            for rect1, rect2 in zip(crossSectionKeyPoints[i-1][0], crossSectionKeyPoints[i][0]):
                rect = []
                for point1, point2 in zip(rect1, rect2):
                    newPoint = []
                    for val1, val2, in zip(point1, point2):
                        newPoint.append(
                            (val2-val1)/(crossSectionKeyPoints[i][1]-crossSectionKeyPoints[i-1][1])*(position-crossSectionKeyPoints[i-1][1]) + val1
                        )
                    rect.append(newPoint)
                crossSection.append(rect)

    return crossSection

# Gets the centroid from rectangles defining the cross section
def getYAvg(rectangles):
    # Gets the total area of the cross section
    totalArea = sum([(rect[1][0]-rect[0][0])*(rect[0][1]-rect[1][1]) for rect in rectangles])
    weightedSum = 0

    # Calculates the weighed sum
    for rect in rectangles:
        weightedSum += (rect[1][0]-rect[0][0])*(rect[0][1]-rect[1][1])*(rect[1][1]+rect[0][1])/2

    return weightedSum/totalArea

# Gets the 2nd moment of area from the cross section based on the rectangles
def get2ndMomentOfArea(rectangles):
    yAvg = getYAvg(rectangles)

    result = 0

    for rect in rectangles:
        result += ((rect[1][0]-rect[0][0])*(rect[0][1]-rect[1][1])**3)/12 # b*h^3/12
        result += (rect[1][0]-rect[0][0])*(rect[0][1]-rect[1][1])*((yAvg-(rect[0][1]+rect[1][1])/2)**2) # A * (yAvg-y)^2
    
    return result

# Gets the positions and widths of the glue tabs at a position
def getGlueTab(position):
    glueTabs = [] # [[y-value, width] ...]

    # Interpolates between the key points with a line
    for tabKeyPoints in glueKeyPoints:
        for j in range(1, len(tabKeyPoints)):
            if position == tabKeyPoints[j][2]:
                glueTabs.append([tabKeyPoints[j][0], tabKeyPoints[j][1]])
                break

            if tabKeyPoints[j][2] > position:
                glueTabs.append([(tabKeyPoints[j][0]-tabKeyPoints[j-1][0])/(tabKeyPoints[j][2]-tabKeyPoints[j-1][2])*(position-tabKeyPoints[j-1][2])+tabKeyPoints[j-1][0],
                            (tabKeyPoints[j][1]-tabKeyPoints[j-1][1])/(tabKeyPoints[j][2]-tabKeyPoints[j-1][2])*(position-tabKeyPoints[j-1][2])+tabKeyPoints[j-1][1]
                            ])
                break
    
    return glueTabs

# Interpolates between the key points
def interpolateType12Failure(keyPoints, position):
    results = [] # [[y-value, thickness, width] ...]
    
    # Interpolates between the key points with a line
    for rectKeyPoints in keyPoints:
        for j in range(1, len(rectKeyPoints)):
            if position == rectKeyPoints[j][3]:
                results.append([rectKeyPoints[j][0], rectKeyPoints[j][1], rectKeyPoints[j][2]])
                break

            if rectKeyPoints[j][3] > position:
                results.append([(rectKeyPoints[j][0]-rectKeyPoints[j-1][0])/(rectKeyPoints[j][3]-rectKeyPoints[j-1][3])*(position-rectKeyPoints[j-1][3])+rectKeyPoints[j-1][0],
                            (rectKeyPoints[j][1]-rectKeyPoints[j-1][1])/(rectKeyPoints[j][3]-rectKeyPoints[j-1][3])*(position-rectKeyPoints[j-1][3])+rectKeyPoints[j-1][1],
                            (rectKeyPoints[j][2]-rectKeyPoints[j-1][2])/(rectKeyPoints[j][3]-rectKeyPoints[j-1][3])*(position-rectKeyPoints[j-1][3])+rectKeyPoints[j-1][2]
                            ])
                break
    
    return results

# Interpolates between the key points
def interpolateType3Failure(keyPoints, position):
    results = [] # [[y-value, thickness, width] ...]
    
    # Interpolates between the key points with a line
    for rectKeyPoints in keyPoints:
        for j in range(1, len(rectKeyPoints)):
            if position == rectKeyPoints[j][2]:
                results.append([rectKeyPoints[j][0], rectKeyPoints[j][1]])
                break

            if rectKeyPoints[j][2] > position:
                results.append([(rectKeyPoints[j][0]-rectKeyPoints[j-1][0])/(rectKeyPoints[j][2]-rectKeyPoints[j-1][2])*(position-rectKeyPoints[j-1][2])+rectKeyPoints[j-1][0],
                            (rectKeyPoints[j][1]-rectKeyPoints[j-1][1])/(rectKeyPoints[j][2]-rectKeyPoints[j-1][2])*(position-rectKeyPoints[j-1][2])+rectKeyPoints[j-1][1],
                            ])
                break
    
    return results

# Plots the cross section at a position
def plotCrossSectionAtPosition(position):
    rectangles = getCrossSection(position)
    yAvg = getYAvg(rectangles)
    I = get2ndMomentOfArea(rectangles)

    glueTabs = getGlueTab(position)
    glueTabs = ["Y-Value: " + str(i[0]) + "  Width: " + str(i[1]) for i in glueTabs]

    fig, ax1 = plt.subplots()
    for rect in rectangles:
        rectPlot = patches.Rectangle((rect[0][0], rect[1][1]), rect[1][0]-rect[0][0] , rect[0][1]-rect[1][1], facecolor = "gray", linewidth=1.5, edgecolor = "black")
        ax1.add_patch(rectPlot)
    
    # Plots the centroid
    ax1.axhline(y=yAvg, color='orange', linestyle='--', linewidth=1)

    fig.text(0.1, .05, "Glue Tabs: "+" | ".join(glueTabs), ha='center')
    fig.text(0.9, .05, "I: "+str(I), ha='center')
    ax1.set_xlabel("X (mm)")
    ax1.set_ylabel("Y (mm)", rotation=90)
    ax1.set_title("Cross Section at Position: "+str(position))
    ax1.set_xlim(min(rect[0][0] for rect in rectangles)-10, max(rect[1][0] for rect in rectangles)+10)
    ax1.set_ylim(min(rect[1][1] for rect in rectangles)-10, max(rect[0][1] for rect in rectangles)+10)
    ax1.set_aspect('equal', adjustable='box')
    plt.show()

# Plots the entire bridge
def plotBridge():
    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection="3d")
    intervalSize = 30
    for position in range(0, bridgeLength+1, intervalSize):
        rectangles = getCrossSection(position)
        for rect in rectangles:
            y = [rect[0][0], rect[1][0], rect[1][0], rect[0][0]]*2
            z = [rect[1][1], rect[1][1], rect[0][1], rect[0][1]]*2
            x = [position]*4+[position+intervalSize]*4
            vertices = np.array([
                [x[i], y[i], z[i]] for i in range(8)
            ])

            faces = [
                [vertices[0], vertices[1], vertices[2], vertices[3]],  # Bottom face
                [vertices[4], vertices[5], vertices[6], vertices[7]],  # Top face
                [vertices[0], vertices[1], vertices[5], vertices[4]],  # Front face
                [vertices[2], vertices[3], vertices[7], vertices[6]],  # Back face
                [vertices[1], vertices[2], vertices[6], vertices[5]],  # Right face
                [vertices[0], vertices[3], vertices[7], vertices[4]]   # Left face
            ]
            cube = Poly3DCollection(faces, facecolors='gray', alpha=1, linewidths=0.1, edgecolors='k')
            ax1.add_collection3d(cube)
    
    ax1.set_xlim([-10, 1710])
    ax1.set_ylim([-200, 100])
    ax1.set_zlim([-10, 210])
    plt.axis('off')
    plt.show()


############################################################################################################

# FORCES AND BENDING MOMENT ANALYSIS

# Gets point loads and position of each one based on train weight:
def calculatePointLoads(trainWeight, trainPosition):
    # trainPosition is the position of the front wheel of the locomotive
    # trainWeight is the total mass of the train
    # Distribution of mass is 1/(1+1.1+1.38*1.1) for the fist car
    # Distribution of mass is 1.1/(1+1.1+1.38*1.1) for the fist car
    # Distribution of mass is 1.38*1.1/(1+1.1+1.38*1.1) for the fist car

    car1 = -1/(1+1.1+1.38*1.1)*trainWeight
    car2 = -1.1/(1+1.1+1.38*1.1)*trainWeight
    car3 = -1.38*1.1/(1+1.1+1.38*1.1)*trainWeight

    if case1:
        car1 = -400/3
        car2 = -400/3
        car3 = -400/3

    # List in format: [[load, position], ...]
    out = [
        [car1/2, trainPosition-856],
        [car1/2, trainPosition-680],
        [car2/2, trainPosition-516],
        [car2/2, trainPosition-340],
        [car3/2, trainPosition-176],
        [car3/2, trainPosition]
    ]

    # Gets rid or all point loads not on the bridge
    ind = 0
    while True:
        if ind > len(out)-1:
            break

        if out[ind][1] < 0 or out[ind][1] > bridgeLength:
            out.pop(ind)
        else:
            ind += 1
    return out

# Cacluates reaction forces
def calculateReactions(pointLoads):
    reactionRight = -sum(i[0]*i[1] for i in pointLoads)/bridgeLength
    reactionLeft = -sum(i[0] for i in pointLoads)-reactionRight

    return [reactionLeft, reactionRight]


# Gets key points on the SFD
def getKeyPointsSFD(trainWeight, trainPosition):
    # Calcualates forces
    pointLoads = calculatePointLoads(trainWeight, trainPosition)
    reactions = calculateReactions(pointLoads)

    # Combines forces into one list
    allForces = [[reactions[0], 0]] + pointLoads + [[reactions[1], bridgeLength]]
    allForces.sort(key = lambda x: x[1])

    # Makes a list of key points on the SFD
    keyPoints = [[0, 0]]
    for i in sorted(allForces, key = lambda x: x[1]):
        keyPoints.append([keyPoints[-1][0]+i[0], i[1]])
    keyPoints.pop(0)

    return keyPoints

# Calculates SFD
def calculateSFD(trainWeight, trainPosition):
    # trainPosition is the position of the front wheel of the locomotive
    # trainWeight is the total mass of the train

    SFD = [] # List of form [shear force ...]

    # Calcualates forces
    pointLoads = calculatePointLoads(trainWeight, trainPosition)
    reactions = calculateReactions(pointLoads)

    # Combines forces into one list
    allForces = [[reactions[0], 0]] + pointLoads + [[reactions[1], bridgeLength]]
    allForces.sort(key = lambda x: x[1])

    # Cumulative sums of the forces for shear
    curIndex = 0
    curShear = 0
    for i in range(0, bridgeLength+1):
        maxCurShear = curShear
        while i > allForces[curIndex][1]:
            curShear += allForces[curIndex][0]
            if abs(curShear) > abs(maxCurShear):
                maxCurShear = curShear
            curIndex += 1
        SFD.append(maxCurShear)

    return SFD

# Gets the max SFD out of all possible train positions
def calculateMaxSFD(trainWeight):
    matrix = []

    # Goes through all possible position
    for i in range(0, bridgeLength+867):
        matrix.append(calculateSFD(trainWeight, i))
    
    maxSFDTension = []
    maxSFDCompression = []

    for i in range(0, len(matrix[0])):
        maxSFDTension.append([max([j[i] for j in matrix]), i])
        maxSFDCompression.append([min([j[i] for j in matrix]), i])

    maxSFD = [max(abs(maxSFDTension[i][0]), abs(maxSFDCompression[i][0])) for i in range(len(maxSFDTension))]
    return [maxSFD[i] if maxSFD[i] > maxSFD[-i-1] else maxSFD[-i-1] for i in range(len(maxSFD))]

# Gets the lower bound on the SFD
def calculateMinSFD(trainWeight):
    return [-i for i in calculateMaxSFD(trainWeight)]

# Gets key points on the BMD
def getKeyPointsBMD(trainWeight, trainPosition):
    # Calcualates forces
    pointLoads = calculatePointLoads(trainWeight, trainPosition)
    reactions = calculateReactions(pointLoads)

    # Combines forces into one list
    allForces = [[reactions[0], 0]] + pointLoads + [[reactions[1], bridgeLength]]
    allForces.sort(key = lambda x: x[1])

    # Makes a list of key points on the SFD
    keyPoints = [[0, 0]]
    for i in sorted(allForces, key = lambda x: x[1]):
        keyPoints.append([keyPoints[-1][0]+i[0], i[1]])
    keyPoints.pop(0)

    # Gets the key points on the BMD
    keyMomentPoints = [[0, 0]]
    for i in range(1,len(keyPoints)):
        keyMomentPoints.append([keyMomentPoints[-1][0] + (keyPoints[i-1][0])*(keyPoints[i][1]-keyPoints[i-1][1]),
                                keyPoints[i][1]])
    
    return keyMomentPoints

# Calculates BMD
def calculateBMD(trainWeight, trainPosition):
    # trainPosition is the position of the front wheel of the locomotive
    # trainWeight is the total mass of the train

    BMD = [] # List of form [moment ...]

    # Calcualates forces
    pointLoads = calculatePointLoads(trainWeight, trainPosition)
    reactions = calculateReactions(pointLoads)

    # Combines forces into one list
    allForces = [[reactions[0], 0]] + pointLoads + [[reactions[1], bridgeLength]]
    allForces.sort(key = lambda x: x[1])

    # Makes a list of key points on the SFD
    keyPoints = [[0, 0]]
    for i in sorted(allForces, key = lambda x: x[1]):
        keyPoints.append([keyPoints[-1][0]+i[0], i[1]])
    keyPoints.pop(0)

    keyMomentPoints = [[0, 0]]
    for i in range(1,len(keyPoints)):
        keyMomentPoints.append([keyMomentPoints[-1][0] + (keyPoints[i-1][0])*(keyPoints[i][1]-keyPoints[i-1][1]),
                                keyPoints[i][1]])

    for i in range(0, bridgeLength+1):        
        # Interpolates between the key points with a line
        for j in range(1, len(keyMomentPoints)):
            if i == keyMomentPoints[j][1]:
                BMD.append(keyMomentPoints[j][0])
                break

            if keyMomentPoints[j][1] > i:
                BMD.append((keyMomentPoints[j][0]-keyMomentPoints[j-1][0])/(keyMomentPoints[j][1]-keyMomentPoints[j-1][1])*(i-keyMomentPoints[j-1][1])+keyMomentPoints[j-1][0])
                break

    return BMD

# Gets the max BMD out of all possible train positions
def calculateMaxBMD(trainWeight):
    matrix = []

    # Goes through all possible position
    for i in range(0, bridgeLength+867):
        matrix.append(calculateBMD(trainWeight, i))
    
    maxBMDPositive = []
    maxBMDNegative = []

    for i in range(0, len(matrix[0])):
        maxBMDPositive.append([max([j[i] for j in matrix]), i])
        maxBMDNegative.append([min([j[i] for j in matrix]), i])
    
    maxBMD = [max(abs(maxBMDPositive[i][0]), abs(maxBMDNegative[i][0])) for i in range(len(maxBMDNegative))]
    return [maxBMD[i] if maxBMD[i] > maxBMD[-i-1] else maxBMD[-i-1] for i in range(len(maxBMD))]


# Calculates the curvature diagram from a BMD
def calculateCurvatureDiagram(BMD):
    curvatureDiagram = [] # [[curvature, position] ...]

    for position in range(len(BMD)):
        # Gets the 2nd Moment of Area
        rectangles = getCrossSection(position)
        SecondMomentOfArea = get2ndMomentOfArea(rectangles)

        # M/(EI)
        curvatureDiagram.append([BMD[position]/(matboardYoungsModulus * SecondMomentOfArea), position])
    
    return curvatureDiagram


############################################################################################################

# TENSION/COMPRESSION FAILURE

# Calculates moment required for tensile failure
def calculateTensileFailureMomentDiagram():
    tensileFailureMoment = [] # [[moment, position]]

    # Goes through all positions
    for position in range(bridgeLength + 1):
        # Gets geometrical properties
        crossSection = getCrossSection(position)
        yAvg = getYAvg(crossSection)
        secondMomentOfArea = get2ndMomentOfArea(crossSection)

        # Calculates the maximum y dist for tension
        maxYTension = abs(yAvg-min([rect[1][1] for rect in crossSection]))

        # Gets the moment corresponding to that distance
        maxM = matboardTensileStrength*secondMomentOfArea/maxYTension

        tensileFailureMoment.append([maxM, position])
    
    return tensileFailureMoment

# Calculates moment required for compressive failure
def calculateCompressiveFailureMomentDiagram():
    compressiveFailureMoment = [] # [[moment, position]]

    # Goes through all positions
    for position in range(bridgeLength + 1):
        # Gets geomtrical properties
        crossSection = getCrossSection(position)
        yAvg = getYAvg(crossSection)
        secondMomentOfArea = get2ndMomentOfArea(crossSection)

        # Calculates the maximum y dist for compression
        maxYCompression = abs(yAvg-max([rect[0][1] for rect in crossSection]))

        # Gets the moment corresponding to that distance
        maxM = matboardCompressiveStrength*secondMomentOfArea/maxYCompression
        
        compressiveFailureMoment.append([maxM, position])
    
    return compressiveFailureMoment


############################################################################################################

# SHEAR FAILURE

# Calculates for shear failure for matboard
def calculateMatboardFailureShearDiagram():
    shearFailure = [] # [[shear, position, height]]

    # Goes through all possible positions
    for position in range(bridgeLength + 1):
        # Gets the geometrical properties
        crossSection = getCrossSection(position)
        yAvg = getYAvg(crossSection)
        secondMomentOfArea = get2ndMomentOfArea(crossSection)

        # Looks at all q values to find the maximum shear + centroid
        maxShear = float("inf")
        height = yAvg
        for yVal in list(np.linspace(min([rect[1][1] for rect in crossSection]), max([rect[0][1] for rect in crossSection]), 500)) + [yAvg]:

            # Gets the width at the specified hieght
            width = 0
            for rect in crossSection:
                if rect[1][1] <= yVal and rect[0][1] >= yVal:
                    width += rect[1][0]-rect[0][0]
            
            # Gets the Q value at the heights
            Q = 0
            for rect in crossSection:
                if rect[1][1] >= yVal and rect[0][1] >= yVal:
                    Q += abs((rect[1][0]-rect[0][0])*(rect[1][1]-rect[0][1]))*((rect[1][1]+rect[0][1])/2-yAvg)
                elif rect[1][1] <= yVal and rect[0][1] >= yVal:
                    Q += abs((rect[1][0]-rect[0][0])*(rect[0][1]-yVal))*((rect[0][1]+yVal)/2-yAvg)
            if Q != 0:
                shear = matboardShearStrength*secondMomentOfArea*width/abs(Q) # tau*I*b/Q = V
                if shear <= maxShear:
                    maxShear = shear
                    height = yVal
        
        shearFailure.append([maxShear, position, height])
    
    return shearFailure

# Calculates for shear failure at glue joints
def calculateGlueFailureShearDiagram():
    shearFailure = [] # [[shear, position, glue tab index]]

    # Goes through all possible positions
    for position in range(bridgeLength + 1):
        # Gets the geometrical properties
        crossSection = getCrossSection(position)
        yAvg = getYAvg(crossSection)
        secondMomentOfArea = get2ndMomentOfArea(crossSection)
        glueTabs = getGlueTab(position)

        # Looks at all glue distances to find the maximum shear
        maxShear = float("inf")
        ind = 0
        truInd = 0
        for yVal, width in glueTabs:
            # Calculates Q at the height
            Q = 0
            for rect in crossSection:
                if rect[1][1] >= yVal and rect[0][1] >= yVal:
                    Q += abs((rect[1][0]-rect[0][0])*(rect[1][1]-rect[0][1]))*((rect[1][1]+rect[0][1])/2-yAvg)
                elif rect[1][1] <= yVal and rect[0][1] >= yVal:
                    Q += abs((rect[1][0]-rect[0][0])*(rect[0][1]-yVal))*((rect[0][1]+yVal)/2-yAvg)
            if Q != 0:
                shear = cementShearStrength*secondMomentOfArea*width/abs(Q) # tau*I*b/Q = V
                if shear <= maxShear:
                    maxShear = shear
                    truInd = ind
            ind += 1
        
        shearFailure.append([maxShear, position, truInd])
    
    return shearFailure

# Calculates for shear failure at a specific glue joint
def calculateGlueFailureShearDiagramSpecific(glueIndex):
    shearFailure = [] # [[shear, position]]

    # Goes through all possible positions
    for position in range(bridgeLength + 1):
        # Gets the geometrical properties
        crossSection = getCrossSection(position)
        yAvg = getYAvg(crossSection)
        secondMomentOfArea = get2ndMomentOfArea(crossSection)
        glueTabs = getGlueTab(position)[glueIndex]

        # Looks at all glue distances to find the maximum shear
        maxShear = float("inf")
        for yVal, width in [glueTabs]:
            # Calculates Q at the height
            Q = 0
            for rect in crossSection:
                if rect[1][1] >= yVal and rect[0][1] >= yVal:
                    Q += abs((rect[1][0]-rect[0][0])*(rect[1][1]-rect[0][1]))*((rect[1][1]+rect[0][1])/2-yAvg)
                elif rect[1][1] <= yVal and rect[0][1] >= yVal:
                    Q += abs((rect[1][0]-rect[0][0])*(rect[0][1]-yVal))*((rect[0][1]+yVal)/2-yAvg)
            if Q != 0:
                shear = cementShearStrength*secondMomentOfArea*width/abs(Q) # tau*I*b/Q = V
                if shear <= maxShear:
                    maxShear = shear
        
        shearFailure.append([maxShear, position])
    
    return shearFailure

############################################################################################################

# BUCKLING FAILURE

# Type 1: 4 Fixed Ends, Constant Stress
def calculateType1FailureMomentDiagramSpecific(rectIndex):
    k = 4

    momentFailure = []
    # Goes through all positions
    for position in range(bridgeLength + 1):
        # Gets the geometry at that position
        y, t, b = interpolateType12Failure(type1FailureKeyPoints, position)[rectIndex]
        crossSection = getCrossSection(position)
        I = get2ndMomentOfArea(crossSection)
        yAvg = getYAvg(crossSection)
        y = y-yAvg

        # Calculates the moment required for failure
        if y <= 0:
            momentFailure.append([float("inf"), position])
            continue
        moment = ((k*pi**2*matboardYoungsModulus)/(12*(1-matboardPoissons**2))*(t/b)**2) *I/y
        momentFailure.append([moment, position])
    
    return momentFailure

def calculateType1FailureMomentDiagram():
    matrix = []

    for i in range(len(type1FailureKeyPoints)):
        matrix.append(calculateType1FailureMomentDiagramSpecific(i))
    
    maxBMD = []

    for i in range(0, len(matrix[0])):
        maxBMD.append([min([j[i][0] for j in matrix]), i])

    return maxBMD

# Type 2: 1 Open End, Constant Stress
def calculateType2FailureMomentDiagramSpecific(rectIndex):
    k = 0.425

    momentFailure = []
    # Goes through all positions
    for position in range(bridgeLength + 1):
        # Gets the geometry at that position
        y, t, b = interpolateType12Failure(type2FailureKeyPoints, position)[rectIndex]
        crossSection = getCrossSection(position)
        I = get2ndMomentOfArea(crossSection)
        yAvg = getYAvg(crossSection)
        y = y-yAvg

        # Calculates the moment required for failure
        if y <= 0:
            momentFailure.append([float("inf"), position])
            continue
        moment = ((k*pi**2*matboardYoungsModulus)/(12*(1-matboardPoissons**2))*(t/b)**2) *I/y
        momentFailure.append([moment, position])
    return momentFailure

def calculateType2FailureMomentDiagram():
    matrix = []

    for i in range(len(type2FailureKeyPoints)):
        matrix.append(calculateType2FailureMomentDiagramSpecific(i))
    
    maxBMD = []

    for i in range(0, len(matrix[0])):
        maxBMD.append([min([j[i][0] for j in matrix]), i])

    return maxBMD

# Type 3: Varying Stress
def calculateType3FailureMomentDiagramSpecific(rectIndex):
    k = 6

    momentFailure = []
    # Goes through all positions
    for position in range(bridgeLength + 1):
        # Gets the geometry at that position
        y, t = interpolateType3Failure(type3FailureKeyPoints, position)[rectIndex]
        crossSection = getCrossSection(position)
        I = get2ndMomentOfArea(crossSection)
        yAvg = getYAvg(crossSection)
        y = y-yAvg

        # Calculate the moment required for failure
        if y <= 0:
            momentFailure.append([float("inf"), position])
            continue
        moment = ((k*pi**2*matboardYoungsModulus)/(12*(1-matboardPoissons**2))*(t/y)**2) *I/y
        momentFailure.append([moment, position])
    return momentFailure

def calculateType3FailureMomentDiagram():
    matrix = []

    for i in range(len(type3FailureKeyPoints)):
        matrix.append(calculateType3FailureMomentDiagramSpecific(i))
    
    maxBMD = []

    for i in range(0, len(matrix[0])):
        maxBMD.append([min([j[i][0] for j in matrix]), i])

    return maxBMD

# Type 4: Shear Stress
@lru_cache(maxsize = 2*bridgeLength)
def getSmallestThickness(lowerBound, upperBound, rectIndex):
    t = float("inf")
    # Finds the smallest thickness out of the key points
    for rect in type4FailureKeyPoints[rectIndex]:
        if rect[2] >= lowerBound and rect[2] <= upperBound and rect[1] < t:
            t = rect[1]

    # Iterates through the key points
    for i in range(len(type4FailureKeyPoints[rectIndex])-1):
        # Gets the point after and the point before
        rect1 = type4FailureKeyPoints[rectIndex][i]
        rect2 = type4FailureKeyPoints[rectIndex][i+1]

        # Gets the smallest thickness if it is not at one of the key points
        if rect1[2] <= lowerBound and rect2[2] >= lowerBound:
            newT = rect1[1]+(rect2[1]-rect1[1])/(rect2[2]-rect1[2])*(lowerBound-rect1[2])
            if newT < t:
                t = newT
        if rect1[2] <= upperBound and rect2[2] >= upperBound:
            newT = rect1[1]+(rect2[1]-rect1[1])/(rect2[2]-rect1[2])*(upperBound-rect1[2])
            if newT < t:
                t = newT
    return t

@lru_cache(maxsize = 2*bridgeLength)
def getLargestHeight(lowerBound, upperBound, rectIndex):
    h = 0
    # Finds the smallest height out of the key points
    for rect in type4FailureKeyPoints[rectIndex]:
        if rect[2] >= lowerBound and rect[2] <= upperBound and rect[0] > h:
            h = rect[0]

    # Iterates through the key points
    for i in range(len(type4FailureKeyPoints[rectIndex])-1):
        # Gets the point after and the point before
        rect1 = type4FailureKeyPoints[rectIndex][i]
        rect2 = type4FailureKeyPoints[rectIndex][i+1]

        # Gets the smallest height if it is not at one of the key points
        if rect1[2] <= lowerBound and rect2[2] >= lowerBound:
            newH = rect1[0]+(rect2[0]-rect1[0])/(rect2[2]-rect1[2])*(lowerBound-rect1[2])
            if newH > h:
                h = newH
        if rect1[2] <= upperBound and rect2[2] >= upperBound:
            newH = rect1[0]+(rect2[0]-rect1[0])/(rect2[2]-rect1[2])*(upperBound-rect1[2])
            if newH > h:
                h = newH
    return h

@lru_cache(maxsize = 2*bridgeLength)
def getSmallestITimesWidthOverQ(position2):
    # Gets geomtrical properties
    crossSection = getCrossSection(position2)
    I = get2ndMomentOfArea(crossSection)
    yAvg = getYAvg(crossSection)

    widthOverQ = float("inf")

    for yVal in list(np.linspace(min([rect[1][1] for rect in crossSection]), max([rect[0][1] for rect in crossSection]), 500)) + [yAvg]:
        # Gets the width at the specified hieght
        width = 0
        for rect in crossSection:
            if rect[1][1] <= yVal and rect[0][1] >= yVal:
                width += rect[1][0]-rect[0][0]
        
        # Gets the Q value at the heights
        Q = 0
        for rect in crossSection:
            if rect[1][1] >= yVal and rect[0][1] >= yVal:
                Q += abs((rect[1][0]-rect[0][0])*(rect[1][1]-rect[0][1]))*((rect[1][1]+rect[0][1])/2-yAvg)
            elif rect[1][1] <= yVal and rect[0][1] >= yVal:
                Q += abs((rect[1][0]-rect[0][0])*(rect[0][1]-yVal))*((rect[0][1]+yVal)/2-yAvg)
        if Q != 0:
            # Checks if it is the smallest
            if I*width/abs(Q) < widthOverQ:
                widthOverQ = I*width/abs(Q)
    return widthOverQ

def calculateType4FailureShearDiagramSpecific(rectIndex):
    k = 5

    shearFailure = []

    # Iterates through all positions
    for position in range(bridgeLength + 1):
        # Finds the closes diaphragms
        lowerBound = [x for x in diaphragmPositions if x < position]
        upperBound = [x for x in diaphragmPositions if x >= position]
        if len(lowerBound) == 0:
            lowerBound = 0
            upperBound = diaphragmPositions[1]
        elif len(upperBound) == 0:
            lowerBound = diaphragmPositions[-2]
            upperBound = diaphragmPositions[-1]
        else:
            lowerBound = lowerBound[-1]
            upperBound = upperBound[0]
        
        a = upperBound-lowerBound

        # Gets the smallest thickness between the two diaphragms
        t = getSmallestThickness(lowerBound, upperBound, rectIndex)
        
        # Gets the largest height between the two diaphragms
        h = getLargestHeight(lowerBound, upperBound, rectIndex)

        # Gets smallest shear between the two diaphragms
        smallestShear = float("inf")
        for position2 in list(np.linspace(lowerBound, upperBound, 1000))+[position]:
            shear = (k*pi**2*matboardYoungsModulus)/(12*(1-matboardPoissons**2))*((t/a)**2+(t/h)**2)*getSmallestITimesWidthOverQ(position2)
            if shear < smallestShear:
                smallestShear = shear

        shearFailure.append([smallestShear, position])
    
    return shearFailure

def calculateType4FailureShearDiagram():
    matrix = []
    for i in range(len(type4FailureKeyPoints)):
        matrix.append(calculateType4FailureShearDiagramSpecific(i))
    
    maxSFD = []

    for i in range(0, len(matrix[0])):
        maxSFD.append([min([j[i][0] for j in matrix]), i])

    return maxSFD


############################################################################################################

# PLOTTING THE RESULTS

# Makes the SFD comparing the max allowable for all failure modes to the actual
@lru_cache(maxsize = 2*bridgeLength)
def plotSFD(trainWeight):
    maxSFD = calculateMaxSFD(trainWeight)
    minSFD = calculateMinSFD(trainWeight)

    # Max shear forces
    matboardFailureShear = [i[0] for i in calculateMatboardFailureShearDiagram()]
    glueFailureShear = [i[0] for i in calculateGlueFailureShearDiagram()]
    type4FailureShear = [i[0] for i in calculateType4FailureShearDiagram()]

    fig, ax1 = plt.subplots()

    # Plots the envelope
    ax1.plot([0]+list(range(len(maxSFD)))+[len(maxSFD)-1], [0]+maxSFD+[0], label='Positive Shear Envelope', color='b')
    ax1.plot([0]+list(range(len(minSFD)))+[len(minSFD)-1], [0]+minSFD+[0], label='Negative Shear Envelope', color='g')

    # Plots the maximum shears for each failure
    ax1.plot(list(range(len(matboardFailureShear))), matboardFailureShear, label='Matboard Failure Max Shear', color='r')
    ax1.plot(list(range(len(glueFailureShear))), glueFailureShear, label='Glue Failure Max Shear', color='m')
    ax1.plot(list(range(len(type4FailureShear))), type4FailureShear, label='Type 4 Buckling Max Shear', color='k')

    ax1.plot(list(range(len(matboardFailureShear))), [-i for i in matboardFailureShear], color='r')
    ax1.plot(list(range(len(glueFailureShear))), [-i for i in glueFailureShear], color='m')
    ax1.plot(list(range(len(type4FailureShear))), [-i for i in type4FailureShear], color='k')

    ax1.set_xlabel("Position (mm)")
    ax1.set_ylabel("Shear Force (N)", rotation=90)
    ax1.set_title("Max Shear Force Diagram (SFD)")
    
    fig.legend()

    plt.show()

# Makes the BMD comparing the max allowable for all failure modes to the actual
@lru_cache(maxsize = 2*bridgeLength)
def plotBMD(trainWeight):
    BMD = calculateMaxBMD(trainWeight)
    
    # Max moments
    tensileFailureMoment = [i[0] for i in calculateTensileFailureMomentDiagram()]
    compressiveFailureMoment = [i[0] for i in calculateCompressiveFailureMomentDiagram()]
    type1FailureMoment = [i[0] for i in calculateType1FailureMomentDiagram()]
    type2FailureMoment = [i[0] for i in calculateType2FailureMomentDiagram()]
    type3FailureMoment = [i[0] for i in calculateType3FailureMomentDiagram()]

    fig, ax1 = plt.subplots()

    # Plots the envelope
    ax1.plot(list(range(len(BMD))), BMD, label='Bending Moment Envelope', color='b')

    # Plots the maximum shears for each failure
    ax1.plot(list(range(len(tensileFailureMoment))), tensileFailureMoment, label='Tensile Failure Max Moment', color='r')
    ax1.plot(list(range(len(compressiveFailureMoment))), compressiveFailureMoment, label='Compressive Failure Max Moment', color='c')
    ax1.plot(list(range(len(type1FailureMoment))), type1FailureMoment, label='Type 1 Buckling Max Moment', color='m')
    ax1.plot(list(range(len(type2FailureMoment))), type2FailureMoment, label='Type 2 Buckling Max Moment', color='y')
    ax1.plot(list(range(len(type3FailureMoment))), type3FailureMoment, label='Type 3 Buckling Max Moment', color='k')

    ax1.set_xlabel("Position (mm)")
    ax1.set_ylabel("Bending Moment (Nmm)", rotation=90)
    ax1.set_title("Max Bending Moment Diagram (BMD)")

    fig.legend()

    plt.show()


############################################################################################################

# FINAL CONCLUSIONS

# Gets the max load the bridge can handle
def calculateMaxLoad():
    SFD = calculateMaxSFD(1)
    BMD = calculateMaxBMD(1)

    # Max shear forces
    matboardFailureShear = [i[0] for i in calculateMatboardFailureShearDiagram()]
    glueFailureShear = [i[0] for i in calculateGlueFailureShearDiagram()]
    type4FailureShear = [i[0] for i in calculateType4FailureShearDiagram()]

    # Max moments
    tensileFailureMoment = [i[0] for i in calculateTensileFailureMomentDiagram()]
    compressiveFailureMoment = [i[0] for i in calculateCompressiveFailureMomentDiagram()]
    type1FailureMoment = [i[0] for i in calculateType1FailureMomentDiagram()]
    type2FailureMoment = [i[0] for i in calculateType2FailureMomentDiagram()]
    type3FailureMoment = [i[0] for i in calculateType3FailureMomentDiagram()]

    # Gets FOS for all points on the bridge
    mults = []
    for i in range(len(SFD)):
        mults.append(min([
            matboardFailureShear[i]/SFD[i] if SFD[i] > 0 else float("inf"),
            glueFailureShear[i]/SFD[i] if SFD[i] > 0 else float("inf"),
            type4FailureShear[i]/SFD[i] if SFD[i] > 0 else float("inf"),
            tensileFailureMoment[i]/BMD[i] if BMD[i] > 0 else float("inf"),
            compressiveFailureMoment[i]/BMD[i] if BMD[i] > 0 else float("inf"),
            type1FailureMoment[i]/BMD[i] if BMD[i] > 0 else float("inf"),
            type2FailureMoment[i]/BMD[i] if BMD[i] > 0 else float("inf"),
            type3FailureMoment[i]/BMD[i] if BMD[i] > 0 else float("inf")
        ]))
    
    # Smallest FOS
    return min(mults)

# Gets the safety factors for each of the failure methods
def calculateSafetyFactors(load = -1):
    if load == -1:
        maxLoad = calculateMaxLoad()
    else:
        maxLoad = load
    SFD = calculateMaxSFD(maxLoad)
    BMD = calculateMaxBMD(maxLoad)

    # Max shear forces
    matboardFailureShear = [i[0] for i in calculateMatboardFailureShearDiagram()]
    matboardFailureShear = min(matboardFailureShear[i]/SFD[i] for i in range(len(SFD)) if SFD[i] > 0)
    glueFailureShear = [i[0] for i in calculateGlueFailureShearDiagram()]
    glueFailureShear = min(glueFailureShear[i]/SFD[i] for i in range(len(SFD)) if SFD[i] > 0)
    type4FailureShear = [i[0] for i in calculateType4FailureShearDiagram()]
    type4FailureShear = min(type4FailureShear[i]/SFD[i] for i in range(len(SFD)) if SFD[i] > 0)

    # Max moments
    tensileFailureMoment = [i[0] for i in calculateTensileFailureMomentDiagram()]
    tensileFailureMoment = min(tensileFailureMoment[i]/BMD[i] for i in range(len(BMD)) if BMD[i] > 0)
    compressiveFailureMoment = [i[0] for i in calculateCompressiveFailureMomentDiagram()]
    compressiveFailureMoment = min(compressiveFailureMoment[i]/BMD[i] for i in range(len(BMD)) if BMD[i] > 0)
    type1FailureMoment = [i[0] for i in calculateType1FailureMomentDiagram()]
    type1FailureMoment = min(type1FailureMoment[i]/BMD[i] for i in range(len(BMD)) if BMD[i] > 0)
    type2FailureMoment = [i[0] for i in calculateType2FailureMomentDiagram()]
    type2FailureMoment = min(type2FailureMoment[i]/BMD[i] for i in range(len(BMD)) if BMD[i] > 0)
    type3FailureMoment = [i[0] for i in calculateType3FailureMomentDiagram()]
    type3FailureMoment = min(type3FailureMoment[i]/BMD[i] for i in range(len(BMD)) if BMD[i] > 0)

    return {
        "Matboard Shear Failure" : matboardFailureShear,
        "Glue Shear Failure" : glueFailureShear,
        "Matboard Tensile Failure" : tensileFailureMoment,
        "Matboard Compressive Failure" : compressiveFailureMoment,
        "Type 1 Failure" : type1FailureMoment,
        "Type 2 Failure" : type2FailureMoment,
        "Type 3 Failure" : type3FailureMoment,
        "Type 4 Failure" : type4FailureShear,
    }



plotBridge()
plotCrossSectionAtPosition(0)

# Gets the maximum load
m = calculateMaxLoad()
print("Max Possible Load: " + f"{m:.4g}")

# Gets the safety factors
sf = calculateSafetyFactors(m)
for i in sf:
    print("Failure Mode:", i, "   FOS: " + f"{sf[i]:.4g}")

# Plots the data
plotSFD(m)
plotBMD(m)

# Plots SFDs for the train moving over the bridge in intervals of 30mm
SFDS = []
for i in range(0, bridgeLength+856, 30):
    SFDS.append(calculateSFD(400, i))

maxSFD = calculateMaxSFD(400)
minSFD = calculateMinSFD(400)
fig, ax1 = plt.subplots()
ax1.set_xlabel("Position (mm)")
ax1.set_ylabel("Shear Force (N)", rotation=90)
ax1.set_title("Shear Force Diagrams (SFD) for Different Train Positions")

for i in SFDS:
    ax1.plot(range(len(i)), i)
ax1.plot(list(range(len(maxSFD)))+[len(maxSFD)-1], maxSFD+[0], lw=3, label="Positive Shear Envelope")
ax1.plot(list(range(len(minSFD)))+[len(minSFD)-1], minSFD+[0], lw=3, label="Negative Shear Envelope")
fig.legend()

    
# Plots BMDs for the train moving across the bridge in intervals of 30mm
BMD = calculateMaxBMD(400)
BMDS = []
for i in range(0, bridgeLength+856, 30):
    BMDS.append(calculateBMD(400, i))

fig2, ax2 = plt.subplots()
ax2.set_xlabel("Position (mm)")
ax2.set_ylabel("Bending Moment (Nmm)", rotation=90)
ax2.set_title("Bending Moment Diagrams (BMD) for Different Train Positions")

for i in BMDS:
    ax2.plot(range(len(i)), i)
ax2.plot(list(range(len(BMD))), BMD, lw=3, label="Bending Moment Envelope")
fig2.legend()

plt.show()
import os, sys, shutil

from OrganicGrowth import OrganicGrowth
from MaxP import MaxP
from Analysis import *

from util import *

#---Common os methods---
getDir = os.path.dirname
joinPath = os.path.join
checkDir = os.path.isdir
checkFile = os.path.isfile
makeDir = os.mkdir
makeDirs = os.makedirs
deleteDir = shutil.rmtree

#---Specify working directories (hard-coded)---
prjRootDir = getDir(getDir(getDir(sys.argv[0]))) # sys.argv[0] == \project\Codes\src\main.py
experimentTitle = "Hurricane_Sandy"
workDir = joinPath(prjRootDir, "Simulation", experimentTitle)
if not checkDir(workDir):
    makeDirs(workDir)
#---End


'''------Input candidates for this experiment-----'''
"""
Hurricane Sandy
QPERCAP & QFAM & QRICH200K & MDGRENT & MHSEVAL should change to negative
72 (NJ) + 64 (NY) coastal tracts = 136
p = 27, 2p = 54, 3p = 81, 4p = 108, 5p = 135
"""
# Source Feature dataset for region growth
maxRegion = [joinPath(prjRootDir, "Data","SandyData.gdb", "SandyArea"), 
                joinPath(prjRootDir, "Data","SandyData.gdb", "SandyArea_NJ"), 
                joinPath(prjRootDir, "Data","SandyData.gdb", "SandyArea_NY")]
# destined Feature dataset for grown region visual analysis 
dstFeaturesDir = joinPath(workDir, "FeaturesForAnalysis")
if not checkDir(dstFeaturesDir):
    makeDir(dstFeaturesDir) 
resultFeature = [joinPath(dstFeaturesDir, "SandyArea.shp"), 
                 joinPath(dstFeaturesDir, "SandyArea_NJCoastal.shp"), 
                 joinPath(dstFeaturesDir, "SandyArea_NYCoastal.shp"), 
                 joinPath(dstFeaturesDir, "SandyArea_NJState_NJCoastal.shp"), 
                 joinPath(dstFeaturesDir, "SandyArea_NYState_NYCoastal.shp"),
                 joinPath(dstFeaturesDir, "SandyArea_NJState_Newwark.shp"),
                 joinPath(dstFeaturesDir, "SandyArea_NYState_StateIsland.shp"),
                 joinPath(dstFeaturesDir, "SandyArea_NJState_Newward_OG_1.shp"),
                 joinPath(dstFeaturesDir, "SandyArea_NJState_Newward_OG_9.shp")]
# text input files directory
inDir = joinPath(workDir, "Input")
if(not checkDir(inDir)):
    makeDirs(inDir)
# exclude features due to no data
excludedFeatureDataPath_Sandy = joinPath(inDir, "excludeFeatures.txt")
# indicator data
indicatorFilePath_Sandy = joinPath(inDir, "Indicators.txt")
indicatorDirsFilePath_Sandy = joinPath(inDir, "IndicatorVunDir.txt")
# base region data
seedRegionData = [joinPath(inDir, "Coastal_all_tracts.txt"), 
                  joinPath(inDir, "Coastal_nj_tracts.txt"), 
                  joinPath(inDir, "Coastal_ny_tracts.txt"),
                  joinPath(inDir, "Newark_tracts.txt"),
                  joinPath(inDir, "StateIsland_tracts.txt"),
                  joinPath(inDir, "Newark_OG_1_tracts.txt"),
                  joinPath(inDir, "Newark_OG_9_tracts.txt")]
'''------End------'''

'''---Experiment setting---'''
# Set title for output folders' and files' names
maxRegionIndex = 1
maxGrownRegionTitle = ["GrowTo_EntireArea", "GrowTo_NJState", "GrowTo_NYState"][maxRegionIndex]
srcFeaturePath = maxRegion[maxRegionIndex]
# Set base region & feature result
seedRegionIndex = 3
baseRegionTitle = ["Coastal_Area_All", "Coastal_NJ", "Coastal_NY", "Newark", "StateIsland", "Newark_OG_1", "Newark_OG_9"][seedRegionIndex]
baseRegionFile = seedRegionData[seedRegionIndex]
dstFeaturePath = resultFeature[5]
# Set excluded data and indicators
excludedFeatureDataPath = excludedFeatureDataPath_Sandy
indicatorFilePath = indicatorFilePath_Sandy
indicatorDirsFilePath = indicatorDirsFilePath_Sandy
# Set Object/FID to a primary key field for result output
FID_keyField = "GEOID"
'''------End------'''

OGBdRN = None

'''---Organic Growth Setting---'''
# Set spatial contiguity relationship for organic growth
spatialContRel = ["BOUNDARY_TOUCHES", "SHARE_A_LINE_SEGMENT_WITH"][0]
'''------End------'''

'''---Max-p Setting---'''
constrainField = "ENUM"
# constrainField = "V5"
constrainList_SandyArea = sorted([3929, 1964, 982, 491, 245, 122, 61, 30, 15, 7])
# constrainList_SandyArea_NJState = sorted([1205, 602, 301, 150, 75, 37, 18, 9, 4])
nNJ, nNY = 1205, 2724
NJRegions, NYRegions = 32, 42
constrainList_SandyArea_NJState = np.sort((np.round(nNJ/np.arange(1.0, NJRegions+0.1, 1.0))-1)).astype(int).tolist()
constrainList_SandyArea_NYState = np.sort((np.round(nNY/np.arange(1.0, NYRegions+0.1, 1.0))-1)).astype(int).tolist()
constrainList = constrainList_SandyArea_NJState
'''------End------'''

#---Set output folders' and files' names---
# Figure & Chart output dir for all study areas
# figOutDir = joinPath(workDir, "All_Fig_charts")
# if not checkDir(figOutDir):
#     makeDir(figOutDir)
# figPrefix = maxGrownRegionTitle[7:]
# Set root directory for all outputs
outDir = joinPath(workDir, maxGrownRegionTitle)
if not checkDir(outDir):
    makeDir(outDir)
# Set sub directory for feature data output
featureInfoDir = joinPath(outDir, "FeatureDataInfo")
if not checkDir(featureInfoDir):
    makeDir(featureInfoDir)
# Set sub directory and final output file for organic growth analysis
organicGrowthFinalOutputFile = joinPath(outDir, "%s_OrganicGrowth_Result.txt" % baseRegionTitle)
organicGrowthPCADir = joinPath(outDir, "%s_OrganicGrowth_PCA" % baseRegionTitle)
# Set sub directory and final output file for max-p analysis
maxpRegionalizationDir = joinPath(outDir, "%s_MaxP_Regionalizations" % maxGrownRegionTitle)
maxpFinalOutputFile = joinPath(outDir, "%s_MaxP_Result.txt" % baseRegionTitle)
maxpPCADir = joinPath(outDir, "%s_MaxP_PCA" % baseRegionTitle)
#------End------

#------Feature dataset related operations------
#define object pointing to feature dataset
censusData = FeatureDataSet(srcFeaturePath, dstFeaturePath, outDir, FID_keyField, excludedFeatureDataPath)
#ask user if new spatial information needed to be derived for PCA
attributeTableFilePath, spatialContiguityFilePath = censusData.checkFeatureInfoDerivation(indicatorFilePath, spatialContRel, featureInfoDir)
#read in attribute values and spatial contiguity
numberOfFeature, fields, attribute2dList = readFieldValueFile(attributeTableFilePath)
spatialContiguityDict = readAdjacencyFile(spatialContiguityFilePath)
#------End------
 
# #---Run Organic Growth---
# organicGrowthExp = OrganicGrowth(censusData, spatialContiguityDict, baseRegionFile)
# organicGrowthExp.grow()
# organicGrowthExp.outputResultsToFeature()
# organicGrowthIdsInRegion = [list(set().union(*organicGrowthExp.layers[:i+1])) for i in xrange(len(organicGrowthExp.layers))]
# # Analysis organic growth results
# SVIDictList, RotatedCompMs, KMO_Bartletts = runPCA(organicGrowthIdsInRegion, attribute2dList, indicatorFilePath, indicatorDirsFilePath, 
#                                                    censusData.ObjIDInitialValue, censusData.FID_keyField, organicGrowthPCADir)
# analysis(organicGrowthIdsInRegion, censusData.dataPath, 
#          SVIDictList, RotatedCompMs, KMO_Bartletts, 
#          indicatorFilePath,  organicGrowthFinalOutputFile, corVLineRegionNum=OGBdRN)
# #---End---
 
#---Run Max-p---
# maxpOut = joinPath(outDir, 'MaxPResults')
maxpExp = MaxP(censusData, FID_keyField, constrainField, indicatorFilePath)
maxpExp.regionalization(constrainList, maxpRegionalizationDir)
maxpExp.getGrownRegion(constrainList, baseRegionFile, maxpRegionalizationDir)
maxpIdsInRegion = maxpExp.grownRegions
# Analysis max-p results
SVIDictList, RotatedCompMs, KMO_Bartletts = runPCA(maxpIdsInRegion, attribute2dList, indicatorFilePath, indicatorDirsFilePath, 
                                                   censusData.ObjIDInitialValue, censusData.FID_keyField, maxpPCADir)
analysis(maxpIdsInRegion, censusData.dataPath, SVIDictList, RotatedCompMs, KMO_Bartletts, indicatorFilePath, maxpFinalOutputFile)
#---End---




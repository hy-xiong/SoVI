import arcpy, os, sys
from time import time as T
from util import *
from util.miscellaneous import read1ColText

#-----Class of oragnic growth method-----
class OrganicGrowth:
    # Input: a number of feature in a Set, adjacency dictionary for each feature
    def __init__(self, featureDataset, adjDict, initialStartFeatureFile):
        # Temp working directory
        self._tempWorkDir = featureDataset._tempWorkDir
        arcpy.env.scratchWorkspace = self._tempWorkDir
        self._scratchFolder = os.path.join(self._tempWorkDir, 'scratch')
        clearDir(self._scratchFolder)
        # Spatial query key and initial value
        self._ObjIDField = featureDataset.ObjIDField
        self._ObjIDInitialValue = featureDataset.ObjIDInitialValue
        # Feature data path
        self._featurePath = featureDataset.dataPath
        # Feature data name
        self._featureName = featureDataset._featureName
        # Number of features
        self._featureAmount = featureDataset.num_of_features
        # Read the input start feature file to find their associated feature objectID in the feature dataset
        startFeaturesList = []  
        tempFieldName, tempIDs = read1ColText(initialStartFeatureFile, True)
        with arcpy.da.SearchCursor(self._featurePath, ["OID@", tempFieldName]) as cursor:  # @UndefinedVariable
            for row in cursor:
                if (row[1] in tempIDs):
                    startFeaturesList.append(row[0])
        # Offset the objID number if the objID number in the input feature dataset starts larger than 0
        if(self._ObjIDInitialValue > 0):
            for i in xrange(len(startFeaturesList)):
                startFeaturesList[i] -= self._ObjIDInitialValue
        # A list of layers in a grown region
        # Each layer contains all features surrounding current grown region
        # Each layer is a Set of features, feature's ID starts from 0
        self.layers = []
        self.layers.append(set(startFeaturesList))
        # Set other initial parameters
        self.currentLayer = 0
        self.increasedLayerNumberList = []
        self.featureCount = len(startFeaturesList)
        self.totalSet = set()
        self.totalSet = self.totalSet | self.layers[0]
        # Spatial Adjacency Dict
        self._adjDict = adjDict
        
    def grow(self, num = 0, layer = 0):
        if (num != 0):
            layer = 9999
        elif (num == 0 and layer == 0):
            layer = 9999
            num = self._featureAmount
        else:
            num = sys.maxsize
            pass
        print "\nStart Organic Growth from current layer " + str(self.currentLayer) + " (Number of features: " + str(self.featureCount) + ")."
        stTime = T()
        for layerIndex in xrange(self.currentLayer + 1, layer + 1):
            # Look for all first order neighbors of current grown region
            currentSet = self.layers[layerIndex - 1]
            if(currentSet == 0):
                baseLayerSet = currentSet
                excludedNeighborSet = currentSet
                newLayerSet = self._findNextFirstOrderNeighborLayer(baseLayerSet, excludedNeighborSet, self._adjDict)
            else:
                prevSet = self.layers[layerIndex - 2]
                baseLayerSet = currentSet
                excludedNeighborSet = currentSet | prevSet
                newLayerSet = self._findNextFirstOrderNeighborLayer(baseLayerSet, excludedNeighborSet, self._adjDict)
            # Update current layer number and total included features
            self.currentLayer += 1
            self.totalSet = self.totalSet | newLayerSet
            # Look for all feature that are not classified as first order neighbors of previous grown region but surrounded by the new grown regions
            enclosedFeatures = self._findWithinFeatures(self.totalSet)
            if(len(enclosedFeatures) > 0):
                self.totalSet = self.totalSet | enclosedFeatures
                newLayerSet = newLayerSet | enclosedFeatures
            # Add the newly found layer to the original region to grow it
            self.layers.append(newLayerSet)
            self.featureCount += len(newLayerSet)
            self.increasedLayerNumberList.append(len(newLayerSet))
            if(self.featureCount >= num):
                break
        print "\nOrganic Growth analysis has completed. " + "\nCurrent layer: " + str(self.currentLayer) + "\nNumber of features: " + str(self.featureCount) + "\nElapsed time: " + timer(stTime, T())
    
    def outputResultsToFeature(self):
        # Check field existence, if not exist, add it
        fieldName = "OGLayerId"
        field = arcpy.ListFields(self._featurePath, fieldName)
        if(len(field) == 0):
            arcpy.AddField_management(self._featurePath, fieldName, "SHORT", 4, "", "", fieldName, "NULLABLE")
        # Write layer number to feature dataset
        featureLayerDict = self._MapfeatureToLayer()
        with arcpy.da.UpdateCursor(self._featurePath, [self._ObjIDField, fieldName]) as cursor: #@UndefinedVariable
            for row in cursor:
                key = row[0] - self._ObjIDInitialValue
                if key in featureLayerDict:
                    row[1] = featureLayerDict[key]
                else:
                    row[1] = -1
                cursor.updateRow(row)
        print "\nRegions generated from organic growth has been output to the " + fieldName + " field in the feature " + self._featurePath
    
    def _findNextFirstOrderNeighborLayer(self, baseLayerSet, excludedNeighborSet,adjDict):
        newLayerSet = set()
        for featureIndex in baseLayerSet:
            for neighborIndex in adjDict[featureIndex]:
                if(neighborIndex not in excludedNeighborSet):
                    newLayerSet.add(neighborIndex)
        return newLayerSet
    
    def _findWithinFeatures(self, baseLayerSet):
        enclosedFeatures = set()
        tempLayer = os.path.join(self._scratchFolder, 'tempLayer' + str(self.currentLayer))
        tempUnionFeature = os.path.join(self._scratchFolder, 'tempUnion' + str(self.currentLayer) + ".shp")
        tempIdentityFeature = os.path.join(self._scratchFolder, 'tempIdentity' + str(self.currentLayer) + ".shp")
        # Uncheck Gap-allowed option in union tool to fill gaps (enclosed areas) in newly found first-order neighbor layer
        whereClause = '\"' + self._ObjIDField + '\" = ' + (' OR \"' + self._ObjIDField + '\" = ').join(str(i + self._ObjIDInitialValue) for i in baseLayerSet)
        arcpy.MakeFeatureLayer_management(self._featurePath, tempLayer, whereClause)
        arcpy.Union_analysis([tempLayer], tempUnionFeature, "ONLY_FID", None, "NO_GAPS")
        # Use identity tool to find enclosed features in newly found first-order neighbor layer
        arcpy.Identity_analysis(tempUnionFeature, self._featurePath, tempIdentityFeature, "ONLY_FID")
        # Field name of temp selected layer and orignal feature data in current generated identity feature
        ObjIDName_tempLayer = "FID_" + self._featureName[:6]
        if(self._featureName[5] == "_"):
            ObjIDName_originalFeature = "FID_" + self._featureName[:5] + "1"
        else:
            ObjIDName_originalFeature = "FID_" + self._featureName[:4] + "_1"
        with arcpy.da.SearchCursor(tempIdentityFeature, [ObjIDName_tempLayer, ObjIDName_originalFeature]) as cursor: #@UndefinedVariable
            for row in cursor:
                if(row[0] == -1 and row[1] != -1):
                    enclosedFeatures.add(row[1] - self._ObjIDInitialValue)
        #Remove scratched data
        arcpy.Delete_management(tempUnionFeature)
        arcpy.Delete_management(tempIdentityFeature)
        return enclosedFeatures
    
    def _MapfeatureToLayer(self):
        featureLayerDict = {}
        for layerNum in xrange(len(self.layers)):
            for featureID in self.layers[layerNum]:
                featureLayerDict[featureID] = layerNum
        return featureLayerDict                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
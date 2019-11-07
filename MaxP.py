'''
Created on Mar 17, 2016

@author: finix429
'''

import pysal as ps
import numpy as np
import pandas as pd
import os, arcpy, shutil
from ACS_Regionalization import *

from util import read1ColText
from collections import defaultdict
from util.miscellaneous import formatStr


class MaxP:
    '''
    classdocs
    '''
    
    def __init__(self, featureDataset, idField, countField, countFieldListFilePath):
        '''
        Constructor of MaxP. Note that isolated feature will not be joined with other features.
        @param featureDataset - instance of FeatureDataSet in data.py of util package
        @param idField - Primiary key field for spatial unit in the shapefile in featureDataset, not FID or ObjectID these kinds, should be sth. like GEOID etc.
        @param countField - single field from the shapefile in featureDataset which will be used as constrain in max-p
        @param countFieldListFilePath - paht to a file constaining a list of fields from the shapefile in featureDataset which will be used as parameters for max-p to compute SSD between regions
        '''
        index = []
        count_values = []
        targetCount_values = []
        srcFeatureSet = featureDataset.dataPath
        targetCountFields = read1ColText(countFieldListFilePath, False)
        # Make sure the feature dataset is a shapefile
        srcFeatureSetName = os.path.basename(srcFeatureSet)
        if '.shp' not in srcFeatureSetName:
            scratchFolder = os.path.join(os.path.dirname(os.path.dirname(srcFeatureSet)),'MaxPSrcData')
            if not os.path.isdir(scratchFolder):
                os.mkdir(scratchFolder)
            scratchShapefileName = '%s.shp' % srcFeatureSetName
            self._scratchShapefile = os.path.join(scratchFolder, scratchShapefileName)
            if arcpy.Exists(self._scratchShapefile):
                arcpy.Delete_management(self._scratchShapefile)
            arcpy.FeatureClassToFeatureClass_conversion(srcFeatureSet, scratchFolder, scratchShapefileName)
        else:
            self._scratchShapefile = srcFeatureSet
        # Generate pandas dataframe for ACS regionalization input
        self._fieldList = [idField, countField]
        self._fieldList.extend(targetCountFields)
        with arcpy.da.SearchCursor(self._scratchShapefile, self._fieldList) as cursor:  #@UndefinedVariable
            for row in cursor:
                index.append(row[0])
                count_values.append(row[1])
                targetCount_values.append(row[2:])
        self._count_pdFrame = pd.DataFrame(count_values, index, [countField])
        self._targetCount_pdFrame = pd.DataFrame(targetCount_values, index, targetCountFields)
        self._targetMOE_pdFrame = np.empty(self._targetCount_pdFrame.shape)
        self._targetMOE_pdFrame.fill(0.0)
        self._targetMOE_pdFrame = pd.DataFrame(self._targetMOE_pdFrame, list(self._targetCount_pdFrame.index.values), list(self._targetCount_pdFrame.columns.values))
        # Get pysal object from shapefile in featureDataset
        self._w = ps.rook_from_shapefile(self._scratchShapefile, idVariable = idField)
        self._shp = ps.open(self._scratchShapefile)
        
        
    def regionalization(self, constrains, outputFolder, clearPreviousResults = False):
        """
        Run Max-p regional based on a list of constrains, each constrain will run an independent Max-p .
        Result will be put in a specific location for re-use purpose.
        @param constrains: a list of constrains
        @param outputFolder: the location where to put regionalization results
        @param clearPreviousResults: if true, clear all results in the specified directory. If dataset is changed or constrain variable is changed, please input true.
        """
        print "Start Max-p for a series of constrains:\n%s " % ",".join(str(constrain) for constrain in constrains)
        if os.path.isdir(outputFolder):
            if clearPreviousResults:
                shutil.rmtree(outputFolder)
                os.mkdir(outputFolder)
        else:
            os.makedirs(outputFolder)
        for constrain in constrains:
            # Run ACS Regionalization as Max-p, use 0. MOE matrix to bypass the MOE evaluation procedure in Max-p
            outputFileSafeName = 'Constrain_%s.txt' % formatStr(constrain, 0)
            outputFile = os.path.join(outputFolder, outputFileSafeName)
            if not os.path.isfile(outputFile):
                print "Start regionalization for constrain %s" % str(constrain)
                np.random.seed(789)  # to ensure we get the same solution each time
                result = ACS_Regions(w=self._w,\
                                      target_est_count = self._targetCount_pdFrame.values,\
                                      target_moe_count = self._targetMOE_pdFrame.values,\
                                      count_est = self._count_pdFrame.values,\
                                      count_th_min = constrain,\
                                      target_th_all = 0.05,\
                                      compactness=self._shp,\
                                      pca = False)
                with open(outputFile, 'w') as fileWriter:
                    fileWriter.write('id in Region: %s\n' % str(result.region_ids))
                    fileWriter.write('regions: %s\n' % str(result.regions))
                    fileWriter.write('regions numbers: %s' % str(len(result.regions)))
                print 'regionalization for constrain', constrain, ' is finished'
                print 'Total time (h:m:s) used:','{0:d}:{1:d}:{2:.2f}'.format(int(divmod(result.time['total'], 3600)[0]), 
                                                                            int(divmod(divmod(result.time['total'], 3600)[1], 60)[0]),
                                                                            divmod(divmod(result.time['total'], 3600)[1], 60)[1])
                print 'Number of regions:', len(result.regions)
            else:
                print 'regionalization for constrain', constrain, 'has been done at\n%s' % outputFile
                with open(outputFile, 'r') as rd:
                    templines = rd.read().split('\n')
                    print templines[2]
        print "All Max-p are done"

    def getGrownRegion(self, constrainList, selectedFeaturesFile, regionalizationOutputFolder):
        """
        Find the grown regions based on the selected features in the selected features file, 
        and append the grown region from each regionalization result to the input shapefile.
        @pararm selectedFeatureFile: the path to the txt file constrain selected features. Only 1 column in this text with a column title in the first row.
        @param regionalizationOutputFolder: the path to the folder constain regionalization outputs
        """
        print "Start appending grown regions to %s as fields in this shapefile" % self._scratchShapefile
        # find the FID of every features
        fieldName, fieldValues = read1ColText(selectedFeaturesFile, True)
        selectedFeatureIDs = []
        grownRegionFieldNamePrefix = "MaxP"
        with arcpy.da.SearchCursor(self._scratchShapefile, [fieldName, "OID@"]) as cursor: #@UndefinedVariable
            for row in cursor:
                if row[0] in fieldValues:
                    selectedFeatureIDs.append(row[1])
        # remove previous attached max-p field in the shapefile
        existsFields = arcpy.ListFields(self._scratchShapefile, "SmallInteger")
        removeFields = []
        for field in existsFields:
            if grownRegionFieldNamePrefix in field:
                removeFields.append(field)
        if len(removeFields) > 0:
            arcpy.DeleteField_management(self._scratchShapefile, removeFields)
        # Read each regionalization results given each constrain
        self.outputConstrains = []
        constrain_idInRegions = defaultdict(list)
        constrain_Regions = defaultdict(list)
        constrainListSet = set(constrainList)
        for regionalizationFile in sorted(os.listdir(regionalizationOutputFolder), key = lambda s: float(s.split('_')[1].split('.')[0])):
            if int(regionalizationFile.split('.')[0][10:]) in constrainListSet:
                constrain = regionalizationFile.split('_')[1].split('.')[0]
                self.outputConstrains.append(constrain)
                with open(os.path.join(regionalizationOutputFolder, regionalizationFile), 'r') as fileReader:
                    idInRegions = fileReader.readline().split(': ')[1].split(', ')
                    idInRegions[0] = idInRegions[0].replace('[','')
                    idInRegions[-1] = idInRegions[-1].replace(']','')
                    idInRegions = [int(float(regionId)) for regionId in idInRegions]
                    regions = fileReader.readline().split(': ')[1].split('], [')
                    regions[0] = regions[0].replace('[','')
                    regions[-1] = regions[-1].replace(']','')
                    regions = [[int(float(featureId)) for featureId in idList.split(', ')] for idList in regions]
                    constrain_idInRegions[constrain] = idInRegions
                    constrain_Regions[constrain] = regions
        # attach max-p results and grown regions to the shapefile and get each Region constrained features
        self.grownRegions = []
        for i in xrange(len(self.outputConstrains) + 1):
            grownRegionFieldName = "MaxP_%d" % i
            regionalizationFieldName = "MaxPRg_%d" % i
            print grownRegionFieldName, regionalizationFieldName
            grownFeaturesList = []
            if i == 0:
                grownFeaturesList = [int(id) for id in selectedFeatureIDs]
            else:
                constrain = self.outputConstrains[i - 1]
                grownRegionList = []
                for selectedFeatureID in selectedFeatureIDs:
                    if constrain_idInRegions[constrain][selectedFeatureID] not in grownRegionList:
                        grownRegionList.append(constrain_idInRegions[constrain][selectedFeatureID])
                for regionID in grownRegionList:
                    grownFeaturesList.extend(constrain_Regions[constrain][regionID])
                arcpy.AddField_management(self._scratchShapefile, regionalizationFieldName, "SHORT")
                with arcpy.da.UpdateCursor(self._scratchShapefile, ["OID@", regionalizationFieldName]) as cursor: #@UndefinedVariable
                    for row in cursor:
                        row[1] = constrain_idInRegions[constrain][row[0]]
                        cursor.updateRow(row)
            arcpy.AddField_management(self._scratchShapefile, grownRegionFieldName, "SHORT")
            where_clause = '\"FID\" = ' + ' OR \"FID\" = '.join(str(id) for id in grownFeaturesList)
            with arcpy.da.UpdateCursor(self._scratchShapefile, grownRegionFieldName, where_clause) as cursor: #@UndefinedVariable
                for row in cursor:
                    row[0] = 1
                    cursor.updateRow(row)
            self.grownRegions.append(grownFeaturesList)
        print "Reigonalization results and grown region for the following order of constrains has been attached to shapefile at %s:\n%s" % (self._scratchShapefile, ', '.join(str(c) for c in self.outputConstrains))
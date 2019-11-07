'''
Created on Oct 26, 2015

@author: hxiong
'''
import arcpy, os, collections
from time import time as T
from miscellaneous import timer
from miscellaneous import read1ColText


#Class referring to a census feature dataset
class FeatureDataSet:
    def __init__(self, srcDataPath, dstDataPath, tempWorkDir, keyField = None, excludedFeatureDataFilePath = ""):
        #Path to feature data
        self._featureName = os.path.basename(dstDataPath)
        #Remove excluded features from dataset and generates a new dataset
        self._originalToNewIDMaps = {}
        if excludedFeatureDataFilePath != "":
            excludedFeatureIDField, excludedFeatureIDs = read1ColText(excludedFeatureDataFilePath, True)
            if(arcpy.Exists(dstDataPath)):
                arcpy.Delete_management(dstDataPath)
            if(type(excludedFeatureIDs[0]) == str):
                whereClause = 'NOT \"' + excludedFeatureIDField + '\" = ' + (' AND NOT \"' + excludedFeatureIDField + '\" = ').join("\'" + str(x) + "\'" for x in excludedFeatureIDs)
            else:
                whereClause = 'NOT \"' + excludedFeatureIDField + '\" = ' + (' AND NOT \"' + excludedFeatureIDField + '\" = ').join(str(x) for x in excludedFeatureIDs)
            arcpy.FeatureClassToFeatureClass_conversion(srcDataPath, os.path.dirname(dstDataPath), os.path.basename(dstDataPath), whereClause)
            print "New feature layer is generated with user specified features excluded.\nLocation at " + dstDataPath + "\n"
        else:
            arcpy.FeatureClassToFeatureClass_conversion(srcDataPath, os.path.dirname(dstDataPath), os.path.basename(dstDataPath))
            print "New feature layer is copied for analysis.\nLocation at " + dstDataPath + "\n" 
        self.dataPath = dstDataPath
        #FID_keyFied mapper
        if(keyField != None):
            self.FID_keyField = {}
            with arcpy.da.SearchCursor(self.dataPath, ["OID@", keyField]) as cursor: #@UndefinedVariable
                for row in cursor:
                    self.FID_keyField[row[0]] = row[1]
        print srcDataPath
        print self.FID_keyField
        #Count # of features in the feature dataset
        arcpy.MakeTableView_management(self.dataPath, "tempTableView")
        self.num_of_features = int(arcpy.GetCount_management("tempTableView").getOutput(0))
        #Get feature ObjectID (primary key)
        self.ObjIDField = arcpy.Describe(self.dataPath).OIDFieldName
        #Get initialValue of feature's ObjectID
        with arcpy.da.SearchCursor(self.dataPath, "OID@") as cursor: #@UndefinedVariable
            firstRow = cursor.next()
            self.ObjIDInitialValue = firstRow[0]
        #Temp working directory 
        self._tempWorkDir = tempWorkDir
        arcpy.env.scratchWorkspace = self._tempWorkDir
        self._scratchFolder = os.path.join(self._tempWorkDir, 'scratch')
        print "Intermediate files are generated at: " + arcpy.env.scratchFolder #@UndefinedVariable
    
    #Ask user which spatial info from feature dataset needs to be updated
    def checkFeatureInfoDerivation(self, fieldNamesFile, spatialContRel, outdir):
        print "Start derive spatial information from input feature dataset for sample size experiment."
        dataCmd = input("Which feature info do you want to update:\n0: Attribute values\n1: Spatial contiguity\n2: Both\n3: None\nIf you use new data, please input 2\nUser input: ")
        print '\n'
        attributeTableFilePath = ''
        spatialContiguityFilePath = ''
        fileList = [f for f in os.listdir(outdir) if os.path.isfile(os.path.join(outdir, f))]
        fieldNames = read1ColText(fieldNamesFile, False)
        for f in fileList:
            if("Attribute" in f):
                attributeTableFilePath = os.path.join(outdir, f)
            if("Contiguity" in f):
                spatialContiguityFilePath = os.path.join(outdir, f)
        if(dataCmd == 0):
            if(attributeTableFilePath != ''):
                os.remove(attributeTableFilePath)
            attributeTableFilePath = self._outfFieldsValue(fieldNames, outdir)
        elif(dataCmd == 1):
            if(spatialContiguityFilePath != ''):
                os.remove(spatialContiguityFilePath)
            spatialContiguityFilePath = self._outfSpatialContiguityMatrix(spatialContRel, outdir)
        elif(dataCmd == 2):
            if(attributeTableFilePath != ''):
                os.remove(attributeTableFilePath)
            attributeTableFilePath = self._outfFieldsValue(fieldNames, outdir)
            if(spatialContiguityFilePath != ''):
                os.remove(spatialContiguityFilePath)
            spatialContiguityFilePath = self._outfSpatialContiguityMatrix(spatialContRel, outdir)
        elif(dataCmd == 3):
            pass
        else:
            raise Exception("Invalid input")            
        return attributeTableFilePath, spatialContiguityFilePath
    
    #Get all values from specified fields
    def _outfFieldsValue(self, field_names, outdir):
        #Inital 2d array for PCA input, with all elements equal to 0
        #Rows = number of features,  Cols = fields used for PCA analysis
        print "Start collecting attributes' values from the table of " + self._featureName
        stTime = T()
        num_of_PCA_vars = len(field_names)
        AllData = [[0 for i in range(num_of_PCA_vars)] for i in range(self.num_of_features)]
        #Get fields value from the feature dataset
        with arcpy.da.SearchCursor(self.dataPath, field_names) as cursor: #@UndefinedVariable
            rowNum = 0
            for row in cursor:
                for colNum in xrange(num_of_PCA_vars):
                    AllData[rowNum][colNum] = row[colNum]
                rowNum += 1
        outfPath = os.path.join(outdir, "Attribute_Table.txt")
        outf = open(outfPath,'w')
        # -----output format-----
        # Number of features: XX
        # Fields: field1    field2    field3    field4    ...
        # feature1_field1_value    feature1_field2_value    feature1_field3_value    feature1_field4_value    ...
        # feature2_field1_value    feature2_field2_value    feature2_field3_value    feature2_field4_value    ...
        # ...
        outString = "Number of features: " + str(self.num_of_features) + "\n"
        outString += "Fields: " + '\t'.join(field_names) + '\n'
        outString += '\n'.join('\t'.join(str(x) for x in y) for y in AllData)
        outf.write(outString)
        outf.close()
        print "All required attributes' values are written to " + outfPath + ". Elapsed Time: " + timer(stTime, T()) + "\n"
        return outfPath
    
    #Get spatial contiguity matrix
    def _outfSpatialContiguityMatrix(self, spatialContRel, outdir):
        print "Start producing spatial contiguity matrix for " + self._featureName +".\nTotally " + str(self.num_of_features) + " number of features."
        stTimeEntireOp = T()
        baselyr = os.path.join(self._scratchFolder, "BaseLayer.lyr")
        arcpy.MakeFeatureLayer_management(self.dataPath, baselyr) 
        outString = ""
        for i in xrange(self.num_of_features):
            if (i % 100 == 0):
                stTime1Feature = T()
            selectedlyr = os.path.join(self._scratchFolder, "SelectedLayer_" + str(i) + ".lyr")
            arcpy.MakeFeatureLayer_management(self.dataPath, selectedlyr, '\"' + self.ObjIDField + '\" = ' + str(i + self.ObjIDInitialValue))
            arcpy.SelectLayerByLocation_management(baselyr, spatialContRel, selectedlyr)
            outString += str(i) + ": ";
            neighborFeatureList = []
            with arcpy.da.SearchCursor(baselyr, ("OID@")) as cursor: #@UndefinedVariable
                for row in cursor:
                    # row only contain 1 element, ObjectID
                    cKey = row[0] - self.ObjIDInitialValue
                    if(i != cKey):
                        neighborFeatureList.append(cKey)
                outString += ', '.join(str(x) for x in neighborFeatureList) + '\n'
            if( i % 100 == 99):
                print "100 features' neighbors are identified (total: " + str(i+1) + "). Elapsed Time: " + timer(stTime1Feature, T())
            arcpy.SelectLayerByAttribute_management(baselyr, "CLEAR_SELECTION")
        outfPath = os.path.join(outdir, "Spatial_Contiguity_" + spatialContRel + ".txt")
        outf = open(outfPath, 'w')
        # -----Output format------
        # 0 (feature id, defaultly starting from 0): neighbor1_id, neighbor2_id, ...
        # 1: neighbor1_id, neighbor2_id, ...
        # ...
        outf.write(outString)
        outf.close()
        print "Spatial Contiguity matrix has been written into " + outfPath + ". Elapsed time: " + timer(stTimeEntireOp, T()) + "\n"
        return outfPath

#Return number of feature, fields, and an numpy array from field value data file
def readFieldValueFile(fieldValueFilePath):
    f = open(fieldValueFilePath, 'r')
    numberOfFeature = f.readline()
    numberOfFeature = numberOfFeature.split(': ')[1]
    fields = f.readline()
    fields = fields.split(': ')[1]
    fields = fields.split()
    temp2dList = []
    for dataline in f:
        temp2dList.append([float(i) for i in dataline.split()])
    f.close()
    return numberOfFeature, fields, temp2dList

#Return a dictionary of spatial adjacency
def readAdjacencyFile(adjacencyFilePath):
    adjDict = collections.defaultdict(list)
    index = 0
    with open(adjacencyFilePath) as f:
        for line in iter(f.readline, ''):
            line = line.split(': ')
            if(line[1] != '\n'):
                lineSplit = line[1].split(', ')
                adjDict[index] = [int(x) for x in lineSplit]
            else:
                adjDict[index] = []
            index += 1
    return adjDict
        
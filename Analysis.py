'''
Created on May 30, 2016

@author: finix429
'''
import os, operator, shutil, copy
import numpy as np
from scipy import stats
from util import PCA, zscore, read1ColText, formatList, formatDict, plotTables, plotLineChart, plotBar

def runPCA(idsInRegion, featureAttributeList, indicatorFilePath, indicatorDirFilePath, initialFIDValue, FID_keyFieldDict, PCAResultDir):
    '''
    This function runs PCA for each region based on each region's feature IDs
    @param idsInRegion: a 2d list containing all features' IDs for each region (Note this id must starts from 0, so it might has offset comparing with the ids in its feature dataset)
    @param featureAttributeList: a 2d list containing all indicators' values
    @param indicatorFilePath: path to the file containing all indicator names
    @param indicatorDirFilePath: path to the file containing the social vulnerability direction for each indicator
    @param FID_keyFieldDict: a dict telling the primary key value based on the FID/ObjectID of a feature in its feature dataset
    @param PCAResultDir: the dir where to write all PCA analysis results
    @return SVIDictList: a list of dict with each dict containing the SVI value for each feature in the base region
    @return RotatedCompMs: a list of 2d list with each 2d list contains the rotated PCA solution for each region
    @return KMO_Bartletts: a 2d list with each sub list containing KMO and Bartletts Sign
    '''
    # Local vars
    baseRegionIDList = idsInRegion[0]
    numOfBaseRegionFeature = len(baseRegionIDList)
    indicatorList = read1ColText(indicatorFilePath, False)
    indicatorDirList = [int(i) for i in read1ColText(indicatorDirFilePath, False)]
    # Output vars
    SVIDictList = []
    RotatedCompMs = []
    KMO_Bartletts = []
    # Start PCA for every grown region in organic growth
    if os.path.isdir(PCAResultDir):
        shutil.rmtree(PCAResultDir)
    os.mkdir(PCAResultDir)
    print 'Start PCA for %d Regions' % len(idsInRegion)
    for regionId in xrange(len(idsInRegion)):
        PCAInput = []
        PCAInputIDList = []
        PCAInputFeatureIDList = []
        for id in idsInRegion[regionId]:
            rowRecord = copy.copy(featureAttributeList[id])
            for col in xrange(len(rowRecord)):
                rowRecord[col] *= indicatorDirList[col]
            PCAInput.append(rowRecord)
            PCAInputIDList.append(id)
            PCAInputFeatureIDList.append(id + initialFIDValue)
#         print "\n Start PCA for Region %d" % regionId
        # Run PCA in SPSS
        StandardizedPCAInput = zscore(PCAInput)
        CorM, NonPositiveCorM, KMO, Bartlett_sig, Comm, VarExpInfo, RotatedVarExpInfo, CompM, RotatedCompM, CompScoreCoeffM, CompScore  = PCA(StandardizedPCAInput, indicatorList, regionId)
        # Compute social vulnerability index for target study area
        baseRegionFeaturePosition = [PCAInputIDList.index(id) for id in baseRegionIDList]
        tempSVIDict = {}
        for i in xrange(len(baseRegionFeaturePosition)):
            SVIValue = sum(CompScore[baseRegionFeaturePosition[i]]) / len(CompScore[baseRegionFeaturePosition[i]])
            keyField = FID_keyFieldDict[baseRegionIDList[i] + initialFIDValue]
            tempSVIDict[keyField] = SVIValue
        # Collect KMO, Barlett Plot_test, rotated factor pattern, and SVI for target study area
        if(not NonPositiveCorM):
            KMO_Bartletts.append([KMO,Bartlett_sig])
        else:
            KMO_Bartletts.append(['N/A', 'N/A'])
        RotatedCompMs.append(RotatedCompM)
        SVIDictList.append(tempSVIDict)
        # Write PCA results
        tempOutFile = os.path.join(PCAResultDir, "Region" + str(regionId) + ".txt")
        featureIdList = [str(featureID) for featureID in PCAInputFeatureIDList]
        compNameList = ['Comp' + str(i + 1) for i in xrange(len(VarExpInfo[0]))]
        with open(tempOutFile, 'w') as f:
            f.write("All feauture IDs are the objectID/FID in the generated new feature dataset with selected features excluded (if there is feature excluded)\n\n")
            f.write(formatList(PCAInputFeatureIDList, 1, "Selected feature row/id number") + "\n\n")
            f.write(formatList(PCAInput, 2, "Raw PCA Input:", [indicatorList, featureIdList]) + "\n\n")
            f.write(formatList(CorM, 2, "Correlation Matrix", [indicatorList, indicatorList]) + '\n\n')
            if(not NonPositiveCorM):
                f.write('KMO: ' + str(KMO) + "\t Bartlett Sig.: " + str(Bartlett_sig) + '\n\n')
            else:
                f.write('This correlation matrix is not positively definite.\nOne the eigenvalues of the correlation matrix must be negative.\nKMO and Bartlett Plot_test cannot be performed\n\n')
            f.write(formatList(Comm, 1, "Communalities", [indicatorList]) + '\n\n')
            f.write(formatList(VarExpInfo, 2, "Varaince Explained of unrotated solution", [compNameList, ["Total", "% Variance", "Cummulative %"]]) + '\n\n')
            f.write(formatList(RotatedVarExpInfo, 2, "Varaince Explained of rotated solution", [compNameList, ["Total", "% Variance", "Cummulative %"]]) + '\n\n')
            f.write(formatList(CompM, 2, "Component loading Matrix", [compNameList, indicatorList]) + '\n\n')
            f.write(formatList(RotatedCompM, 2, "Rotated Component loading Matrix", [compNameList, indicatorList]) + '\n\n')
            f.write(formatList(CompScoreCoeffM, 2, "Component Score Coefficient Matrix", [compNameList, indicatorList]) + '\n\n')
            f.write(formatList(CompScore[:numOfBaseRegionFeature], 2, "Component Scores",[compNameList,featureIdList]) + '\n\n')        
            f.write(formatDict(tempSVIDict, "Social Vulnerability Indices"))
    print 'PCA for all regions (# = %d) are done.' % len(idsInRegion)
    return SVIDictList, RotatedCompMs, KMO_Bartletts

def analysis(idsInRegion, analysisFeatureFilePath, SVIDictList, RotatedCompMs, KMO_Bartletts, 
             indicatorFilePath, analysisResultFilePath, corVLineRegionNum=None):
    '''
    This function analysis the PCA results and writes the final output
    @param idsInRegion: a 2d list containing all features' IDs for each region (Note this id must starts from 0, so it might has offset comparing with the ids in its feature dataset)
    @param analysisFeatureFilePath: path to which feature is used for analysis
    @param SVIDictList: a list of dict with each dict containing the SVI value for each feature in the base region
    @param RotatedCompMs: a list of 2d list with each 2d list contains the rotated PCA solution for each region
    @param KMO_Bartletts: a 2d list with each sub list containing KMO and Bartletts Sign
    @param indicatorFilePath: path to the file containing all indicator names
    @param analysisResultFilePath: path to the file where to put all analysis results
#     @param fig_outDir: directory to output all figures and charts
#     @param figNamePrefix: the prefix of figure Name for each output fig & chart in this analysis
    '''
    indicatorList = read1ColText(indicatorFilePath, False)
    # create SVI value list with same order of feature ID
    idList = sorted(SVIDictList[0].keys())
    SVIValuesList = [[SVIDict[id] for id in idList] for SVIDict in SVIDictList]
    SVIValuesList_zscore = stats.zscore(SVIValuesList, ddof=1, axis=1).tolist()
    # Compute SVI stability, get spearman & pearson correlation and statistical significance (Absolute value + ranking value)
    SpearmanCorMatrixSize = len(SVIValuesList)
    SpearmanCorMatrix = [[1.0 for i in xrange(SpearmanCorMatrixSize)] for i in xrange(SpearmanCorMatrixSize)]
    SpearmanCorCheckMatrix = [['' for i in xrange(SpearmanCorMatrixSize)] for i in xrange(SpearmanCorMatrixSize)]
    PearsonCorMatrix = copy.deepcopy(SpearmanCorMatrix)
    PearsonCorCheckMatrix = copy.deepcopy(SpearmanCorCheckMatrix)
    p_crit = .05
    for i in xrange(SpearmanCorMatrixSize):
        for j in xrange(i, SpearmanCorMatrixSize):
            SpearmanCor = stats.spearmanr(SVIValuesList[i], SVIValuesList[j])
            SpearmanCorMatrix[i][j] = SpearmanCor[0]
            SpearmanCorMatrix[j][i] = SpearmanCor[0]
            if SpearmanCor[1] < p_crit:
                SpearmanCorCheckMatrix[i][j] = 'Y'
                SpearmanCorCheckMatrix[j][i] = 'Y'
            PearsonCor = stats.pearsonr(SVIValuesList_zscore[i], SVIValuesList_zscore[j])
            PearsonCorMatrix[i][j] = PearsonCor[0]
            PearsonCorMatrix[j][i] = PearsonCor[0]
            if PearsonCor[1] < p_crit:
                PearsonCorCheckMatrix[i][j] = 'Y'
                PearsonCorCheckMatrix[j][i] = 'Y'
    # Given a way of classifying absolute SoVI score, finding out how many census tracts stay in their classified values during the entire growing process
    classifyInterval = [-1.5, -0.5, 0.5, 1.5]
    classifyIntervalVisualText = ['=< -1.5', '-1.5 to -0.5', '-0.5 to 0.5', '0.5 to 1.5', '> 1.5']
    SVIValueClassified = [[-1 for j in xrange(len(SVIValuesList_zscore[i]))] for i in xrange(len(SVIValuesList_zscore))]
    stayInClassifyIntervalUnitCount = {i: 0 for i in xrange(len(classifyInterval)+1)}
    for i in xrange(len(SVIValuesList_zscore)):
        for j in xrange(len(SVIValuesList_zscore[i])):
            score_class = SVIValueClassified[i][j]
            for k in xrange(len(classifyInterval)):
                if SVIValuesList_zscore[i][j] <= classifyInterval[k]:
                    score_class = k
                    break
            if score_class == SVIValueClassified[i][j]:
                score_class = len(classifyInterval)
            SVIValueClassified[i][j] = score_class
    SVIValueClassified = np.array(SVIValueClassified).T
    for analysisUnitRow in SVIValueClassified:
        uniqueValues = np.unique(analysisUnitRow)
        if uniqueValues.shape[0] == 1:
            stayInClassifyIntervalUnitCount[uniqueValues[0]] += 1
    stayInClassifyIntervalUnitPct = {i: stayInClassifyIntervalUnitCount[i]*100.0/len(idList) for i in stayInClassifyIntervalUnitCount}
    # Compute mean and std dev. for each census unit's ranking value during sample size growth (ranking value)
    SVIRankValuesList = []
    for SVIValues in SVIValuesList:
        SVIRankValues = [operator.itemgetter(0)(t) for t in sorted(enumerate(SVIValues,1), key=operator.itemgetter(1), reverse = True)]
        SVIRankValuesList.append(SVIRankValues)
    SVIRankValuesList_np = np.array(SVIRankValuesList)
    SVIRankValuesList_stdevByCol = np.std(SVIRankValuesList_np, axis = 0, ddof = 1).tolist()
    SVIRankValuesList_meanByCol = np.mean(SVIRankValuesList_np, axis = 0).tolist()
    SVIRankValuesList_ZscoreByCol = stats.zscore(SVIRankValuesList_np, axis = 0, ddof = 1).tolist()
    SVI_baseRegionRank_zscore = SVIRankValuesList_ZscoreByCol[0]
    SVI_baseRegionRank = SVIRankValuesList_np[0].tolist()
    SVIs_deviationFromMean_visual = [SVIRankValuesList_meanByCol, SVIRankValuesList_stdevByCol, SVI_baseRegionRank, SVI_baseRegionRank_zscore]
    SVIs_deviationFromMean_visual = np.transpose(np.array(SVIs_deviationFromMean_visual)).tolist()
    # % of variables on 1st component, loading value > loadingCut(.5) is evaluated as the variable loads on the component
    loadingCut = .5
    pctVarsOn1stComponentList = []
    for RotatedCompM in RotatedCompMs:
        Counter_1stComp = 0
        for varLoadings in RotatedCompM:
            if (abs(varLoadings[0]) > loadingCut):
                Counter_1stComp += 1
        pctVarsOn1stComp = Counter_1stComp * 1.0 / len(indicatorList)
        pctVarsOn1stComponentList.append(pctVarsOn1stComp)
    # Average Component loading, loading > loadingCut is used to compute, 
    # Also count how many vars with loading > loadingCut in each component, variable with all loading < loadingCut and cross-loading var
    # Also shows which variables loads on which component
    AvgCompLoadings_NVarHigherThanLCPerComp = [] # 2d - row: per region solution; col: Comp # avgloading, Comp # N var > loadingCut, ...
    AvgCompLoadings = [] # 2d - row: per region solution; col: Comp #
    NVarHigherThanLCPerComp = [] # 2d - row: per region solution; col: Comp #
    NVarLessThanLC = [] # 1d - per region solution
    NVarCrossloading = [] # 1d - per region solution
    VarOnComp = [] # 2d - row: per region solution; col: variables; value: the component number this variable loads on and its direction
    VarOnCompChartData = [] # a list of 2d lists: each elementary list contains the which variable loads on which component
    maxCompN = 0
    for solution_index in xrange(len(RotatedCompMs)):
        RotatedCompM = RotatedCompMs[solution_index]
        tempCompN = len(RotatedCompM[0])
        if(maxCompN < tempCompN):
            maxCompN = tempCompN
        tempAvgCompLoading = [0. for x in xrange(tempCompN)]
        tempNVarHigherThanLC = [0 for x in xrange(tempCompN)]
        tempNVarLessThanLCCounter= 0
        tempCrossloadingCounter = 0
        tempVarOnComp = []
        tempVarOnCompChartData = []
        for variable_index in xrange(len(RotatedCompM)):
            varLoadings = RotatedCompM[variable_index]
            tempVarHigherThanLCList = []
            tempVarOnCompChartData_row = []
            for comp_index in xrange(len(varLoadings)):
                if(abs(varLoadings[comp_index]) > loadingCut):
                    tempNVarHigherThanLC[comp_index] += 1
                    tempAvgCompLoading[comp_index] += abs(varLoadings[comp_index])
                    tempVarHigherThanLCList.append('%d+' % (comp_index + 1) if varLoadings[comp_index] > 0 else '%d-' % (comp_index + 1))
                    tempVarOnCompChartData_row.append(1 if varLoadings[comp_index] > 0 else -1)
                else:
                    tempVarOnCompChartData_row.append(0)
            if(len(tempVarHigherThanLCList) < 1):
                tempNVarLessThanLCCounter += 1
            elif(len(tempVarHigherThanLCList) > 1):
                tempCrossloadingCounter += 1
            else:
                pass
            tempVarOnComp.append(', '.join(tempVarHigherThanLCList))
            tempVarOnCompChartData.append(tempVarOnCompChartData_row)
        for comp_id in xrange(tempCompN):
            if(tempNVarHigherThanLC[comp_id] > 0):
                tempAvgCompLoading[comp_id] = tempAvgCompLoading[comp_id] / tempNVarHigherThanLC[comp_id]
        AvgCompLoadings.append(tempAvgCompLoading)
        NVarHigherThanLCPerComp.append(tempNVarHigherThanLC)
        NVarLessThanLC.append(tempNVarLessThanLCCounter)
        NVarCrossloading.append(tempCrossloadingCounter)
        VarOnComp.append(tempVarOnComp)
        VarOnCompChartData.append(tempVarOnCompChartData)
    for solution_index in xrange(len(RotatedCompMs)):
        tempCompN = len(AvgCompLoadings[solution_index])
        if(tempCompN < maxCompN):
            AvgCompLoadings[solution_index].extend(['N\A'] * (maxCompN - tempCompN))
            NVarHigherThanLCPerComp[solution_index].extend(['N\A'] * (maxCompN - tempCompN))
    for solution_index in xrange(len(RotatedCompMs)):
        temp_AvgCompLoadings_NVarHigherThanLCPerComp = []
        for comp_id in xrange(maxCompN):
            temp_AvgCompLoadings_NVarHigherThanLCPerComp.append(AvgCompLoadings[solution_index][comp_id])
            temp_AvgCompLoadings_NVarHigherThanLCPerComp.append(NVarHigherThanLCPerComp[solution_index][comp_id])
        AvgCompLoadings_NVarHigherThanLCPerComp.append(temp_AvgCompLoadings_NVarHigherThanLCPerComp)
    # Membership of census unit in first quantile (10% most vulnerable) (ranking value)
    sortedSVIIDLists = [[id_SVI_Pair[0] for id_SVI_Pair in sortedSVIList] for sortedSVIList in (sorted(SVIDict.items(), key = operator.itemgetter(1), reverse = True) for SVIDict in SVIDictList)]
    censusUnitsNum = len(sortedSVIIDLists[0])
    quantilePct = 0.1
    firstQuantileCensusUnitsNum = int(censusUnitsNum * quantilePct)
    firstQuantileCensusUnitsList = []
    if(firstQuantileCensusUnitsNum > 0):
        for sortedSVI_index in sortedSVIIDLists:
            firstQuantileCensusUnitsList.append(sortedSVI_index[:firstQuantileCensusUnitsNum])
        firstQuantileCensusUnitsSimilarityMatrix = [[0. for i in xrange(len(firstQuantileCensusUnitsList))] for i in xrange(len(firstQuantileCensusUnitsList))]
        for i in xrange(len(firstQuantileCensusUnitsList)):
            for j in xrange(i, len(firstQuantileCensusUnitsList)):
                if(i == j):
                    firstQuantileCensusUnitsSimilarityMatrix[i][j] = 1.
                else:
                    SimilarCensusUnits = 0
                    for censusId_i in firstQuantileCensusUnitsList[i]:
                        for censusId_j in firstQuantileCensusUnitsList[j]:
                            if (censusId_i == censusId_j):
                                SimilarCensusUnits += 1
                                break
                    similarity = SimilarCensusUnits * 1.0 / firstQuantileCensusUnitsNum
                    if(similarity == 0.):
                        firstQuantileCensusUnitsSimilarityMatrix[i][j] = ''
                        firstQuantileCensusUnitsSimilarityMatrix[j][i] = ''
                    else:
                        firstQuantileCensusUnitsSimilarityMatrix[i][j] = similarity
                        firstQuantileCensusUnitsSimilarityMatrix[j][i] = similarity       
        firstQuantileMembershipDict = {}
        for firstQuantileCensusUnits in firstQuantileCensusUnitsList:
            for censusId in firstQuantileCensusUnits:
                if censusId in firstQuantileMembershipDict:
                    firstQuantileMembershipDict[censusId] += 1
                else:
                    firstQuantileMembershipDict[censusId] = 1
        for censusId, occuranceNum in firstQuantileMembershipDict.items():
            firstQuantileMembershipDict[censusId] = occuranceNum * 1.0 / len(firstQuantileCensusUnitsList)
        firstQuantileMembershipOrdered = sorted(firstQuantileMembershipDict.items(), key = operator.itemgetter(1), reverse = True)
    # Output Analysis Results
    if(os.path.isfile(analysisResultFilePath)):
        os.remove(analysisResultFilePath)
    # output txt & figures
    with open(analysisResultFilePath, 'w') as f:
        #axis
        regionIDs = [('Region%d' % i) for i in xrange(0,len(idsInRegion))]
        # general info and input data info     
        f.write("All feauture IDs are the objectID/FID in the generated new feature dataset with selected features excluded (if there is feature excluded)\n\n")
        f.write('Indicator List: %s\n\n' % ', '.join(indicatorList))
        f.write('Initial input dataset: %s\n\n' % analysisFeatureFilePath)
        # Each region's sample size + N:P ratio
        grownNum = [len(ids) for ids in idsInRegion]
        grownInc = [0]
        grownInc.extend([grownNum[i] - grownNum[i - 1] for i in xrange(1, len(grownNum))])
        featureN = len(indicatorList)
        pRatio = [ (n * 1.) / featureN for n in grownNum]
        f.write(formatList([grownInc, grownNum, pRatio] , 2, 'Sample size', [regionIDs, ["Inc", "Num", "N:P"]]) + '\n\n')
        # For every extracted component in each region: show # of variables having loading value higher than loading cut, and average loading of variables having loading value higher than loading cut
        f.write(formatList(AvgCompLoadings_NVarHigherThanLCPerComp, 2, 
                           "Avg loading of vars with loading > {0:.2f} & # of vars with loading > {0:.2f} per Comp".format(loadingCut), 
                           [[('Comp%d_AvgLoad' % (i/2 + 1)) if i % 2 == 0 else ('Comp%d_>%.2f_vars' % ((i-1)/2 + 1, loadingCut)) for i in xrange(2 * maxCompN)], regionIDs]) + '\n\n')
        # For every region, show # of variables having all loading values (on all components lower than loading cut), # of cross laoding variables, and KMO & Bartlett Plot_test
        NVarLessThanLC_Crossloading_KMO_Barrlett = [[0. for x in xrange(4)] for x in xrange(len(AvgCompLoadings_NVarHigherThanLCPerComp))]
        for solution_id in xrange(len(AvgCompLoadings_NVarHigherThanLCPerComp)):
            NVarLessThanLC_Crossloading_KMO_Barrlett[solution_id][0] = NVarLessThanLC[solution_id]
            NVarLessThanLC_Crossloading_KMO_Barrlett[solution_id][1] = NVarCrossloading[solution_id]
            NVarLessThanLC_Crossloading_KMO_Barrlett[solution_id][2] = KMO_Bartletts[solution_id][0]
            NVarLessThanLC_Crossloading_KMO_Barrlett[solution_id][3] = KMO_Bartletts[solution_id][1]
        f.write(formatList(NVarLessThanLC_Crossloading_KMO_Barrlett, 2, 
                           "# of vars < %.2f, # of vars cross loading, KMO, Barlett Plot_test" % loadingCut, 
                           [["N vars with all loading < %.2f" % loadingCut,"N vars cross loading", "KMO", "Bartlett Plot_test"], regionIDs]) + '\n\n')
        # For every region, show which variable loads on which component
        f.write(formatList(VarOnComp, 2, 
                           "Variable loading relationship with component (value is component number, +/- meaning loading value is positive/negative)", 
                           [indicatorList, regionIDs]) + '\n\n')
        # show spearman correlation between each region
        f.write(formatList(SpearmanCorMatrix, 2, 
                           "Spearman Correlation between each pair of SVIs derived from grown regions", 
                           [regionIDs, regionIDs]) +'\n\n')
        # show spearman correlation significance between each region
        f.write(formatList(SpearmanCorCheckMatrix, 2, 
                           "Sig. < .05 Spearman Correlation between each pair of SVIs derived from grown regions", 
                           [regionIDs, regionIDs]) +'\n\n')
        # show pearson correlation between each region
        f.write(formatList(PearsonCorMatrix, 2, 
                           "Pearson Correlation between each pair of SVIs derived from grown regions", 
                           [regionIDs, regionIDs]) +'\n\n')
        # show pearson correlation significance between each region
        f.write(formatList(PearsonCorCheckMatrix, 2, 
                           "Sig. < .05 Pearson Correlation between each pair of SVIs derived from grown regions", 
                           [regionIDs, regionIDs]) +'\n\n')
        # show how many census units stays in their classified intervals
        f.write(formatList([[stayInClassifyIntervalUnitCount[i] for i in xrange(len(classifyIntervalVisualText))],
                            [stayInClassifyIntervalUnitPct[i] for i in xrange(len(classifyIntervalVisualText))]], 2, 
                           "Census units stay in its own classified interval during the whole growing process", 
                           [classifyIntervalVisualText, ['# of census units', '% of census units']]) +'\n\n')
        # For each pair of region, show % of most vulnerable census units (first 10%) appear in both region (baseline is 10% most vulnerability census units, not all census units)
        f.write(formatList(firstQuantileCensusUnitsSimilarityMatrix, 2, 
                          "Percentage of same census units found in two SVIs' first quantile list (10% most vulnerable):\n" + 
                          "Number of first quantile census units: " + str(firstQuantileCensusUnitsNum) + "\nZero values are not displayed", 
                          [regionIDs, regionIDs], 1) +'\n\n')
        # Show how SVI rank varies across different regions, and how baseline deviate from the mean rank
        f.write(formatList(SVIs_deviationFromMean_visual, 2, 
                           "SVI rank deviation from mean during growth. (Grown regions=%d)" % len(idsInRegion), 
                           [["Mean rank during growth", "Std. dev. rank during growth", "Rank derived from base region", "Z-score of rank derived from base region"], idList]) + '\n\n')
        # Show how often each most vulnerable (first 10%) census unit are found as most vulnerable in all regions
        f.write("First quantile census units membership:\n%s\n\n" % '\n'.join('{0:s}: {1:.2f}%'.format(censusId_pct[0], censusId_pct[1] * 100) for censusId_pct in firstQuantileMembershipOrdered))
        f.write(formatList(pctVarsOn1stComponentList, 1, "Pct of variables on 1st component:", [regionIDs]))
#     # plot rotated component pattern for each region, and highlight values larger than loading cut for PCA component pattern visualization
#     binaryTables_higherThanLoadingCut = []
#     tableAxes = []
#     maxTableNX = []
#     maxTableNY = []
#     nTables = len(RotatedCompMs)
#     for rotatedCompM in RotatedCompMs:
#         binaryTables_higherThanLoadingCut.append([[1 if cmpLoad > loadingCut else 0 for cmpLoad in var] for var in rotatedCompM])
#         tableAxes.append([['Comp%d' % (i+1) for i in xrange(len(rotatedCompM[0]))], indicatorList])
#         maxTableNX.append(len(rotatedCompM))
#         maxTableNY.append(len(rotatedCompM[0]))
#     figSizeXY = max(maxTableNX) * max(maxTableNY) * len(RotatedCompMs) * 0.6
#     tableTitles = [('Region%d (N=%d, N:P=%.2f)' % (i, grownNum[i], pRatio[i])) for i in xrange(len(idsInRegion))]
#     figDest = '%s_rotatedComps.png' % analysisResultFilePath[:-4]
#     plotTables(binaryTables_higherThanLoadingCut, tableAxes, tableTitles, figDest, value_l = RotatedCompMs, figSize = (figSizeXY, figSizeXY))
    # return parameters for plotting
    # plot spearman and pearson correlation between the base region and all grown regions
    X = grownNum
    Ys = [SpearmanCorMatrix[0], PearsonCorMatrix[0]][:1]
    Ycolors = ['b', 'g'][:1]
    Ylabels = ['Spearman Correlation', 'Pearson Correlation'][:1]
    Ycolor_check = [SpearmanCorCheckMatrix[0], PearsonCorCheckMatrix[0]][:1]
    markerColorIndices = []
    for i in xrange(len(Ys)):
        colorIndex = {Ycolors[i]: [], 'r': []}
        for j in xrange(len(Ys[i])):
            if Ycolor_check[i][j] != 'Y':
                colorIndex['r'].append(j)
            else:
                colorIndex[Ycolors[i]].append(j)
        markerColorIndices.append(colorIndex)
    figDest = '%s_correlation.png' % analysisResultFilePath[:-4]
    figTitle = 'Spearman correlation and Pearson correlation of SoVI score between base region and all grown regions'
    if corVLineRegionNum:
        plotLineChart(X, Ys, Ylabels, Ycolors, markerColorIndices, ['Sample Size', 'Correlation'],figDest, ylim=[0., 1.], figTitle=figTitle, fontSize=15.0, verticleLineXs=[X[corVLineRegionNum]])
    else:
        plotLineChart(X, Ys, Ylabels, Ycolors, markerColorIndices, ['Sample Size', 'Correlation'],figDest, ylim=[0., 1.], figTitle=figTitle, fontSize=15.0)
    CorrelationPlotXY = [grownNum, SpearmanCorCheckMatrix[0], markerColorIndices]
    # plot how many census units stays in their own classified intervals
    XTickLabels = classifyIntervalVisualText
    yInd = sorted(stayInClassifyIntervalUnitPct.keys())
    Y = [stayInClassifyIntervalUnitPct[k] for k in yInd]
    fConv = lambda x: '%.2f%%' % x
    figDest = '%s_censusUnitsStayInClassifiedInterval.png' % analysisResultFilePath[:-4]
    barText = ['%s, %d' % (fConv(stayInClassifyIntervalUnitPct[k]), stayInClassifyIntervalUnitCount[k]) for k in yInd]
    plotBar(XTickLabels, [Y], [''], 'b', ['std to mean', '% of census unit in base region'], figDest, barText=[barText], ylim=[0.,100.], yTickFormatConv=fConv, 
            figTitle='Census Units having consistent classified vulnerability level among all region\'s results')
    CensusUnitsStayInIntervalsPlotParam = [classifyIntervalVisualText, stayInClassifyIntervalUnitPct]
    # plot census units stay in first quantile during entire growing process
    xIntervals = [0.2, 0.4, 0.6, 0.8, 1]
    xTickLabels = ['0%-20%', '20%-40%', '40%-60%', '60%-80%', '80%-100%']
    xTickCount = [[s, 0] for s in xTickLabels]
    for e in firstQuantileMembershipOrdered:
        for i in xrange(len(xIntervals)):
            if e[1] <= xIntervals[i]:
                xTickCount[i][1] += 1
                break
    Y = [e[1] for e in xTickCount]
    figDest = '%s_mostVulnerableCensusUnits.png' % analysisResultFilePath[:-4]
    barText = ['%d' % v for v in Y]
    plotBar(xTickLabels, [Y], [''], 'b',
            ['Occurrence rate (# of regions: %d)' % len(regionIDs), '# of census units (# of census units in base region: %d)' % len(idList)], 
            figDest, barText=[barText], ylim=[0, int(max(Y)*1.2)+1],
            figTitle='Occurrence rate of census units found as most vulnerable (top %.0f%%) in all region\'s results' % (quantilePct*100.0))
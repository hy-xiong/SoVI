'''
Created on Jan 29, 2016
 
@author: hxiong
'''
import spss, spssaux
import numpy as np
from scipy.stats import mstats
import warnings
 
def PCA(StandardizedPCAInput, varList, regionId):
    """ Use SPSS python api to perform PCA
         
        Arguments:
        PCAInput - 2d python list for PCA input
        varList - a list of variables for each columns in the PCA input
         
        Returns:
        CorrelationMatrix - Correlation matrix
        KMO - Kaiser-Mayer-Olkin value
        Bartlett_sig - Significance value of Bartlett's Sphericity Test
        Communalities - Communalities of extracted components
        VarExplainedInfo - Variance explained from unrotated solution, including absolute variance, % of variance, and cummulative %
        RotatedVarExplainedInfo - Rotated variance explained from unrotated solution, including absolute variance, % of variance, and cummulative %
        ComponentMatrix - Unrotated component loading matrix
        RotatedComponentMatrix - Rotated component loading matrix
        ComponentScoreCoefficientMatrix - Component score coefficient derived from rotated solution
        ComponentScore - Component score derived from score coefficient
    """
     
    # SPSS command & dataset setup 
    spss.Submit("NEW FILE")
    with spss.DataStep():
        datasetObj = spss.Dataset()
        for var in varList:
            datasetObj.varlist.append(var)
        for row in StandardizedPCAInput:
            datasetObj.cases.append(row)
    if regionId == 18:
        debugFileOutputDir = r'C:\Users\hxiong\Dropbox\Haoyi Vulnerability\Simulation\Hurricane_Sandy'
        np.savetxt(np.array(StandardizedPCAInput), r'%s\PCAInut_r%d' % regionId, fmt='%.7f')
    spssPCASyntax = """FACTOR 
    /VARIABLES {0}
    /MISSING LISTWISE 
    /ANALYSIS {0}
    /PRINT UNIVARIATE INITIAL CORRELATION KMO EXTRACTION ROTATION FSCORE 
    /CRITERIA MINEIGEN(1) ITERATE(25) 
    /EXTRACTION PC 
    /CRITERIA ITERATE(100) 
    /ROTATION VARIMAX 
    /SAVE REG(ALL) 
    /METHOD=CORRELATION.""".format(' '.join(varList))
    spss.SetOutput("off")
    varNum = len(varList)
    # Create XML output from SPSS
    tag = spssaux.CreateXMLOutput(spssPCASyntax, omsid = 'Factor Analysis')
    # Get correlation matrix
    CorrelationMatrix = spssaux.getValuesFromXmlWorkspace(tag, 'Correlation Matrix', cellAttrib="number")
    CorrelationMatrix = _spssOutputTableConversion(CorrelationMatrix, varNum, varNum)
    # Get KMO and Bartlett Plot_test sig.
    KMO_and_Bartlett = spssaux.getValuesFromXmlWorkspace(tag, 'KMO and Bartlett Test', cellAttrib= "number")
    KMO_and_Bartlett = _spssOutputTableConversion(KMO_and_Bartlett, 1) 
    NonpositiveDefiniteCorM = False
    KMO = 0.
    Bartlett_sig = 0.
    if(len(KMO_and_Bartlett) == 0):
        NonpositiveDefiniteCorM = True
    else:
        KMO = KMO_and_Bartlett[0]
        Bartlett_sig = KMO_and_Bartlett[3]
    # Get Communalities
    Communalities = spssaux.getValuesFromXmlWorkspace(tag, 'Communalities', colCategory = "Extraction", cellAttrib="number")
    Communalities = _spssOutputTableConversion(Communalities, 1)
    # Get variances explained in unrotated solution
    VarExplained = spss.EvaluateXPath(tag[0], "/outputTree", """//pivotTable//category[@text="Extraction Sums of Squared Loadings"]/dimension["Statistics"]/category[@text="Total"]/cell/@number""")
    PctVarExplained = spss.EvaluateXPath(tag[0], "/outputTree", """//pivotTable//category[@text="Extraction Sums of Squared Loadings"]/dimension["Statistics"]/category[@text="% of Variance"]/cell/@number""")
    CummulativePctVarExplained = spss.EvaluateXPath(tag[0], "/outputTree", """//pivotTable//category[@text="Extraction Sums of Squared Loadings"]/dimension["Statistics"]/category[@text="Cumulative %"]/cell/@number""")
    VarExplained = _spssOutputTableConversion(VarExplained, 1)
    PctVarExplained = _spssOutputTableConversion(PctVarExplained, 1)
    CummulativePctVarExplained = _spssOutputTableConversion(CummulativePctVarExplained, 1)
    VarExplainedInfo = [VarExplained, PctVarExplained, CummulativePctVarExplained]
    # Get variances explained in rotated solution
    RotatedVarExplained = spss.EvaluateXPath(tag[0], "/outputTree", """//pivotTable//category[@text="Rotation Sums of Squared Loadings"]/dimension["Statistics"]/category[@text="Total"]/cell/@number""")
    RotatedPctVarExplained = spss.EvaluateXPath(tag[0], "/outputTree", """//pivotTable//category[@text="Rotation Sums of Squared Loadings"]/dimension["Statistics"]/category[@text="% of Variance"]/cell/@number""")
    RotatedCummulativePctVarExplained = spss.EvaluateXPath(tag[0], "/outputTree", """//pivotTable//category[@text="Rotation Sums of Squared Loadings"]/dimension["Statistics"]/category[@text="Cumulative %"]/cell/@number""")
    RotatedVarExplained = _spssOutputTableConversion(RotatedVarExplained, 1)
    RotatedPctVarExplained = _spssOutputTableConversion(RotatedPctVarExplained, 1)
    RotatedCummulativePctVarExplained = _spssOutputTableConversion(RotatedCummulativePctVarExplained, 1)
    RotatedVarExplainedInfo = [RotatedVarExplained, RotatedPctVarExplained, RotatedCummulativePctVarExplained]
    # Get number of extracted components
    if(len(VarExplained) != len(RotatedVarExplained)):
        w = "Region %d: unrotated and rotated solution finds different number of component based on Kaiser Criterion." % regionId
        warnings.warn(w, RuntimeWarning)
    CompNum = len(VarExplained)
    ComponentScoreColumnIndex = [varNum + i  for i in xrange(CompNum)]
    # Get component matrix
    ComponentMatrix = spssaux.getValuesFromXmlWorkspace(tag, 'Factor Matrix', cellAttrib= "number")
    ComponentMatrix = _spssOutputTableConversion(ComponentMatrix, CompNum, varNum)
    # Get rotated component matrix
    RotatedComponentMatrix = spssaux.getValuesFromXmlWorkspace(tag, 'Rotated Factor Matrix', cellAttrib= "number")
    RotatedComponentMatrix = _spssOutputTableConversion(RotatedComponentMatrix, CompNum, varNum)
    # Get component score coefficient matrix
    ComponentScoreCoefficientMatrix = spssaux.getValuesFromXmlWorkspace(tag, 'Factor Score Coefficient Matrix', cellAttrib= "number")
    ComponentScoreCoefficientMatrix = _spssOutputTableConversion(ComponentScoreCoefficientMatrix, CompNum, varNum)
    # Get component score
    dataCursor = spss.Cursor(ComponentScoreColumnIndex)
    ComponentScore = dataCursor.fetchall()
    dataCursor.close()
    return CorrelationMatrix, NonpositiveDefiniteCorM, KMO, Bartlett_sig, Communalities, VarExplainedInfo, RotatedVarExplainedInfo, ComponentMatrix, RotatedComponentMatrix, ComponentScoreCoefficientMatrix, ComponentScore
 
def zscore(PCAInput):
    """Standardize all input variables (columns)"""
    PCAInput_np = np.array(PCAInput)
    Transposed_PCAInput_np = PCAInput_np.transpose()
    Standardized_Transposed_PCAInput_np = mstats.zscore(Transposed_PCAInput_np, 1, 1)
    Standardized_PCAInput_np = Standardized_Transposed_PCAInput_np.transpose()
    return Standardized_PCAInput_np.tolist()
 
def _spssOutputTableConversion(spssOutputTable, spssCols = 1, spssRows = 0):
    """Convert the 1d unicode value list from SPSS XML outputs into correction dimensional numeric value list"""
    for x in xrange(len(spssOutputTable)):
        spssOutputTable[x] = float(spssOutputTable[x])
    if(spssCols == 1):
        return spssOutputTable
    else:
        Converted2dTable = []
        for spssRow in xrange(spssRows):
            Converted2dTable.append(spssOutputTable[spssRow * spssCols : (spssRow + 1) * spssCols])
        return Converted2dTable
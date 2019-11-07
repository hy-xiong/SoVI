'''
Created on Oct 26, 2015

@author: hxiong
'''

import os, math
import numpy as np
from matplotlib import pyplot as plt
plt.ioff()

def timer(stTime, edTime):
    """Convert elapsed time from seconds to h:m:s format"""
    hour, remain = divmod(edTime - stTime, 3600)
    minute, second = divmod(remain, 60)
    return "{0:d}:{1:d}:{2:.2f}".format(int(hour),int(minute), second)

def plotTables(l, axis, titles, figDest, value_l = None, fontSize = 10.0, figSize = None):
    """
    plot a list of 2d list (classified values) into one single figure using multiple plot. 
    l - a list of 2d lists, each 2d list contains classified values for color visualization
    title - list of titles for each table
    axis - a list of axis for each table. For each table's label axis, axis[0] is column title, axis[1] is row title
    figDest - the destinate path to save figure
    value_l - the corresponding real value tables for tables in l
    """
    nTables = len(l)
    # find two integer that multiplies larger than nTables, but difference of those two integer should be less than 2
    for nrow in xrange(1, nTables):
        ncol = nTables * 1.0 / nrow
        if abs(ncol - nrow) <= 2:
            ncol = int(math.ceil(ncol))
            break
    # plot all tables in multiple plots manner
    plt.rcParams['font.size']=fontSize
    fig, axs = plt.subplots(nrow, ncol)
#     manager = plt.get_current_fig_manager()
#     manager.window.showMaximized()
    if figSize:
        fig.set_size_inches(figSize)
    for irow in xrange(nrow):
        for icol in xrange(ncol):
            axs[irow][icol].xaxis.set_visible(False)
            axs[irow][icol].yaxis.set_visible(False)
            axs[irow][icol].set_frame_on(False)
    for i in xrange(len(l)):
        i_x, i_y = divmod(i, ncol)
        axs[i_x][i_y].set_title(titles[i], fontsize = 12)
        cell_color = []
        for vRow in l[i]:
            cell_color_row = []
            for v in vRow:
                if v == 1:
                    cell_color_row.append('r')
                elif v == -1:
                    cell_color_row.append('b')
                else:
                    cell_color_row.append('w')
            cell_color.append(cell_color_row)
        if value_l != None:
            cell_text = [['%.2f' % v for v in row] for row in value_l[i]]
        else:
            cell_text = [['' for v in row] for row in l[i]]
        axs[i_x][i_y].table(cellText = cell_text, cellColours = cell_color, colLabels = axis[i][0], rowLabels = axis[i][1], loc = 'center')
#         axs[i_x][i_y].tables[0].set_fontsize(fontSize)
    fig.savefig(figDest)

def plotLineChart(X, Ys, Ylabels, Ycolors, markerColorIndices, axisLabels, figDest, ylim=None,
                  figTitle=None, figSize=(25, 15), markerSize=10.0, fontSize=20.0, verticleLineXs=None):
    """
    plot a line chart for multiple lines
    X - x-axis values
    Ys - Y-axis values for each line
    Ycolors - colors for each line
    markerColorIndices - a list of dicts. Each dict contains <key-color, value - list of point index> for each line
    axisLabels - [0] for xlabel, [1] for ylabel
    figDest - the destinate path to save figure
    """
#     plt.rc('axes', titlesize=fontSize)     # fontsize of the axes title
#     plt.rc('axes', labelsize=fontSize)    # fontsize of the x and y labels
#     plt.rc('xtick', labelsize=fontSize)    # fontsize of the tick labels
#     plt.rc('ytick', labelsize=fontSize)    # fontsize of the tick labels
#     plt.rc('legend', fontsize=fontSize)    # legend fontsize
    plt.rcParams['font.size']=fontSize
    fig, ax = plt.subplots(1, 1)
#     manager = plt.get_current_fig_manager()
#     manager.window.showMaximized()
    if figTitle:
        ax.set_title(figTitle)
    if figSize:
        fig.set_size_inches(figSize)
    for i in xrange(len(Ys)):
        ax.plot(X, Ys[i], '%s-' % Ycolors[i], label=Ylabels[i])
        for ptColor, ptIndex in markerColorIndices[i].iteritems():
            x = [X[j] for j in ptIndex]
            y = [Ys[i][j] for j in ptIndex]
            ax.plot(x, y, 'o', c=ptColor, markersize=markerSize, label='_nolegend_')
    if ylim:
        ax.set_ylim(*ylim)
        vlYlim = ylim
    else:
        vlYlim = ax.get_ylim()
    if verticleLineXs:
        for i in xrange(len(Ys)):
            ax.plot(np.repeat(verticleLineXs[i], 100), np.linspace(vlYlim[0], vlYlim[1], num=100), '%s--' % Ycolors[i], label='_nolegend_')
    ax.legend(loc=0)
    ax.set_xlabel(axisLabels[0])
    ax.set_ylabel(axisLabels[1])
    fig.savefig(figDest)

# import numpy as np
# X = range(1, 31)
# Ys = [np.random.random(30), np.random.random(30)]
# Ycolors = ['b', 'g']
# Ylabels = ['line1', 'line2']
# markerColorIndices = []
# for i in xrange(len(Ys)):
#     colorIndex = {Ycolors[i]: [], 'r': []}
#     for j in xrange(len(Ys[i])):
#         if Ys[i][j] < 0.5:
#             colorIndex['r'].append(j)
#         else:
#             colorIndex[Ycolors[i]].append(j)
#     markerColorIndices.append(colorIndex)
# figDest = r'C:\Users\hxiong\Desktop\test.png'
# figTitle = 'test'
# plotLineChart(X, Ys, Ylabels, Ycolors, markerColorIndices, ['x', 'y'], figDest, figTitle=figTitle, fontSize=20.0)

def plotBar(XTickLabels, Ys, Ylabels, Ycolors, axisLabels, figDest, 
            barText=None, ylim=None, yTickFormatConv=None, figTitle=None, figSize=(25, 15), fontSize=20.0):
    """
    plot a bar chart
    XTickLabels - x-axis tick labels
    Y - Y-axis values
    axisLabels - [0] for xlabel, [1] for ylabel, both coule be None
    figDest - the destinate path to save figure
    barText - a list of text needed to be displayed on the top of each bar
    """
    plt.rcParams['font.size']=fontSize
    fig, ax = plt.subplots(1,1)
#     manager = plt.get_current_fig_manager()
#     manager.window.showMaximized()
    if figTitle:
        ax.set_title(figTitle)
    if figSize:
        fig.set_size_inches(figSize)
    width = 0.4
    x =np.arange(1, len(XTickLabels)+1)
    rects = []
    for i in xrange(len(Ys)):
        rects.append(ax.bar(x + width*i, Ys[i], width, color=Ycolors[i], label=Ylabels[i]))
#     ax.set_xticks(x + width/2.0 * len(Ys))
    ax.set_xticks(x)
    ax.set_xticklabels(XTickLabels)
    if axisLabels[0] != None:
        ax.set_xlabel(axisLabels[0]) 
    if axisLabels[1] != None:
        ax.set_ylabel(axisLabels[1])
    if ylim:
        ax.set_ylim(*ylim)
    if yTickFormatConv:
        yticks = ax.get_yticks()
        ax.set_yticklabels([yTickFormatConv(v) for v in yticks])
    if barText:
        for i in xrange(len(rects)):
            for j in xrange(len(rects[i])):
                ax.text(rects[i][j].get_x()+width/2.0, 1.05*rects[i][j].get_height(), barText[i][j], ha='center', va='bottom')
    fig.savefig(figDest)

# XTickLabels = ['a1a', 'b1b', 'c1c']
# Ys = [[0.12, 0.05, 0.31], [0.3, 0.4, 0.3]]
# Ylabels = ['t1', 't2']
# Ycolors = ['b', 'r']
# plotBar(XTickLabels, Ys, Ylabels, Ycolors, ['x', 'y'], 
#         r'C:\Users\hxiong\Desktop\test.png', 
#         barText=[XTickLabels]*len(Ys), yTickFormatConv=lambda x: '%.2f%%' % (x*100.0), ylim=[0,1], figTitle='test')

def formatList(l, dim, title = None, axis = [], precision = 3):
    """
    Format 1d and 2d List for print and write
    l - input list
    dim - total number of dimensions of input list
    title - list title
    axis - a list of label axis. For 1d list, axis[0] is the column Title. For 2d list, axis[0] is column title, axis[1] is row title
    precision - precision number for decimal values
    """
    # Display parameter
    spaceAmongCols = 2
    # Set title
    if (title == None):
        out_s = "The {0}d array is:\n".format(dim)
    else:
        out_s = "{0}:\n".format(title)
    if(axis and len(axis) != dim):
        raise Exception("Error, number of label axis does not match input list.\nNumber of label axis = {0}.\nInput list dimensions={1}".format(len(axis), dim))
    # Get data width for each row
    if(dim == 1):
        out_strList = [[formatStr(v, precision) for v in l]]
        if(axis):
            out_strList.insert(0, [str(s) for s in axis[0]])
    elif(dim == 2):
        out_strList = [[formatStr(v, precision) for v in l_row] for l_row in l]
        if(axis):
            out_strList_firsRow = [str(s) for s in axis[0]]
            out_strList_firsRow.insert(0, '')
            out_strList.insert(0, out_strList_firsRow)
            for rowIndex in xrange(1, len(out_strList)):
                out_strList[rowIndex].insert(0, str(axis[1][rowIndex - 1]))
    else:
        raise Exception("Cannot print {0} dimensional list".format(dim))  
    # Compute column width for every entire column
    out_strList_T = map(list, zip(*out_strList))
    colWidth = [max(len(s) for s in out_strList_T_row) for out_strList_T_row in out_strList_T]
    # Generate String output
    out_s += '\n'.join((' ' * spaceAmongCols).join('{0:>{width}s}'.format(out_strList_row[i], width = colWidth[i]) for i in xrange(len(out_strList_row))) for out_strList_row in out_strList)
    return out_s

def formatDict(d, title = None):
    """Format dict for print and write"""
    if (title == None):
        s = "The dict is:\n"
    else:
        s = title + ":\n"
    
    s += '\n'.join(str(key) + ": " + str(item) for key, item in d.iteritems())
    return s


def clearDir (directory):
    """ Clear all files in a directory """
    for f in os.listdir(directory):
        filePath = os.path.join(directory, f)
        os.unlink(filePath)

def read1ColText(f, firstRowTitle):
    """Read text file with only single column of data
        f - file path
        firstRowTitle - whether the first row is the column title or not
    """
    with open(f, 'r') as rf:
        content = rf.read()
    lines = content.split('\n')
    lines = [i for i in lines if i != ""]
    if(firstRowTitle):
        return lines[0], lines[1:]
    else:
        return lines

def formatStr(value, precision, width = 0):
    """Format input based on its type and value (right aligned)
        value - function input, could be int, float or str
        precision - the desired precision for float input
        width - the width of the entire return string
    """
    if('int' in str(type(value))):
        base = 'd'
    elif('str' in str(type(value))):
        base = 's'
    else:
        base = 'f'
    if(base == 'f'):
        precision = '.{0}'.format(precision)
    else:
        precision = ''
    if(width == 0):
        width = ''
    return '{0:>{width}{precision}{base}}'.format(value, width = width, precision = precision, base = base)

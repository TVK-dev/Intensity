import numpy as np
import pandas as pd


def list_files(directory):
    lp=[]
    ln=[]
    for file in os.listdir(directory):    
        if file.endswith(".csv"):
            lp.append(directory + file)
            namefile = file.replace('.csv','')
            ln.append(namefile)
    return lp, ln
    
def csv_to_lists(fn):
    Intensity = []
    Borders = []
    Type = []
    with open(fn, 'r') as f:
        header = f.readline()
        for line in f:
            stringfile = line.replace(',','.') 
            if stringfile.find(';') < 0:
                a = stringfile
                b = '0'
                c = '0'
            else:    
                a, b, c = stringfile.split(';')            
            Intensity.append(float(a))
            Borders.append(int(b))
            Type.append(int(c))
    return Intensity, Borders, Type

def cut_events(Intencity, Borders, Typ, delta = 0):
    listevents = []
    lborbers = []
    for i in range(len(Intencity)):
        if Typ[i] == 4:
            levent = []
            lborder = []
            for j in range(i, 0, -1):
                if Borders[j] == 2:
                    break
            for m in range(i, len(Intencity)):
                if Borders[m] == 3:
                    break
            j = j - delta
            m = m + delta
            if j < 0:
                j = 0
            if m > len(Intencity) - 1:
                m = len(Intencity) - 1
            for k in range(j, m + 1):
                levent.append(Intencity[k])
                lborder.append(Borders[k])
            listevents.append(levent) 
            lborbers.append(lborder) 
    return listevents, lborbers


def find_max(MasDat, MA):  
    j = 0
    i = 0
    numPoint = len(MasDat)
    LocalMaxIndex = []
    LocalMax = []            
    
    while j < numPoint-1 and i < numPoint-1:
        for i in range(j + 1, numPoint):
            if MasDat[i] > MA[i]:
                break

        for j in range(i + 1, numPoint):
            if  MasDat[j] < MA[j]:
                break
        maxm = -10000000#
        for k in range(i, j + 1):
            if MasDat[k] > maxm:
                maxm = MasDat[k]
                kk = k

        LocalMaxIndex.append(kk)
        LocalMax.append(MasDat[kk]) 
    return  LocalMaxIndex, LocalMax
         

def find_min(MasDat, MA):      
    numPoint = len(MasDat)
    LocalMinIndex = []
    LocalMin = []            
    LocalMinIndex.append(0)
    LocalMin.append(MasDat[0])            
    j = 0
    i = 0
    while j < numPoint-1 and i < numPoint-1:
        for i in range(j + 1, numPoint):
            if MasDat[i] < MA[i]:
                break
        for j in range(i + 1, numPoint):
            if MasDat[j] > MA[j]:
                break
        mmin = 10000000#
        for k in range(i, j + 1):
            if MasDat[k] < mmin: 
                mmin = MasDat[k]
                kk = k
              
        LocalMinIndex.append(kk)
        LocalMin.append(MasDat[kk])  
    LocalMinIndex.append(numPoint-1)
    LocalMin.append(MasDat[numPoint-1])                
    return  LocalMinIndex, LocalMin
    
    
def find_baseline(MasDat, ListMin, ListMinIndex):
    Baseline = [0 for i in range(len(MasDat))]
    for i in range(1, len(ListMin)):
        x1 = ListMinIndex[i - 1]
        x2 = ListMinIndex[i]
        y1 = ListMin[i - 1]
        y2 = ListMin[i]
        
        if x2 - x1 != 0:       
            for j in range(x1, x2 + 1):
                Baseline[j] = y1 + (j - x1) * (y2 - y1) / (x2 - x1)    
    return  Baseline


def find_MA(MasDat, koefMA = 0.9):
    MA = []    
    nextvalue = MasDat[0]
    for i in range(1, len(MasDat) + 1):
        nextvalue = nextvalue * koefMA + (1 - koefMA) * MasDat[i - 1]
        MA.append(nextvalue) 
    return MA    

def find_MAb(MasDat, koefMA = 0.9):
    MA = [0 for i in range(len(MasDat))]    
    nextvalue = MasDat[len(MasDat) - 1]
    for i in range(len(MasDat) - 1, -1, -1):
        nextvalue = nextvalue * koefMA + (1 - koefMA) * MasDat[i]
        MA[i] = nextvalue 
    return MA    

def calc_perc(MasDat, Baseline, ListMax):
    listAbs = []
    listPerc = []
    for i in range(len(ListMax)):
        ab = MasDat[ListMax[i]] - Baseline[ListMax[i]]
        #print(MasDat[ListMax[i]])
        if Baseline[ListMax[i]] != 0:
            perc = (MasDat[ListMax[i]] - Baseline[ListMax[i]])/Baseline[ListMax[i]]
        listAbs.append(ab)
        listPerc.append(100*perc)
    return listAbs, listPerc 

def calc_borders(Intensity):
    borders = []
    MA = find_MA(Intensity)
    MAb = find_MAb(Intensity)    
    LocalMinIndex, LocalMin = find_min(Intensity, MA)
    Baseline = find_baseline(Intensity, LocalMin, LocalMinIndex)
    LocalMaxIndex, LocalMax = find_max(Intensity, MA)  
    for i in range(len(LocalMaxIndex)):
        j =  LocalMaxIndex[i] - 1
        border = []
        while MA[j] < Intensity[j]:
            j = j - 1
            if j < 1:
                break
        border.append(j) 
        border.append(LocalMaxIndex[i])
        j =  LocalMaxIndex[i] + 1
        while MAb[j] < Intensity[j]:
            j = j + 1
            if j > len(Intensity) - 1:
                break
        border.append(j) 
        borders.append(border)
    return borders 


def file_to_events(filename):
    Intensity, BordersRazm, Type = csv_to_lists(filename) 
    MA = find_MA(Intensity)
    LocalMinIndex, LocalMin = find_min(Intensity, MA)
    Baseline = find_baseline(Intensity, LocalMin, LocalMinIndex)
    LocalMaxIndex, LocalMax = find_max(Intensity, MA)
    listAbs, listPerc = calc_perc(Intensity, Baseline, LocalMaxIndex)   
    borders = calc_borders(Intensity)
    return LocalMaxIndex, listAbs, listPerc, borders

def fill_graph(l):
    ll = 0
    for i in range(len(l)):
        if l[i] != 0:
            ll = l[i]
        else:
            l[i] = ll
    return l

def calc_phase(y1, y2):
    MA = intproc.find_MA(y1)
    LocalMinIndex1, LocalMin1 = intproc.find_min(y1, MA)
    LocalMaxIndex1, LocalMax1 = intproc.find_max(y1, MA)
    MA = intproc.find_MA(y2)
    LocalMinIndex2, LocalMin2 = intproc.find_min(y2, MA)
    LocalMaxIndex2, LocalMax2 = intproc.find_max(y2, MA)

    lph = [0]*len(y1)
    lf = [0]*len(y1)
    ldeg = []
    lfreq = []
    jj = 0
    jj1 = 0
    for i in range(min(len(LocalMaxIndex2), len(LocalMaxIndex1))):
    # цена градуса по опорной кривой
        if LocalMaxIndex1[i] - LocalMaxIndex1[i + 1] != 0:
            degpr  = 360/(LocalMaxIndex1[i] - LocalMaxIndex1[i + 1])
    # отставание/опережение
        mmin1 = 10e6  #минимум слева
        mmin2 = 10e6
        for j in range(len(LocalMaxIndex2)):
            dist =  abs(LocalMaxIndex2[j] - LocalMaxIndex1[i])
            if dist < mmin1 and LocalMaxIndex2[j] > LocalMaxIndex1[i]:    
                mmin1 =  dist
                jj = j
            if dist < mmin2 and LocalMaxIndex2[j] <= LocalMaxIndex1[i]:    
                mmin2 =  dist
                jj1 = j
    
        if abs(LocalMaxIndex2[jj] - LocalMaxIndex1[i])<abs(LocalMaxIndex2[jj1] - LocalMaxIndex1[i]):    
           deg = (LocalMaxIndex2[jj] - LocalMaxIndex1[i])*degpr
        else:
           deg = (LocalMaxIndex2[jj1] - LocalMaxIndex1[i])*degpr 
        
        if LocalMaxIndex2[jj] - LocalMinIndex2[jj+1] !=0:
             freq = abs((LocalMaxIndex1[i] - LocalMinIndex1[i + 1])/(LocalMaxIndex2[jj] - LocalMinIndex2[jj + 1]))
        ldeg.append(deg)
        lfreq.append(freq)
        lph[LocalMaxIndex2[i]] = deg
        lf[LocalMaxIndex2[i]] = freq
    
    lph = fill_graph(lph)        
    lf = fill_graph(lf) 
    return lf, lph


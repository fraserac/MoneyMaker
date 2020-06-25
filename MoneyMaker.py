# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 17:57:07 2020

@author: Wolowizard
MoneyMaker
"""
import numpy as np
import itertools as it
import sys
import math


title = "BrentFord (A) team vs West Brom"
dictWin = {"Sky bet, W1" : 50/20, "VirginBet, W2,":70/5, "William Hill, W3": 17/5}
dictLose = {"Sky bet, L1" : 190/10, "VirginBet, L2,":1/10, "William Hill, L3": 21/10}
dictDraw = {"Sky bet, D1" : 32/10, "VirginBet, D2,":19/4, "William Hill, D3": 11/5}

maxStake = 150.00#
res = maxStake/100


resVec = np.linspace(0, maxStake, 20)

winVec = [dictWin["Sky bet, W1"], dictWin["VirginBet, W2,"], dictWin["William Hill, W3"]]   #Contains odds in the form of coefficient of stake
loseVec = [dictLose["Sky bet, L1"], dictLose["VirginBet, L2,"], dictLose["William Hill, L3"]]
drawVec =[dictDraw["Sky bet, D1"], dictDraw["VirginBet, D2,"], dictDraw["William Hill, D3"]]
Si = 0.0
Sj = 0.0
Sk = 0.0
permW = 0.0
permL = 0.0
permD = 0.0
X=[] # vector of earnings
outputDict = {}
bigM = np.zeros(len(winVec)**3, dtype=object).reshape(len(winVec), len(winVec)**2)

permList = list((it.product('012', repeat=3)))

for ii in range(len(permList)):
    permW = winVec[int(permList[ii][0])]# 
    permL = loseVec[int(permList[ii][1])]
    permD = drawVec[int(permList[ii][2])]
    
    
    for jj in range(len(winVec)**2):
        for kk in range(len(winVec)):
            bigM[kk, jj] = -1*np.ones(len(winVec)**2).reshape(len(winVec), len(winVec))
            bigM[kk, jj][0,0] = permW
            bigM[kk, jj][1,1] = permL
            bigM[kk, jj][2,2] = permD
    
            

sep = " " 

for i in resVec:
    for j in resVec:
        for k in resVec:
            Si = i
            Sj = j 
            Sk = k
            Smat = [Si, Sj, Sk]
            for h in range(len(winVec)):
                for g in range(len(winVec)**2):
                    X = np.matmul(bigM[h,g], Smat)
            
            for vals in range (len(X)):
                if(all(aaa >= -100 for aaa in X)):
                    d = [h,g]
                    e = sep.join([str(elem) for elem in d])
                    outputDict[e] = [Smat, X, np.linalg.norm(X)]
        
        
curVal = 0.0 
curInd = [0,0]       

for keys in outputDict:
    if outputDict[keys][2] > curVal:
        curVal = outputDict[keys][2]
        curInd = [keys[0], keys[2]]

ind = sep.join([str(eleme) for eleme in curInd])       
indNo = [int(ind[0]), int(ind[2])]


if ind == '0 0':
    print("No viable candidate", ind)
    sys.exit()



print("Winning odds plus optimal stake:", bigM[indNo[0], indNo[1]][0,0] , "£", outputDict[ind][0][0]) # last term is S vector
print("Losing odds plus optimal stake:", bigM[indNo[0], indNo[1]][1,1], "£",outputDict[ind][0][1])             
print("Drawing odds plus optimal stake:", bigM[indNo[0], indNo[1]][2,2], "£",outputDict[ind][0][2] )

valListWin = list(dictWin.values())
valListLose = list(dictLose.values())
valListDraw = list(dictDraw.values())
keyListWin = list(dictWin.keys())
keyListLose = list(dictLose.keys())
keyListDraw = list(dictDraw.keys())


for item  in valListWin:
        if math.isclose(item, bigM[indNo[0], indNo[1]][0,0], rel_tol= 1e-1):
            winComp = keyListWin[valListWin.index(item)]
            
for item  in valListLose:
        if math.isclose(item, bigM[indNo[0], indNo[1]][1,1], rel_tol= 1e-1):
            loseComp = keyListLose[valListLose.index(item)]

for item  in valListDraw:
        if math.isclose(item,bigM[indNo[0], indNo[1]][2,2], rel_tol= 1e-1):
            drawComp = keyListDraw[valListDraw.index(item)]

print("Companies for win, loss, draw bets for game ", title, ": ", winComp, loseComp, drawComp )
print("Net loss/gain: ", "£", X-np.sum(outputDict[ind][0]))
            
            


#iterate through each output element and find vector with max magnitude and
#associated data. 
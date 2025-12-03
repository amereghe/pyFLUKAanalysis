'''
A. Mereghetti, 16/01/2023
A dirty script to crunch data saved in fort.23
'''

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import xlsxwriter
import pylandau
import pickle
import gc
from copy import deepcopy

class BeamPart():
    def __init__(self,partName,Ek,lPerNucleon=False):
        self.partName=None; self.mm=None; self.AA=None; self.ZZ=None
        self.pc=None; self.betaGammaRel=None; self.betaRel=None; self.gammaRel=None; self.Ek=None
        #
        self.setPartProp(partName)
        self.setPartEnergy(Ek,lPerNucleon)
        self.echoProperties()

    def setPartProp(self,partName):
        print("Identifying beam particle...")
        if (partName.upper().startswith("P")):
            self.mm=938.2720813 # p (FLUKA value) [MeV/c2]
            self.AA=1; self.ZZ=1
        elif (partName.upper().startswith("C")):
            self.mm=11.1748642E+03 # C (FLUKA value) [MeV/c2]
            self.AA=12; self.ZZ=6
        else:
            print("...unknown particle %s!"%(partName))
            exit(1)
        self.partName=partName

    def setPartEnergy(self,Ek,lPerNucleon=False):
        print("Setting energy properties of beam particle...")
        self.Ek=Ek
        if (lPerNucleon): self.Ek=self.Ek*self.AA
        self.pc=np.sqrt(self.Ek*(self.Ek+2*self.mm))
        self.betaGammaRel=self.pc/self.mm
        self.gammaRel=self.Ek/self.mm+1
        self.betaRel=self.betaGammaRel/self.gammaRel
        
    def echoProperties(self):
        if (self.partName.upper().startswith("P")):
            print("Beam particle: PROTON!")
        elif (self.partName.upper().startswith("C")):
            print("Beam particle: C-12!")
        print("""...beam data:
* m=%g MeV/c2; pc=%g MeV/c; Ek=%g MeV;
* A,Z=%d,%d;
* betaRel,gammaRel=%g,%g;
"""%(self.mm,self.pc,self.Ek,self.AA,self.ZZ,self.betaRel,self.gammaRel))

class Fort23Data():
    def __init__(self):
        '''FLUKA units'''
        self.IJ=np.array([])
        self.XX=np.array([]) # [cm]
        self.YY=np.array([]) # [cm]
        self.ZZ=np.array([]) # [cm]
        self.TXX=np.array([]) # []
        self.TYY=np.array([]) # []
        self.TZZ=np.array([]) # []
        self.ETRACK=np.array([]) # total energy [GeV]
        self.PTRACK=np.array([]) # momentum [GeV/c]
        self.XP=np.array([]) # [rad]
        self.YP=np.array([]) # [rad]
        self.DE=np.array([]) # [MeV]
        self.DP=np.array([]) # [MeV/c]
        self.TH=np.array([]) # [rad]
        self.DEFL=np.array([]) # delta total energy [GeV]
        self.DPFL=np.array([]) # delta momentum [GeV/c]
        self.DTFL=np.array([]) # [rad]
        self.sFacts={} # scaling factors []
        self.sUnits={} # units

    def __len__(self):
        return len(self.IJ)

#     def __del__(self):
    def DelAll(self):
        del self.IJ
        del self.XX
        del self.YY
        del self.ZZ
        del self.TXX
        del self.TYY
        del self.TZZ
        del self.ETRACK
        del self.PTRACK
        del self.XP
        del self.YP
        del self.DE
        del self.DP
        del self.TH
        del self.DEFL
        del self.DPFL
        del self.DTFL

    def setFacts(self,facts={},units={},lDebug=True):
        kFs=np.unique(list(facts.keys())+list(units.keys()))
        for kF in kFs:
            if (lDebug):
                print("...setting factor %E to %s with unit %s;"%(facts[kF],kF,units[kF]))
            self.sFacts[kF]=facts[kF]
            self.sUnits[kF]=units[kF]

    @staticmethod
    def FromFile(fileName,maxRows=None):
        new=Fort23Data()
        print("parsing file",fileName,"...")
        if ( not os.path.exists(fileName) ):
            print("...file does not exist!")
        else:
            # parse file
            if (maxRows is not None):
                dumpData=np.loadtxt(fileName,max_rows=maxRows)
            else:
                dumpData=np.loadtxt(fileName)
            print("...acquired",len(dumpData),"lines!")
            # assign data
            print("storing data...")
            new.IJ=dumpData[:,0]
            new.XX=dumpData[:,1];  new.YY=dumpData[:,2];  new.ZZ=dumpData[:,3]
            new.TXX=dumpData[:,4]; new.TYY=dumpData[:,5]; new.TZZ=dumpData[:,6]
            new.ETRACK=dumpData[:,7]; new.PTRACK=dumpData[:,8]
            if (dumpData.shape[1]>9):
                new.DEFL=dumpData[:,9 ]
                new.DPFL=dumpData[:,10]
                new.DTFL=dumpData[:,11]
        print("...done.")
        return new

    def RoundZ(self,Z0s,prec=None):
        print("rounding Z values...")
        uZs=np.unique(self.ZZ)
        print("...unique %d Z values BEFORE rounding:"%(len(uZs)))
        print(uZs)
        if (prec is None):
            npZ0s=np.array(Z0s)
            prec=np.min(npZ0s[npZ0s>0.0])*0.995
            # print(prec)
        for ii,Z0 in enumerate(Z0s):
            if (Z0!=0.0):
                indices=np.where(np.absolute(self.ZZ/Z0-1)<=prec)[0]
            else:
                indices=np.where(np.absolute(self.ZZ)<=prec)[0]
            self.ZZ[indices]=Z0
        uZs=np.unique(self.ZZ)
        print("...unique %d Z values AFTER rounding:"%(len(uZs)))
        print(uZs)
        print("...done.")

    def RemoveDieing(self,lDebug=True):
        print("removing particles dieing along their path...")
        uZs=np.unique(self.ZZ); nZs=len(uZs)
        print("...unique values of Z:",uZs)
        iKills=[]; ii=0
        while (ii<len(self)):
            tmpZZ=self.ZZ[ii:ii+nZs]
            if (np.array_equal(tmpZZ,uZs)):
                ii=ii+nZs
            else:
                print("...killing entry at line %d: ZZ=%g"%(ii+1,tmpZZ[0]))
                iKills.append(ii)
                ii=ii+1
        if (len(iKills)>0):
            print("...removing %d lines out of %d!"%(len(iKills),len(self)))
            self.IJ=np.delete(self.IJ,iKills)
            self.XX=np.delete(self.XX,iKills); self.YY=np.delete(self.YY,iKills); self.ZZ=np.delete(self.ZZ,iKills)
            self.TXX=np.delete(self.TXX,iKills); self.TYY=np.delete(self.TYY,iKills); self.TZZ=np.delete(self.TZZ,iKills)
            self.ETRACK=np.delete(self.ETRACK,iKills); self.PTRACK=np.delete(self.PTRACK,iKills)
            if (len(self.DEFL)>0): self.DEFL=np.delete(self.DEFL,iKills)
            if (len(self.DPFL)>0): self.DPFL=np.delete(self.DPFL,iKills)
            if (len(self.DTFL)>0): self.DTFL=np.delete(self.DTFL,iKills)
        else:
            print("...all %d lines survive!"%(len(self)))
        nRem=len(np.where(self.ZZ==uZs[0])[0])
        print("...remaining with %d particles;"%(nRem))
        lKill=False
        for ii,ZZ in enumerate(uZs):
            if (ii==0): continue
            nCurr=len(np.where(self.ZZ==ZZ)[0])
            if (nCurr!=nRem):
                print("...%d lines at Z=%g!"%(nCurr,ZZ))
                lKill=True
        if (lKill): exit(1)
        print("...done.")
                
    def GetAngles(self,lDebug=True):
        print("computing angles...")
        self.XP=self.TXX/self.TZZ
        self.YP=self.TYY/self.TZZ
        self.TH=np.sqrt(self.TXX**2+self.TYY**2)/self.TZZ
        if (lDebug): print("...min(TH):",np.min(self.TH))
        if (lDebug): print("...max(TH):",np.max(self.TH))
        print("...done.")

    def GetDeltaEnergy(self,pc0=None,Et0=None,Ek0=None,m0=None,lDebug=True):
        '''
        - Et0, Ek0 and m0 in MeV
        - pc0 in MeV/c
        '''
        print("computing delta energy and momentum...")
        if (pc0 is not None and Et0 is not None):
            print("...info supplied by user: pc0=%g MeV/c and Et0=%g MeV..."%(pc0,Et0))
        elif (m0 is not None and Ek0 is not None):
            print("...info supplied by user: using m0=%g MeV/c2 and Ek0=%g MeV..."%(m0,Ek0))
            Et0=Ek0+m0
            pc0=np.sqrt(Ek0*(Ek0+2*m0))
        else:
            print("...user should provide either pc0,Et0 or Ek0,m0")
            exit(1)
        print("...using pc0=%g MeV/c and Et0=%g MeV..."%(pc0,Et0))
        self.DE=self.ETRACK*1000-Et0
        self.DP=self.PTRACK*1000-pc0
        if (lDebug): print("...min(DE):",np.min(self.DE))
        if (lDebug): print("...max(DE):",np.max(self.DE))
        print("...done.")

    def GetDifferentialDataSet(self,iRef=None):
        print("differentiating populations...")
        if (iRef is not None): print("...iRef=%d..."%(iRef))
        new=deepcopy(self)
        # remove IJ, to avoid repetitions
        del new.IJ
        # remove transverse coordinates, to avoid misunderstandings
        del new.XX; del new.YY
        # remove direction cosine, to avoid misunderstandings
        del new.TXX; del new.TYY; del new.TZZ
        # remove absolute energies and momenta, to avoid misunderstandings
        del new.ETRACK; del new.PTRACK
        # only XP, YP, TH, DE, DP remain
        gc.collect()
        
        uZs=np.unique(new.ZZ)
        for ii,ZZ in reversed(list(enumerate(uZs))):
            indices=np.where(new.ZZ==ZZ)[0]
            print(indices)
            if (iRef is None):
                if (ii==0):
                    # first layer: just keep values
                    continue
                print(indices-1) # diff wrt previous layer
                new.XP[indices]=new.XP[indices]-new.XP[indices-1]
                new.YP[indices]=new.YP[indices]-new.YP[indices-1]
                # new.TH[indices]=np.abs(new.TH[indices]-new.TH[indices-1])
                new.DE[indices]=new.DE[indices]-new.DE[indices-1]
                new.DP[indices]=new.DP[indices]-new.DP[indices-1]
                # if (len(new.DEFL)>0): new.DEFL[indices]=new.DEFL[indices]-new.DEFL[indices-1]
                # if (len(new.DPFL)>0): new.DPFL[indices]=new.DPFL[indices]-new.DPFL[indices-1]
                # if (len(new.DTFL)>0): new.DTFL[indices]=np.absolute(new.DTFL[indices]-new.DTFL[indices-1])
            elif (iRef < 0):
                print(indices+iRef) # diff wrt n-th layer before
                new.XP[indices]=new.XP[indices]-new.XP[indices+iRef]
                new.YP[indices]=new.YP[indices]-new.YP[indices+iRef]
                # new.TH[indices]=np.abs(new.TH[indices]-new.TH[indices+iRef])
                new.DE[indices]=new.DE[indices]-new.DE[indices+iRef]
                new.DP[indices]=new.DP[indices]-new.DP[indices+iRef]
                # if (len(new.DEFL)>0): new.DEFL[indices]=new.DEFL[indices]-new.DEFL[indices+iRef]
                # if (len(new.DPFL)>0): new.DPFL[indices]=new.DPFL[indices]-new.DPFL[indices+iRef]
                # if (len(new.DTFL)>0): new.DTFL[indices]=np.absolute(new.DTFL[indices]-new.DTFL[indices+iRef])
            else:
                print("...shift=%d..."%(ii-iRef)) # diff wrt specific layer
                print(indices-(ii-iRef))
                new.XP[indices]=new.XP[indices]-new.XP[indices-(ii-iRef)]
                new.YP[indices]=new.YP[indices]-new.YP[indices-(ii-iRef)]
                # new.TH[indices]=np.abs(new.TH[indices]-new.TH[indices-(ii-iRef)])
                new.DE[indices]=new.DE[indices]-new.DE[indices-(ii-iRef)]
                new.DP[indices]=new.DP[indices]-new.DP[indices-(ii-iRef)]
                # if (len(new.DEFL)>0): new.DEFL[indices]=new.DEFL[indices]-new.DEFL[indices-(ii-iRef)]
                # if (len(new.DPFL)>0): new.DPFL[indices]=new.DPFL[indices]-new.DPFL[indices-(ii-iRef)]
                # if (len(new.DTFL)>0): new.DTFL[indices]=np.absolute(new.DTFL[indices]-new.DTFL[indices-(ii-iRef)])
        new.TH=np.sqrt(new.XP**2+new.YP**2)
        print("...done.")
        return new

    def LayerUnscattered(self,Zs=None):
        if (Zs is None): Zs=np.unique(self.ZZ)
        if (type(Zs) is not list): Zs=[Zs]
        for ii,ZZ in enumerate(Zs):
            indices=np.where(self.ZZ==ZZ)[0]
            nEntries=len(indices)
            print("layer at ZZ=%g - number of entries: %d"%(ZZ,nEntries))
#             print(len(self.TH[indices]==0.0),self.TH[indices]==0.0)
#             print(len(self.TH[indices[self.TH[indices]==0.0]]),self.TH[indices[self.TH[indices]==0.0]])
#             print(len(np.where(self.TH[indices]==0.0)[0]),np.where(self.TH[indices]==0.0)[0])
            noAngVar=len(np.where(self.TH[indices]==0.0)[0])
            print(" ...no angle variation: %d entries, i.e. %g%%"%(noAngVar,noAngVar/nEntries*100))
            noEn_Var=len(np.where(self.DE[indices]==0.0)[0])
            print(" ...no energy variation: %d entries, i.e. %g%%"%(noEn_Var,noEn_Var/nEntries*100))
            noAngEn_Var=len(np.where((self.TH[indices]==0.0) & (self.DE[indices]==0.0))[0])
#             print(np.where((self.TH[indices]==0.0) & (self.DE[indices]==0.0))[0])
            print(" ...no angle and energy variation: %d entries, i.e. %g%%"%(noAngEn_Var,noAngEn_Var/nEntries*100))
            if (len(self.DTFL)>0):
                noAngVar=len(np.where(self.DTFL[indices]==0.0)[0])
                print(" ...no angle variation (FLUKA inside): %d entries, i.e. %g%%"%(noAngVar,noAngVar/nEntries*100))
            if (len(self.DEFL)>0):
                noEn_Var=len(np.where(self.DEFL[indices]==0.0)[0])
                print(" ...no energy variation (FLUKA inside): %d entries, i.e. %g%%"%(noEn_Var,noEn_Var/nEntries*100))
            if (len(self.DTFL)>0 and len(self.DEFL)>0):
                noAngEn_Var=len(np.where((self.DTFL[indices]==0.0) & (self.DEFL[indices]==0.0))[0])
                print(" ...no angle and energy variation (FLUKA inside): %d entries, i.e. %g%%"%(noAngEn_Var,noAngEn_Var/nEntries*100))
                

    def LayerStats(self,myKeys=["XP","YP","TH"]):
        uZs=np.unique(self.ZZ)
        nZeros=np.zeros(len(uZs))
        for ii,ZZ in enumerate(uZs):
            indices=np.where(self.ZZ==ZZ)[0]
            nEntries=len(indices)
            print("layer at ZZ=%g - number of entries: %d"%(ZZ,nEntries))
            #
            for myKey in myKeys:
                mySet=self.__dict__[myKey][indices]
                nZeros[ii]=len(mySet[(mySet==0.0)])
                print("- %s: min=%g; max=%g; mean=%g; std=%g; 0-vals: %d (i.e. %g%%)"%(myKey,
                    np.min(mySet),np.max(mySet),np.mean(mySet),
                    np.std(mySet),nZeros[ii],nZeros[ii]/nEntries*100
                ))
        return nZeros/nEntries

    def GetMeanEnergy(self,myPercents=[2,98]):
        print("computing mean energy loss of (cut) population...")
        uZs=np.unique(self.ZZ)
        aves=np.empty((2,len(uZs))); aves[:]=np.nan
        for ii,ZZ in enumerate(uZs):
            indices=np.where(self.ZZ==ZZ)[0]
            aves[0,ii]=np.mean(CutPopulation(self.DE[indices],myPercents=myPercents))
            aves[1,ii]=np.mean(CutPopulation(self.DP[indices],myPercents=myPercents))
        print("...done.")
        return np.nan_to_num(-aves,posinf=0.0,neginf=0.0)

    def ShowScatter(self,Zs=None,LabZs=None,
                    myWhat=[["XX","XP"],["XX","YY"],["ETRACK","PTRACK"],["YY","YP"],["XP","YP"],["DE","DP"]],
                    myFacts=[[10,1000],[10,10],[1000,1000],[10,1000],[1000,1000],[1,1]],
                    myUnits=[["mm","mrad"],["mm","mm"],["MeV","MeV/c"],["mm","mrad"],["mrad","mrad"],["MeV","MeV/c"]],
                    figName=None,nCols=3):
        'just plotting'
        if (Zs is None): Zs=np.unique(self.ZZ)
        if (LabZs is None): LabZs=["Z=%g cm"% ZZ for ZZ in Zs ]
        nRows=round(len(myWhat)/nCols+0.001)
        fig, axs = plt.subplots(nRows,nCols,figsize=(5*nCols,5*nRows))
        for iZ in range(len(Zs)-1,-1,-1):
            # from most scattered to least scattered
            indices=np.where(self.ZZ==Zs[iZ])[0] ; kk=0
            for ii in range(nRows):
                for jj in range(nCols):
                    print("... %s vs %s for ZZ=%g cm..."%(myWhat[kk][1],myWhat[kk][0],Zs[iZ]))
                    if (nRows==1):
                        if (nCols==1):
                            myAx=axs
                        else:
                            myAx=axs[jj]
                    else:
                        myAx=axs[ii,jj]
                    myAx.plot(self.__dict__[myWhat[kk][0]][indices]*myFacts[kk][0],
                              self.__dict__[myWhat[kk][1]][indices]*myFacts[kk][1],
                              ".",label=LabZs[iZ])
                    if (iZ==0):
                        # last population
                        myAx.grid()
                        myAx.set(xlabel="%s [%s]"%(myWhat[kk][0],myUnits[kk][0]),
                                 ylabel="%s [%s]"%(myWhat[kk][1],myUnits[kk][1]))
                    kk=kk+1
                    if (kk==len(myWhat)): # last plot
                        myAx.legend()
                        if (kk%2==1): break # in case of odd plots
        plt.tight_layout()
        if ( figName is not None ):
            oFileName="%s_scatter.png"%(figName)
            print("saving plot to figure %s ..."%(oFileName))
            plt.savefig(oFileName)
        else:
            plt.show()
        plt.close(fig)

    def _HistMe(self,ZZ,LabZ,nBins=200,
                myWhat="TH",myFact=1E3,myUnit="mrad",
                myPc=[2,98], myCuts=None, lSkipZeros=True,
                figTitle=None,figName=None,picName=None,fitModels=None):
        'computing histograms and plotting'
        print("...hist of %s at ZZ=%g..."%(myWhat,ZZ))
        indices=np.where(self.ZZ==ZZ)[0]
        if (myWhat in self.sFacts.keys()):
            mySet=self.__dict__[myWhat][indices]*self.sFacts[myWhat]
            mUnit=self.sUnits[myWhat]
        else:
            mySet=self.__dict__[myWhat][indices]*myFact
            mUnit=myUnit
        origLen=len(mySet)
        print("   ...original data-set size: %d"%(origLen))
        # kill zeros
        if (lSkipZeros):
            print("   ...filtering out values exactly equal to 0...")
            currLen=len(mySet)
            mySet=mySet[(mySet!=0.0)]
            nRemoved=currLen-len(mySet)
            print("   ...removed %d entries, equivalent to %g%%!"%(nRemoved,nRemoved/origLen*100))
            print("   ...data-set size reduced to: %d"%(len(mySet)))
        # use percentiles
        if (myPc is not None):
            print("   ...percentiles: [%g:%g]"%(myPc[0],myPc[1]))
            currLen=len(mySet)
            mySet=CutPopulation(mySet,myPercents=myPc)
            nRemoved=currLen-len(mySet)
            print("   ...removed %d entries, equivalent to %g%%!"%(nRemoved,nRemoved/origLen*100))
            print("   ...data-set size reduced to: %d"%(len(mySet)))
        print("   ...range: [%g:%g]"%(min(mySet),max(mySet)))
        # compute histograms
        print("   ...actually histogramming...")
        hist,binEdges=np.histogram(mySet,bins=nBins)
        mids=binEdges[:-1]+0.5*(binEdges[1]-binEdges[0])
        # compute stats
        print("   ...computing basic stats...")
        myStats={
            "mean": np.nanmean(mySet),
            "median": np.nanmedian(mySet),
            "std": np.nanstd(mySet),
            "maxC": np.max(hist),
            "XofMax": mids[np.argmax(hist)]
        }

        # fits
        if (fitModels is None):
            fitParams=None; fitVals=None; fitMids=None
        else:
            if (type(fitModels) is not list): fitModels=[fitModels]
            if (myCuts is None): myCuts=[ None for ii in fitModels ]
            fitParams=[]; fitVals=[]; fitMids=[]
            for fitModel,myCut in zip(fitModels,myCuts):
                print("...fitting hist of %s at ZZ=%g with model %s..."%(myWhat,ZZ,fitModel))
                if (myCut is None):
                    currHist=hist; currMids=mids
                else:
                    jndices=np.where((myCut[0]<mids) & (mids<myCut[1]))[0]
                    currHist=hist[jndices]; currMids=mids[jndices]
                myFitPars,myFitVals=myFit(currHist,currMids,fitModel)
                fitParams.append(myFitPars); fitVals.append(myFitVals)
                fitMids.append(currMids)

        # plot
        histLabel="%s - %s"%(myWhat,LabZ)
        nCols=1; nRows=1
        fig, axs = plt.subplots(nRows,nCols,figsize=(10*nCols,5*nRows))
        axs = showHist(mids,hist,histLabel,fitModels,fitParams,fitMids,fitVals,mUnit,axs,myC="k")
        axs.grid(); axs.legend()
        axs.set(xlabel="%s [%s]"%(myWhat,mUnit),ylabel="counts")
        if (figTitle is not None):
            fig.suptitle(figTitle)
        plt.tight_layout()
        
        return hist, mids, fitParams, fitVals, fitMids, myStats

    def HistMe(self,Zs=None,LabZs=None,
               myPercents=[2,98],nBins=200,myDensity=False,
               myCuts=None, lSkipZeros=False,
               myWhat=[["XX","YY"],["XP","YP"],["DE","DP"]],
               myFacts=[[10,10],[1000,1000],[1,1]],
               myUnits=[["mm","mm"],["mrad","mrad"],["MeV","MeV/c"]],
               ZsZeros=[0.0],
               figTitle=None,
               figName=None,picName=None,fitModel=None):
        'computing histograms and plotting'
        if (Zs is None): Zs=np.unique(self.ZZ)
        if (LabZs is None): LabZs=["Z=%g cm"% ZZ for ZZ in Zs ]
        nCols=len(myWhat[0])
        nRows=len(myWhat)
        # computing histograms
        hist={}; mids={}
        for iZ in range(len(Zs)):
            hist[LabZs[iZ]]={} ; mids[LabZs[iZ]]={}
            indices=np.where(self.ZZ==Zs[iZ])[0]
            for ii in range(nRows):
                for jj in range(nCols):
                    print("...hist of %s at ZZ=%g..."%(myWhat[ii][jj],Zs[iZ]))
                    mySet=self.__dict__[myWhat[ii][jj]][indices]*myFacts[ii][jj]
                    origLen=len(mySet)
                    print("   ...original data-set size: %d"%(origLen))
                    # use percentiles
                    if (myPercents is not None):
                        myPc=[0.0,0.0]
                        if (type(myPercents[0]) is list):
                            myPc[0]=myPercents[ii][jj][0]; myPc[1]=myPercents[ii][jj][1]
                        else:
                            myPc[0]=myPercents[0]; myPc[1]=myPercents[1]
                        print("   ...considering only population in percentiles [%g:%g] ..."%(myPc[0],myPc[1]))
                        mySet=CutPopulation(mySet,myPercents=myPc)
                        print("   ...data-set size reduced to: %d"%(len(mySet)))
                    # use cuts
                    if (myCuts is not None):
                        myC=[min(mySet),max(mySet)]
                        if (type(myPercents[0]) is list):
                            myC[0]=myCuts[ii][jj][0]; myC[1]=myCuts[ii][jj][1]
                        else:
                            myC[0]=myCuts[0]; myC[1]=myCuts[1]
                        print("   ...keeping population in range [%g:%g] ..."%(myC[0],myC[1]))
                        currLen=len(mySet)
                        mySet=mySet[(myC[0]<=mySet) & (mySet<=myC[1])]
                        nRemoved=currLen-len(mySet)
                        print("   ...removed %d entries, equivalent to %g%%!"%(nRemoved,nRemoved/origLen*100))
                        print("   ...data-set size reduced to: %d"%(len(mySet)))
                    # kill zeros
                    if (lSkipZeros and Zs[iZ] not in ZsZeros ):
                        print("   ...filtering out values exactly equal to 0...")
                        currLen=len(mySet)
                        mySet=mySet[(mySet!=0.0)]
                        nRemoved=currLen-len(mySet)
                        print("   ...removed %d entries, equivalent to %g%%!"%(nRemoved,nRemoved/origLen*100))
                        print("   ...data-set size reduced to: %d"%(len(mySet)))
                    print("   ...range: [%g:%g]"%(min(mySet),max(mySet)))
                    hist[LabZs[iZ]][myWhat[ii][jj]],binEdges=np.histogram(mySet,bins=nBins,density=myDensity)
                    mids[LabZs[iZ]][myWhat[ii][jj]]=binEdges[:-1]+0.5*(binEdges[1]-binEdges[0])
        # fitting
        if (fitModel is None):
            fitParameters=None; fitVals=None
        else:
            fitParameters={myWhat[ii][jj]:[] for ii in range(nRows) for jj in range(nCols)}
            fitVals={}
            for iZ in range(len(Zs)):
                fitVals[LabZs[iZ]]={} 
                for ii in range(nRows):
                    for jj in range(nCols):
                        if (Zs[iZ] not in ZsZeros ):
                            if (isinstance(fitModel,str)):
                                myFitModel=fitModel
                            elif (isinstance(fitModel,list)):
                                myFitModel=fitModel[ii][jj]
                            print("...fitting hist of %s at ZZ=%g with model %s..."%(myWhat[ii][jj],Zs[iZ],myFitModel))
                            myFitPars,myFitVals=myFit(hist[LabZs[iZ]][myWhat[ii][jj]],mids[LabZs[iZ]][myWhat[ii][jj]],myFitModel)
                        else:
                            myFitPars=[]; myFitVals=[]
                        fitParameters[myWhat[ii][jj]].append(myFitPars)
                        fitVals[LabZs[iZ]][myWhat[ii][jj]]=myFitVals

                        
        # plotting 
        fig, axs = plt.subplots(nRows,nCols,figsize=(5*nCols,5*nRows))
        for iZ in range(len(Zs)-1,-1,-1):
            # from most scattered to least scattered
            for ii in range(nRows):
                for jj in range(nCols):
                    if (nRows==1):
                        if (nCols==1):
                            myAx=axs
                        else:
                            myAx=axs[jj]
                    else:
                        myAx=axs[ii,jj]
                    if (myDensity):
                        yLab=r"pdf [%s$^{-1}$]"%(myUnits[ii][jj])
                    else:
                        yLab="counts"
                    myAx.step(mids[LabZs[iZ]][myWhat[ii][jj]],
                              hist[LabZs[iZ]][myWhat[ii][jj]],
                              where='mid',label=LabZs[iZ])
                    if (fitModel is not None and Zs[iZ] not in ZsZeros ):
                        if (isinstance(fitModel,str)):
                            myFitModel=fitModel
                        elif (isinstance(fitModel,list)):
                            myFitModel=fitModel[ii][jj]
                        if (myFitModel.upper().startswith("GAUSS")):
                            myMu=fitParameters[myWhat[ii][jj]][iZ][2]
                            mySg=fitParameters[myWhat[ii][jj]][iZ][0]
                            myLab=r"GAUSS fit($\mu$=%g %s,$\sigma$=%g %s)"%(
                                      myMu,myUnits[ii][jj],
                                      mySg,myUnits[ii][jj])
                        elif (myFitModel.upper().startswith("LGAUSS")):
                            mySg=fitParameters[myWhat[ii][jj]][iZ][0]
                            myLab=r"lGAUSS fit($\sigma$=%g %s)"%(
                                      mySg,myUnits[ii][jj])
                        elif (myFitModel.upper().startswith("RUTH")):
                            myLab=r"RUTH fit"
                        elif (myFitModel.upper().startswith("LANGAU")):
                            myMu=fitParameters[myWhat[ii][jj]][iZ][0]
                            mySg=fitParameters[myWhat[ii][jj]][iZ][2]
                            myLab=r"LANGAU fit($\mu$=%g %s,$\sigma$=%g %s)"%(
                                      myMu,myUnits[ii][jj],
                                      mySg,myUnits[ii][jj])
                        myAx.plot(mids[LabZs[iZ]][myWhat[ii][jj]],
                                  fitVals[LabZs[iZ]][myWhat[ii][jj]],
                                  'k--',label=myLab)
                    if (iZ==0):
                        # last population
                        myAx.grid()
                        myAx.set(xlabel="%s [%s]"%(myWhat[ii][jj],myUnits[ii][jj]),
                                 ylabel=yLab)
                    myAx.legend()
        if (figTitle is not None):
            fig.suptitle(figTitle)
        plt.tight_layout()
        if ( figName is not None ):
            oFileName="%s_hist.png"%(figName)
            print("saving plot to figure %s ..."%(oFileName))
            plt.savefig(oFileName)
        else:
            plt.show()
        plt.close(fig)
        # pickle
        if (picName is not None):
            with open("%s.pickle"%(picName),'wb') as ff:
                print("...saving histograms to pickle file...")
                pickles={"hist":hist,"mids":mids}
                pickle.dump(pickles,ff)
        return hist, mids, fitParameters, fitVals

def CutPopulation(what,myPercents=[2,98]):
    percentiles=np.percentile(what,myPercents)
    return what[np.where((percentiles[0]<what) & (what<percentiles[1]))[0]]

def Gauss(x,sig=1,A=1,x0=0):
    return A/(np.sqrt(2*np.pi)*sig)*np.exp(-(x-x0)**2/(2*sig**2))

def lGauss(x,sig=1,A=1):
    return x/sig*A*np.exp(-x**2/(2*sig**2))

def RuthScat(x,A,b,c,d):
    return A/(((x-b)/c)**4)+d

def ExpDecay(x,A,lam,x0):
    return A*np.exp(-(x-x0)/lam)

def ExpDecaySimple(x,lam):
    return 100*np.exp(-x/lam)

def myFit(hist,mids,fitModel):
    if (fitModel.upper().startswith("GAUSS")):
        if (np.isnan(hist).all()):
            # avoid errors when fitting histograms of tracking through vacuum
            fitParameters=np.array([0,1,np.mean(mids)])
            fitVals=np.empty(len(hist)); fitVals[:]=np.nan
        else:
            A=np.max(hist)
            mean=np.mean(mids)
            sigma=np.sqrt( np.sum(hist*(mids-mean)**2) / np.sum(hist) )
            param_bounds=([0,0,-np.inf],[np.inf,np.inf,np.inf])
            print("   ...starting vals: sig=%g,A=%g,mean=%g"%(sigma,A,mean))
            fitParameters,covariance=curve_fit( Gauss, mids, hist, bounds=param_bounds, p0=(sigma,A,mean), maxfev=50000 )
            print("   ...fit vals: sig=%g,A=%g,mean=%g"%(fitParameters[0],fitParameters[1],fitParameters[2]))
            fitVals=Gauss(mids,*fitParameters)
    elif (fitModel.upper().startswith("LGAUSS")):
        if (np.isnan(hist).all()):
            # avoid errors when fitting histograms of tracking through vacuum
            fitParameters=np.array([0,1])
            fitVals=np.empty(len(hist)); fitVals[:]=np.nan
        else:
            A=np.max(hist)
            sigma=np.sqrt( np.sum(hist*mids**2) / np.sum(hist) )
            param_bounds=([0,0],[np.inf,np.inf])
            print("   ...starting vals: sig=%g,A=%g"%(sigma,A))
            fitParameters,covariance=curve_fit( lGauss, mids, hist, bounds=param_bounds, p0=(sigma,A), maxfev=50000 )
            print("   ...fit vals: sig=%g,A=%g"%(fitParameters[0],fitParameters[1]))
            fitVals=lGauss(mids,*fitParameters)
    elif (fitModel.upper().startswith("RUTH")):
        if (np.isnan(hist).all()):
            # avoid errors when fitting histograms of tracking through vacuum
            fitParameters=np.zeros(4)
            fitVals=np.empty(len(hist)); fitVals[:]=np.nan
        else:
            A=np.max(hist); b=0.0; c=1.0; d=0.0
            print("   ...starting vals: A=%g,b=%g,c=%g,d=%g"%(A,b,c,d))
            param_bounds=([0,-np.inf,-np.inf,-np.inf],[np.inf,np.inf,np.inf,np.inf])
            fitParameters,covariance=curve_fit( RuthScat, mids, hist, bounds=param_bounds, p0=(A,b,c,d), maxfev=50000 )
            print("   ...fit vals: A=%g,b=%g,c=%g,d=%g"%(fitParameters[0],fitParameters[1],fitParameters[2],fitParameters[3]))
            fitVals=RuthScat(mids,*fitParameters)
    elif ( fitModel.upper().startswith("LANGAU") ):
        myHist=deepcopy(hist); myMids=deepcopy(mids); lFlip=False
        if (min(myMids)<0.0 and max(myMids)<0.0):
            print("   ...flipping sign;")
            myHist=myHist[::-1]
            myMids=-myMids[::-1]
            lFlip=True
        if (np.isnan(myHist).all()):
            # avoid errors when fitting histograms of tracking through vacuum
            fitParameters=np.array([np.mean(myMids),0,0,1])
            fitVals=np.empty(len(myHist)); fitVals[:]=np.nan
        else:
            mpv = myMids[np.argmax(myHist)]
            sigma = np.sqrt( np.sum(myHist*(myMids-mpv)**2) / np.sum(myHist) )
            # eta = sigma
            eta = sigma/50
            A = np.max(myHist)
            param_bounds=([mpv*0.8,0,0,0],[mpv*1.2,10,4*sigma,1000*A])
            print("   ...starting vals: mpv=%g,eta=%g,sigma=%g,A=%g"%(mpv,eta,sigma,A))
            fitParameters,covariance=curve_fit( pylandau.langau, myMids, myHist, bounds=param_bounds, p0=(mpv, eta, sigma, A ), maxfev=50000 )
            print("   ...fit vals: mpv=%g,eta=%g,sigma=%g,A=%g"%(fitParameters[0],fitParameters[1],fitParameters[2],fitParameters[3]))
            fitVals=pylandau.langau(myMids, *fitParameters)
            if (lFlip): fitVals=fitVals[::-1]
    elif ( fitModel.upper().startswith("LANDAU") ):
        myHist=deepcopy(hist); myMids=deepcopy(mids); lFlip=False
        if (min(myMids)<0.0 and max(myMids)<0.0):
            print("   ...flipping sign;")
            myHist=myHist[::-1]
            myMids=-myMids[::-1]
            lFlip=True
        if (np.isnan(myHist).all()):
            # avoid errors when fitting histograms of tracking through vacuum
            fitParameters=np.array([np.mean(myMids),0,0,1])
            fitVals=np.empty(len(myHist)); fitVals[:]=np.nan
        else:
            mpv = myMids[np.argmax(myHist)]
            eta = 1
            A = np.max(myHist)
            param_bounds=([mpv*0.8,0,0],[mpv*1.2,100,1000*A])
            print("   ...starting vals: mpv=%g,eta=%g,A=%g"%(mpv,eta,A))
            fitParameters,covariance=curve_fit( pylandau.landau, myMids, myHist, bounds=param_bounds, p0=(mpv, eta, A ), maxfev=50000 )
            print("   ...fit vals: mpv=%g,eta=%g,A=%g"%(fitParameters[0],fitParameters[1],fitParameters[2]))
            fitVals=pylandau.landau(myMids, *fitParameters)
            if (lFlip): fitVals=fitVals[::-1]
    return fitParameters, fitVals

def LynchDahl(x,betaRel,pc,Z,X0,lLog=True):
    Ys=13.6/(betaRel*pc)*Z*np.sqrt(x/X0)
    if (lLog):
        Ys=Ys*(1+0.038*np.log(x/X0*Z**2/betaRel**2))
    return Ys

def LynchDahlBetter(x,betaRel,pc,z,Z,A,rho,alpha=1/137.,FF=0.98):
    print("better Lynch and Dahl")
    ChiA2=2.007E-5*Z**(2/3.)/pc**2*(1+3.34*(Z*z*alpha/betaRel)**2) # screening angle
    ChiC2=0.157*(Z*(Z+1)*x*rho/A)*(z/(pc*betaRel))**2 # characteristic angles
    Omega=ChiC2/ChiA2
    vv=0.5*Omega/(1-FF)
    Ys=ChiC2/(1+FF**2)*((1+vv)/vv*np.log(1+vv)-1)
    print("ChiA:",np.sqrt(ChiA2))
    print("ChiC:",np.sqrt(ChiC2))
    print("Omega:",Omega)
    print("vv:",vv)
    print("Ys:",Ys)
    return Ys

def AnalyseSetAndCore(mySet,nBinsAll=2000,nBinsCore=200,myPercents=[2,98],myDensity=True,myModel="Gauss"):
    # histogram on entire popoulation
    hist,binEdges=np.histogram(mySet,bins=nBinsAll,density=myDensity)
    mids=binEdges[:-1]+0.5*(binEdges[1]-binEdges[0])
    # histogram and fit on core popoulation
    coreHist, coreMids, parameters, fit_vals = AnalyseCore(mySet,nBins=nBinsCore,myPercents=myPercents,myDensity=myDensity,myModel=myModel)
    if ( myModel=="Langau" ):
        mids=-mids[::-1]; hist=hist[::-1]
    return hist, mids, coreHist, coreMids, parameters, fit_vals
def AnalyseCore(mySet,nBins=200,myPercents=[2,98],myDensity=True,myModel=None):
    # histogram
    hist,binEdges=np.histogram(CutPopulation(mySet,myPercents=myPercents),bins=nBins,density=myDensity)
    mids=binEdges[:-1]+0.5*(binEdges[1]-binEdges[0])
    # fit
    if ( myModel is None ):
        parameters=np.array([0,1,np.mean(mids)])
        fit_vals=np.empty(len(hist)); fit_vals[:]=np.nan
    elif ( myModel=="Gauss" ):
        if (np.isnan(hist).all()):
            # avoid errors when fitting histograms of tracking through vacuum
            parameters=np.array([0,1,np.mean(mids)])
            fit_vals=np.empty(len(hist)); fit_vals[:]=np.nan
        else:
            param_bounds=([0,0,-np.inf],[np.inf,np.inf,np.inf])
            parameters, covariance = curve_fit( Gauss, mids, hist, bounds=param_bounds, maxfev=5000 )
            fit_vals=Gauss(mids,parameters[0],parameters[1],parameters[2])
    elif ( myModel=="Langau" ):
        mids=-mids[::-1]; hist=hist[::-1]
        if (np.isnan(hist).all()):
            # avoid errors when fitting histograms of tracking through vacuum
            parameters=np.array([np.mean(mids),0,0,1])
            fit_vals=np.empty(len(hist)); fit_vals[:]=np.nan
        else:
            # param_bounds=([0,0,-np.inf],[np.inf,np.inf,np.inf])
            # parameters, covariance = curve_fit( Gauss, mids, hist, bounds=param_bounds, maxfev=5000 )
            # fit_vals=Gauss(mids,parameters[0],parameters[1],parameters[2])
            myMaxY=np.max(hist);
            # myInd=np.argmax(hist); myMaxX=mids[myInd]
            # myInd=np.argmin(hist); myMinX=mids[myInd]
            # mpv, eta, sigma, A = myMaxX, 0.1, 1E-4, myMaxY  # protons
            # mpv, eta, sigma, A = (myMaxX+myMinX)/2., 0.1, 1E-2, myMaxY  # carbon
            # mpv = np.average( mids, weights=hist )
            mpv = mids[np.argmax(hist)]
            eta = 0.1
            sigma = np.sqrt( np.sum(hist*(mids-mpv)**2) / np.sum(hist) )
            A = myMaxY
            param_bounds=([0,0,0,0],[np.inf,np.inf,np.inf,np.inf])
            parameters, covariance = curve_fit( pylandau.langau, mids, hist, bounds=param_bounds, p0=(mpv, eta, sigma, A ), maxfev=5000 )
            fit_vals=pylandau.langau(mids, *parameters)
            print(*parameters)
    else:
        print("unknown fitting model:",myModel)
        return
    return hist, mids, parameters, fit_vals

def PlotSetAndCore(hist, mids, coreHist, coreMids, parameters, fit_vals, ax, axLabel="xp", axUnit="rad", myTitle="", nLevel=3,myModel="Gauss"):
    ax.step(mids,hist,where='mid',label="all points")
    PlotCore(coreHist, coreMids, parameters, fit_vals, ax, axLabel=axLabel, axUnit=axUnit, myTitle=myTitle, nLevel=nLevel,myModel=myModel)
def PlotCore(hist, mids, parameters, fit_vals, ax, axLabel="xp", axUnit="rad", myTitle="", nLevel=3, myModel="Gauss"):
    ax.step(mids,hist,where='mid', label="beam core")
    ax.plot(mids,fit_vals,"r-",label="Fit")
    ax.legend()
    ax.grid()
    if ( myModel=="Gauss" ):
        sig=parameters[0]
        mean=parameters[2]
    elif ( myModel=="Langau" ):
        sig=parameters[2]
        mean=parameters[0]
    else:
        print("unknown fitting model:",myModel)
        return
    ax.set_title("{case} - $\sigma$={sig:.3E} {unit} - $\mu$={mean:.3E} {unit}".format(sig=sig,unit=axUnit,mean=mean,case=myTitle))
    ax.set(xlabel="{xWhat} [{xUnit}]".format(xWhat=axLabel,xUnit=axUnit), ylabel=r"pdf [({yUnit})$^{{-1}}$]".format(yUnit=axUnit))
    ax.set_xlim([-nLevel*sig+mean,nLevel*sig+mean])

def Acquire(fileNameAll,fileNameOnlyDDS,Ek0):
    # full data set
    DataSetAll=Fort23Data.FromFile(fileNameAll)
    DataSetAll.GetDeltaEnergy(Ek0)
    # only DDS data set
    DataSetOnlyDDS=Fort23Data.FromFile(fileNameOnlyDDS)
    DataSetOnlyDDS.GetDeltaEnergy(Ek0)
    return DataSetAll, DataSetOnlyDDS

def ExtendMatrix(MyMat,AddMat):
    NewMat=np.zeros((MyMat.shape[0]+AddMat.shape[0],max([MyMat.shape[1],AddMat.shape[1]])))
    NewMat[0:MyMat.shape[0],0:MyMat.shape[1]]=MyMat
    NewMat[MyMat.shape[0]:MyMat.shape[0]+AddMat.shape[0],0:AddMat.shape[1]]=AddMat
    return NewMat

def AnalyzeAngles(DataSetAll,DataSetOnlyDDS):
    print("analysing angle sets...")
    DataSetAll.GetAngles()
    DataSetOnlyDDS.GetAngles()
    sigmas=np.empty((2,4)); sigmas[:]=np.nan
    # analyse and plot
    fig, axs = plt.subplots(3,2,figsize=(15,10))
    # - ony vacuum window
    indices=np.where(DataSetAll.ZZ==-96)[0]
    hist, mids, coreHist, coreMids, parameters, fit_vals = AnalyseSetAndCore(DataSetAll.XP[indices])
    PlotSetAndCore(hist, mids, coreHist, coreMids, parameters, fit_vals, axs[0,0], axLabel="xp", axUnit="rad", myTitle="XP after Vacuum Window")
    sigmas[0,0]=parameters[0]
    hist, mids, coreHist, coreMids, parameters, fit_vals = AnalyseSetAndCore(DataSetAll.YP[indices])
    PlotSetAndCore(hist, mids, coreHist, coreMids, parameters, fit_vals, axs[0,1], axLabel="yp", axUnit="rad", myTitle="YP after Vacuum Window")
    sigmas[1,0]=parameters[0]
    # - ony DDS
    hist, mids, coreHist, coreMids, parameters, fit_vals = AnalyseSetAndCore(DataSetOnlyDDS.XP)
    PlotSetAndCore(hist, mids, coreHist, coreMids, parameters, fit_vals, axs[1,0], axLabel="xp", axUnit="rad", myTitle="XP after DDS")
    sigmas[0,1]=parameters[0]
    hist, mids, coreHist, coreMids, parameters, fit_vals = AnalyseSetAndCore(DataSetOnlyDDS.YP)
    PlotSetAndCore(hist, mids, coreHist, coreMids, parameters, fit_vals, axs[1,1], axLabel="yp", axUnit="rad", myTitle="YP after DDS")
    sigmas[1,1]=parameters[0]
    # - both
    indices=np.where(DataSetAll.ZZ==-60)[0]
    hist, mids, coreHist, coreMids, parameters, fit_vals = AnalyseSetAndCore(DataSetAll.XP[indices])
    PlotSetAndCore(hist, mids, coreHist, coreMids, parameters, fit_vals, axs[2,0], axLabel="xp", axUnit="rad", myTitle="XP after all")
    sigmas[0,2]=parameters[0]
    hist, mids, coreHist, coreMids, parameters, fit_vals = AnalyseSetAndCore(DataSetAll.YP[indices])
    PlotSetAndCore(hist, mids, coreHist, coreMids, parameters, fit_vals, axs[2,1], axLabel="yp", axUnit="rad", myTitle="YP after all")
    sigmas[1,2]=parameters[0]
    # - energy 
    # hist, mids, coreHist, coreMids, parameters, fit_vals = AnalyseSetAndCore(DataSetOnlyDDS.DE)
    # PlotSetAndCore(hist, mids, coreHist, coreMids, parameters, fit_vals, axs[1,0], axLabel=r"$\DeltaE$", axUnit="MeV")
    # hist, mids, coreHist, coreMids, parameters, fit_vals = AnalyseSetAndCore(DataSetOnlyDDS.DP)
    # PlotSetAndCore(hist, mids, coreHist, coreMids, parameters, fit_vals, axs[1,1], axLabel=r"$\Deltap$", axUnit="MeV/c")
    plt.tight_layout()
    plt.show()
    sigmas[:,3]=np.sqrt(sigmas[:,0]**2+sigmas[:,1]**2)
    print("...done.")
    return sigmas

def WriteToXLSX(wb,Ek0,pc0,RR,myCase,what,sheetNames=["xp","yp"],headers=["VW [rad]", "DDS [rad]", "both [rad]", "combined [rad]"],iRow=1):
    if (iRow==1):
        # first row of data: new file, write header
        WriteToXLSXheaders(wb,sheetNames=sheetNames,headers=headers)
    WriteToXLSXvals(wb,Ek0,pc0,RR,myCase,what,sheetNames=sheetNames,headers=headers,iRow=iRow)
def WriteToXLSXvals(wb,Ek0,pc0,RR,myCase,what,sheetNames=["xp","yp"],headers=["VW [rad]", "DDS [rad]", "both [rad]", "combined [rad]"],iRow=1):
    for iSheet,sheetName in enumerate(sheetNames):
        ws=wb.get_worksheet_by_name(sheetName)
        ws.write_row(iRow,0,[RR,Ek0,pc0,myCase])
        ws.write_row(iRow,4,what[iSheet,:])
def WriteToXLSXheaders(wb,sheetNames=["xp","yp"],headers=["VW [rad]", "DDS [rad]", "both [rad]", "combined [rad]"]):
    for iSheet,sheetName in enumerate(sheetNames):
        ws=wb.add_worksheet(sheetName)
        ws.write_row(0,0,["R [mm]","Ek [MeV]","Dp [MeV/c]","Case"])
        ws.write_row(0,4,headers)
    
def MainOnlyAngles():
    # vacuum
    fileNameAll="p320mm/MCS_VacWindow_DDS_exp_p320_All001_fort.23.gz"
    fileNameOnlyDDS="p320mm/MCS_VacWindow_DDS_exp_p320_onlyDDS001_fort.23.gz"
    # air
    fileNameAll="p320mm/MCS_VacWindow_DDS_exp_p320_Air_All001_fort.23.gz"
    fileNameOnlyDDS="p320mm/MCS_VacWindow_DDS_exp_p320_Air_onlyDDS001_fort.23.gz"
    xlsxFileName="test.xlsx"
    Ek0=228.57 # p@320mm [MeV]
    DataSetAll, DataSetOnlyDDS = Acquire(fileNameAll,fileNameOnlyDDS,Ek0)
    sigmas=AnalyzeAngles(DataSetAll,DataSetOnlyDDS)
    WriteToXLSX(xlsxFileName,Ek0,sigmas)

def AnalyseGaussianEnergies(DataSet,Zs,labels,particle="PROTON",figName=None,picName=None,myModel=None):
    print("analysing energies...")
    if (myModel is None):
        print("...no fitting...")
    else:
        print("...with %s as fitting model..."%(myModel))
    # analyze data
    iSigs=np.arange(len(Zs))
    peaks=np.empty((2,len(Zs))); peaks[:]=np.nan
    fitParams=np.empty((2,4*len(Zs))); peaks[:]=np.nan
    if (picName is not None): pickles={"DE":{},"DP":{}}
    # analyse and plot: position
    figEne, axs = plt.subplots(len(Zs),2,figsize=(5*len(Zs),10))
    for thick,myTitle,iSig in zip(Zs,labels,iSigs):
        print(thick,myTitle,iSig)
        indices=np.where(DataSet.ZZ==thick)[0]
        #   DE
        if (particle=="PROTON"):
            hist, mids, coreHist, coreMids, parameters, fit_vals = AnalyseSetAndCore(DataSet.DE[indices]*10,myModel=myModel)
            PlotSetAndCore(hist, mids, coreHist, coreMids, parameters, fit_vals, axs[iSig,0], axLabel="$\Delta$E", axUnit="100 keV", myTitle=myTitle, nLevel=30, myModel=myModel)
        elif (particle=="CARBON"):
            hist, mids, coreHist, coreMids, parameters, fit_vals = AnalyseSetAndCore(DataSet.DE[indices],myModel=myModel)
            PlotSetAndCore(hist, mids, coreHist, coreMids, parameters, fit_vals, axs[iSig,0], axLabel="$\Delta$E", axUnit="MeV", myTitle=myTitle, nLevel=5, myModel=myModel)
        peaks[0,iSig]=parameters[0]
        fitParams[0,iSig*4:(iSig+1)*4]=parameters
        if (picName is not None):
            pickles["DE"][myTitle]={"hist":coreHist,"mids":coreMids,"fit_vals":fit_vals,"parameters":parameters}
        #   DP
        if (particle=="PROTON"):
            hist, mids, coreHist, coreMids, parameters, fit_vals = AnalyseSetAndCore(DataSet.DP[indices]*10,myModel=myModel)
            PlotSetAndCore(hist, mids, coreHist, coreMids, parameters, fit_vals, axs[iSig,1], axLabel="$\Delta$p", axUnit="100 keV/c", myTitle=myTitle, nLevel=30, myModel=myModel)
        elif (particle=="CARBON"):
            hist, mids, coreHist, coreMids, parameters, fit_vals = AnalyseSetAndCore(DataSet.DP[indices],myModel=myModel)
            PlotSetAndCore(hist, mids, coreHist, coreMids, parameters, fit_vals, axs[iSig,1], axLabel="$\Delta$p", axUnit="MeV/c", myTitle=myTitle, nLevel=5, myModel=myModel)
        peaks[1,iSig]=parameters[0]
        fitParams[1,iSig*4:(iSig+1)*4]=parameters
        if (picName is not None):
            pickles["DP"][myTitle]={"hist":coreHist,"mids":coreMids,"fit_vals":fit_vals,"parameters":parameters}
    plt.tight_layout()
    # show results
    if ( figName is not None ):
        oFileName="%s_ENEs.png"%(figName)
        print("saving plot to figure %s..."%(oFileName))
        plt.savefig(oFileName)
    else:
        plt.show()
    plt.close(figEne)
    # save pickle
    if (picName is not None):
        with open("%s_ENEs.pickle"%(picName),'wb') as ff:
            print("saving pickle file...")
            pickle.dump(pickles,ff)
    print("...done.")
    return peaks, fitParams
    
def AnalyseGaussianPositions(DataSet,Zs,labels,figName=None,picName=None,myPercents=[2,98],myModel=None):
    print("analysing positions...")
    if (myModel is None):
        print("...no fitting...")
    else:
        print("...with %s as fitting model..."%(myModel))
    # analyze data
    iSigs=np.arange(len(Zs))
    sigmas=np.empty((2,len(Zs))); sigmas[:]=np.nan
    if (picName is not None): pickles={"XX":{},"YY":{}}
    # analyse and plot: position
    figPos, axsPos = plt.subplots(len(Zs),2,figsize=(5*len(Zs),10))
    for thick,myTitle,iSig in zip(Zs,labels,iSigs):
        print(thick,myTitle,iSig)
        indices=np.where(DataSet.ZZ==thick)[0]
        #   x
        hist, mids, parameters, fit_vals = AnalyseCore(DataSet.XX[indices],myPercents=myPercents,myModel=myModel)
        PlotCore(hist, mids, parameters, fit_vals, axsPos[iSig,0], axLabel="x", axUnit="cm", myTitle="X "+myTitle)
        sigmas[0,iSig]=parameters[0]
        if (picName is not None):
            pickles["XX"][myTitle]={"hist":hist,"mids":mids,"fit_vals":fit_vals,"parameters":parameters}
        #   y
        hist, mids, parameters, fit_vals = AnalyseCore(DataSet.YY[indices],myPercents=myPercents)
        PlotCore(hist, mids, parameters, fit_vals, axsPos[iSig,1], axLabel="y", axUnit="cm", myTitle="Y "+myTitle)
        sigmas[1,iSig]=parameters[0]
        if (picName is not None):
            pickles["YY"][myTitle]={"hist":hist,"mids":mids,"fit_vals":fit_vals,"parameters":parameters}
    plt.tight_layout()
    # show results
    if ( figName is not None ):
        oFileName="%s_POSs.png"%(figName)
        print("saving plot to figure %s..."%(oFileName))
        plt.savefig(oFileName)
    else:
        plt.show()
    plt.close(figPos)
    # save pickle
    if (picName is not None):
        with open("%s_POSs.pickle"%(picName),'wb') as ff:
            print("saving pickle file...")
            pickle.dump(pickles,ff)
    # 
    print("...done.")
    return sigmas

def AnalyseGaussianAngles(DataSet,Zs,labels,figName=None,picName=None,myPercents=[2,98],myModel=None):
    print("analysing angles...")
    if (myModel is None):
        print("...no fitting...")
    else:
        print("...with %s as fitting model..."%(myModel))
    # analyze data
    iSigs=np.arange(len(Zs))
    sigmas=np.empty((2,len(Zs))); sigmas[:]=np.nan
    if (picName is not None): pickles={"XP":{},"YP":{}}
    # analyse and plot: angles
    figAng, axsAng = plt.subplots(len(Zs),2,figsize=(5*len(Zs),10))
    for thick,myTitle,iSig in zip(Zs,labels,iSigs):
        print(thick,myTitle,iSig)
        indices=np.where(DataSet.ZZ==thick)[0]
        #   xp
        hist, mids, parameters, fit_vals = AnalyseCore(DataSet.XP[indices],myPercents=myPercents,myModel=myModel)
        PlotCore(hist, mids, parameters, fit_vals, axsAng[iSig,0], axLabel="xp", axUnit="rad", myTitle="XP "+myTitle)
        sigmas[0,iSig]=parameters[0]
        if (picName is not None):
            pickles["XP"][myTitle]={"hist":hist,"mids":mids,"fit_vals":fit_vals,"parameters":parameters}
        #   yp
        hist, mids, parameters, fit_vals = AnalyseCore(DataSet.YP[indices],myPercents=myPercents,myModel=myModel)
        PlotCore(hist, mids, parameters, fit_vals, axsAng[iSig,1], axLabel="yp", axUnit="rad", myTitle="YP "+myTitle)
        sigmas[1,iSig]=parameters[0]
        if (picName is not None):
            pickles["YP"][myTitle]={"hist":hist,"mids":mids,"fit_vals":fit_vals,"parameters":parameters}
    plt.tight_layout()
    # show results
    if ( figName is not None ):
        oFileName="%s_ANGs.png"%(figName)
        print("saving plot to figure %s..."%(oFileName))
        plt.savefig(oFileName)
    else:
        plt.show()
    plt.close(figAng)
    # save pickle
    if (picName is not None):
        with open("%s_ANGs.pickle"%(picName),'wb') as ff:
            print("saving pickle file...")
            pickle.dump(pickles,ff)
    print("...done.")
    return sigmas

def ReturnMassEnergy(myPart,myRange):
    if ( myPart.startswith("P") ):
        MPdata=np.loadtxt("./R1",skiprows=1)
        Ek=MPdata[np.where(MPdata[:,3]==myRange)[0],2]
        m0=938.2720813 # p (FLUKA value) [MeV/c2]
    elif ( myPart.startswith("C") ):
        MPdata=np.loadtxt("./R1_RIFI",skiprows=1)
        Ek=12*MPdata[np.where(MPdata[:,3]==myRange)[0],1]
        m0=11.1748642E3 # C (FLUKA value) [MeV/c2]
    else:
        print("...unknown particle %s!"%(myPart))
        exit(1)
    return m0, Ek

def LoadPickles(pickles=None,
                particle="PROTON",   # CARBON, PROTON
                ranges=[30],         # [mm]
                media=["Air","Vac"], # "Air","Vac"
                RiFis=[False] ):     # RiFis are in (True, only CARBON) or not (False, CARBON and PROTON)
    if (pickles is None): pickles={}
    if (particle not in pickles.keys()): pickles[particle]={}
    # syntax: pickles["PROTON"/"CARBON"]["030"]["NOrifi"/"WITHrifi"]["Vac"/"AIR"]["FULL"/"DIFF"]["XX"/"YY"/"XP"/"YP"/"DE"/"DP"]["hist"/"mids"/"fit_vals"/"parameters"]
    for myRange in ranges:
        sRange="%03i"%(myRange)
        if (sRange not in pickles[particle].keys()): pickles[particle][sRange]={}
        for RiFi in RiFis:
            if (RiFi):
                sRiFi="WITHrifi"
            else:
                sRiFi="NOrifi"
            if (sRiFi not in pickles[particle][sRange].keys()): pickles[particle][sRange][sRiFi]={}
            for medium in media:
                fort23FileName, xlsxFileName, figName, picName, figNameDiff, picNameDiff = InflatePaths(particle,myRange,medium,RiFi=RiFi)
                # check path to pickle files
                myPath=os.path.dirname(picName)
                if not os.path.isdir(myPath):
                    print("...path to pickle folder %s does NOT exist - skipping..."%(myPath))
                    continue
                if (medium not in pickles[particle][sRange][sRiFi].keys()): pickles[particle][sRange][sRiFi][medium]={}
                for mySet1 in ["POSs","ANGs","ENEs"]:
                    for mySet2 in ["","diff_"]:
                        myFile="%s_%s%s.pickle"%(picName,mySet2,mySet1)
                        if not os.path.isfile(myFile):
                            print("...path to pickle file %s does NOT exist - skipping..."%(myFile))
                            continue
                        print("...loading data from pickle file %s ..."%(myFile))
                        sSet2="FULL"
                        if (mySet2=="diff_"): sSet2="DIFF"
                        if (sSet2 not in pickles[particle][sRange][sRiFi][medium].keys()): pickles[particle][sRange][sRiFi][medium][sSet2]={}
                        with open(myFile,"rb") as ff:
                            tmp=pickle.load(ff)
                        for key,value in tmp.items():
                            pickles[particle][sRange][sRiFi][medium][sSet2][key]=value
    return pickles

def singleLayerAnalysis(myDataSet,ZZ,LabZZ,figName=None,picName=None,myTitle=None,whats=None,cuts=None,fitModels=None,nBins=None):

    hist={}; mids={}; fitParams={}; fitVals={}; fitMids={}; myStats={}
    
    myFigTitle=myTitle
    if (myFigTitle is not None): myFigTitle="%s - %s"%(myTitle,LabZZ)
    myFigName=figName
    if (myFigName is not None): myFigName="%s"%(figName)
    myPicName=picName
    if (myPicName is not None): myPicName="%s"%(picName)
    myDataSet.LayerUnscattered(Zs=ZZ)
    for myWhat in whats:
        if (myWhat in ["TH","XP","YP","DTFL"]):
            myPc=[0,99]
        elif (myWhat in ["XX","YY"]):
            myPc=[1,99]
        elif (myWhat in ["DE","DEFL"]):
            myPc=[2,100]
        hist[myWhat], mids[myWhat], fitParams[myWhat], fitVals[myWhat], fitMids[myWhat], myStats[myWhat] = \
           myDataSet._HistMe(ZZ,LabZZ,nBins=nBins[myWhat],myPc=myPc,myWhat=myWhat, 
                       myCuts=cuts[myWhat], fitModels=fitModels[myWhat],
                       figTitle=myFigTitle,figName=myFigName,picName=myPicName)
        
    return hist, mids, fitParams, fitVals, fitMids, myStats

def UnscatterAnalysis(DataSet,myKey=["TH"],myTitle=None):
    Zs=np.unique(DataSet.ZZ)
    # get differential data set (consecutive layers)
    DiffDataSet=DataSet.GetDifferentialDataSet()
    # 0-interactors (consecutive layers)
    nSurv=DiffDataSet.LayerStats(myKey)
    myXs=np.diff(Zs); myYs=nSurv[1:]
    # 0-interactors (cumulative distance)
    DiffDataSet_start=DataSet.GetDifferentialDataSet(iRef=0)
    nSurv_start=DiffDataSet_start.LayerStats(myKey)
    myXs=np.append(myXs,np.array(Zs))*1E4
    myYs=np.append(myYs,np.array(nSurv_start))*100
    indices=np.argsort(myXs)
    myXs=myXs[indices]; myYs=myYs[indices]
    u,indices=np.unique(myXs,return_index=True)
    myXs=myXs[indices]; myYs=myYs[indices]
    # - fit
    A=100; x0=0
    lam=0.15 # GOLD
    # lam=10 # KAPTON
    # lam=1 # Light GOLD
    indices=np.where(myYs>0.0)[0]
    print(indices,myXs[indices],myYs[indices])
    if (len(indices)>5):
        fitParameters,covariance=curve_fit( ExpDecay, myXs[indices], myYs[indices], p0=(A,lam,x0), maxfev=50000 )
        fitVals=ExpDecay(myXs,*fitParameters)
        lambdaFit=fitParameters[1]
    else:
        fitParameters,covariance=curve_fit( ExpDecaySimple, myXs[indices], myYs[indices], p0=(lam), maxfev=50000 )
        fitVals=ExpDecaySimple(myXs,*fitParameters)
        lambdaFit=fitParameters[0]
    
    # plot
    fig, myAx = plt.subplots()
    myAx.plot(myXs,myYs,"bo",label="no interactions (all)")
    myAx.plot(np.array(Zs)*1E4,np.array(nSurv_start)*100,"r.",label="no interactions (wrt entrance)")
    myAx.plot(myXs,fitVals,"g-",label=r"fit: $\lambda$=%g $\mu$m"%(lambdaFit))
    myAx.grid()
    myAx.set(ylabel="[%]",xlabel=r"thickness [$\mu$m]")
    myAx.legend()
    plt.yscale("log")
    fig.suptitle(myTitle)
    plt.tight_layout()
    plt.show()
    print("fit parameters:",fitParameters)
    print("covariance matrix:",covariance)
    
def StandardAnalysis(fort23FileName,Zs,LabZs,BeamPart,X0,figName=None,picName=None,myTitle=None):
    # get data set
    DataSet=Fort23Data.FromFile(fort23FileName)
    # DataSet.RoundZ(Z0s=Zs,prec=2E-5) # KAPTON and LGOLD, SS
    # DataSet.RoundZ(Zs,prec=1.2E-6) # GOLD, SS
    DataSet.RoundZ(Zs)
    # exit()
    DataSet.RemoveDieing()
    # get angles
    DataSet.GetAngles()
    # get energy losses
    DataSet.GetDeltaEnergy(BeamPart.Ek,BeamPart.mm)
    # get differential data set (consecutive layers)
    DiffDataSet=DataSet.GetDifferentialDataSet()
    # 0-interactors (consecutive layers)
    nSurv=DiffDataSet.LayerStats()
    myXs=np.diff(Zs); myYs=nSurv[1:]
    # 0-interactors (cumulative distance)
    DiffDataSet_start=DataSet.GetDifferentialDataSet(iRef=0)
    nSurv_start=DiffDataSet_start.LayerStats()
    myXs=np.append(myXs,np.array(Zs))*1E4
    myYs=np.append(myYs,np.array(nSurv_start))*100
    indices=np.argsort(myXs)
    myXs=myXs[indices]; myYs=myYs[indices]
    u,indices=np.unique(myXs,return_index=True)
    myXs=myXs[indices]; myYs=myYs[indices]
    # - fit
    A=100; x0=0
    lam=0.1 # GOLD
    # lam=10 # KAPTON
    # lam=1 # Light GOLD
    indices=np.where(myYs>0.0)[0]
    print(indices,myXs[indices],myYs[indices])
    if (len(indices)>3):
        fitParameters,covariance=curve_fit( ExpDecay, myXs[indices], myYs[indices], p0=(A,lam,x0), maxfev=50000 )
        fitVals=ExpDecay(myXs,*fitParameters)
        lambdaFit=fitParameters[1]
    else:
        fitParameters,covariance=curve_fit( ExpDecaySimple, myXs[indices], myYs[indices], p0=(lam), maxfev=50000 )
        fitVals=ExpDecaySimple(myXs,*fitParameters)
        lambdaFit=fitParameters[0]
    
    # plot
    fig, myAx = plt.subplots()
    myAx.plot(myXs,myYs,"bo",label="no interactions (all)")
    myAx.plot(np.array(Zs)*1E4,np.array(nSurv_start)*100,"r.",label="no interactions (wrt entrance)")
    myAx.plot(myXs,fitVals,"g-",label=r"fit: $\lambda$=%g $\mu$m"%(lambdaFit))
    myAx.grid()
    myAx.set(ylabel="[%]",xlabel=r"thickness [$\mu$m]")
    myAx.legend()
    plt.yscale("log")
    fig.suptitle(myTitle)
    plt.tight_layout()
    plt.show()
    print(fitParameters)
    print(covariance)
    # exit()

    # --------------------------------------------------------
    # analysis of parsed distributions
    # # - scatter plot
    # print("Scatter plot of original distribution...")
    # DataSet.ShowScatter(Zs=Zs,LabZs=LabZs,figName="%s_orig"%(figName))
    
    # - histograms
    print("Histograms of original distribution...")
    myFigTitle=myTitle
    if (myFigTitle is not None): myFigTitle="Histograms of original distributions - %s"%(myTitle)
    myFigName=figName
    if (myFigName is not None): myFigName="%s_orig"%(figName)
    myPicName=picName
    if (myPicName is not None): myPicName="%s_orig"%(picName)
    DataSet.HistMe(Zs=Zs,LabZs=LabZs,
                   myWhat=[["XX","XP","DE"],["YY","YP","DP"]],
                   myFacts=[[1E7,1E6,1],[1E7,1E6,1]],
                   myUnits=[["nm",r"$\mu$rad","MeV"],["nm",r"$\mu$rad","MeV/c"]],
                   figTitle=myFigTitle,figName=myFigName,picName=myPicName)

    # - histograms with fits
    myFigTitle=myTitle
    if (myFigTitle is not None): myFigTitle="Fits of original distributions - %s"%(myTitle)
    myFigName=figName
    if (myFigName is not None): myFigName="%s_orig_fits"%(figName)
    myPicName=picName
    if (myPicName is not None): myPicName="%s_orig_fits"%(picName)
    CumulHist, CumulMids, CumulFitParameters = \
               DataSet.HistMe(Zs=Zs,LabZs=LabZs,
                              lSkipZeros=True,
                              myWhat=[["XP","YP","TH"]],
                              myFacts=[[1E6,1E6,1E6]],
                              myUnits=[[r"$\mu$rad",r"$\mu$rad",r"$\mu$rad"]],
                              fitModel=[["Gauss","Gauss","lGauss"]],
                              figTitle=myFigTitle,figName=myFigName,picName=myPicName)

    # --------------------------------------------------------
    # analysis of differential distributions
    # # - scatter plot
    # print("Scatter plot of diff distribution...")
    # DiffDataSet.ShowScatter(Zs=Zs,LabZs=LabZs,
    #             myWhat=[["XP","YP"],["DE","DP"]],
    #             myFacts=[[1000,1000],[1,1]],
    #             myUnits=[["mrad","mrad"],["MeV","MeV/c"]],
    #             figName="%s_diff"%(figName),nCols=2)
    # - histograms with fits
    myFigTitle=myTitle
    if (myFigTitle is not None): myFigTitle="Fits of single layer angular distributions - %s"%(myTitle)
    myFigName=figName
    if (myFigName is not None): myFigName="%s_singLayer_fits"%(figName)
    myPicName=picName
    if (myPicName is not None): myPicName="%s_singLayer_fits"%(picName)
    print("Histograms of diff distribution...")
    DiffHist, DiffMids, DiffFitParameters = \
        DiffDataSet.HistMe(Zs=Zs,LabZs=LabZs,
                           lSkipZeros=True,
                           myWhat=[["XP","YP","TH"]],
                           myFacts=[[1E6,1E6,1E6]],
                           myUnits=[[r"$\mu$rad",r"$\mu$rad",r"$\mu$rad"]],
                           fitModel=[["Gauss","Gauss","lGauss"]],
                           figTitle=myFigTitle,figName=myFigName,picName=myPicName)
    
#     DiffDataSet.HistMe(Zs=Zs[1:3],LabZs=LabZs[1:3],
#                        lSkipZeros=True,
#                        myWhat=[["XP","YP"]],
#                        myFacts=[[1E6,1E6]],
#                        myUnits=[[r"$\mu$rad",r"$\mu$rad"]],
#                        fitModel="Gauss")
#                 figName="%s_diff_XPYP"%(figName),picName="%s_diff_XPYP"%(picName),
#     DiffDataSet.HistMe(Zs=Zs,LabZs=LabZs,
#                 myWhat=[["DE","DP"]],
#                 myFacts=[[1,1]],
#                 myUnits=[["MeV","MeV/c"]])
# #                 figName="%s_diff_DEDP"%(figName),picName="%s_diff_DEDP"%(picName),
# #                 fitModel="LanGau")

    # compare against predictions by Lynch & Dahl
    myWhat=["XP","YP","TH"]
    myFacts=[1E6,1E6,1E6]
    myUnits=[r"$\mu$rad",r"$\mu$rad",r"$\mu$rad"]
    nCols=len(myWhat)
    fig, axs = plt.subplots(1,nCols,figsize=(5*nCols,5))
    for iWhat in range(len(myWhat)):
        if (nCols==1):
            myAx=axs
        else:
            myAx=axs[iWhat]
        FitParam=[CumulFitParameters[myWhat[iWhat]][ii][0] for ii in range(len(Zs))]
        myAx.plot(np.array(Zs)*1E4,FitParam,"o",label="Cumulative scatter")
        FitParam=[DiffFitParameters[myWhat[iWhat]][ii][0] for ii in range(len(Zs))]
        myAx.plot(np.concatenate((np.array(Zs)[0],np.diff(Zs)))*1E4,FitParam,"o",label="Single-layer scatter")
        myXs=np.linspace(Zs[0],Zs[-1],num=200)
        myAx.plot(myXs*1E4,LynchDahl(myXs,BeamPart.betaRel,BeamPart.pc,BeamPart.Z,X0,lLog=True)*myFacts[iWhat],
                  "r--",label="Lynch&Dahl")
        myAx.plot(myXs*1E4,LynchDahl(myXs,BeamPart.betaRel,BeamPart.pc,BeamPart.Z,X0,lLog=False)*myFacts[iWhat],
                  "m--",label="Lynch&Dahl (no Log)")
        # myAx.plot(myXs*1E4,LynchDahlBetter(myXs,betaRel,pc0,Z,Zmat,Amat,rho)*myFacts[iWhat],"g--",label="Lynch&Dahl (better)")
        myAx.grid()
        myAx.set(ylabel="%s [%s]"%(myWhat[iWhat],myUnits[iWhat]),
                     xlabel=r"thickness [$\mu$m]"
        )
        myAx.legend()
    plt.tight_layout()
    plt.show()

    print("...done;")

def CompareHistograms(mids,hist,histLabels,myWhat,myUnit,fitModels=None,fitParams=None,fitMids=None,fitVals=None,axs=None):
    if (type(mids[0]) is not list): mids=[mids]
    if (type(hist[0]) is not list): hist=[hist]
    if (type(histLabels) is string): histLabels=[histLabels]
    if (fitModels is None):
        fitModels=[ None for ii in range(len(hist)) ]
    elif (type(fitModels[0][0]) is not list):
        fitModels=[fitModels]
    if (fitParams is None):
        fitParams=[ None for ii in range(len(hist)) ]
    elif (type(fitParams[0][0]) is not list):
        fitParams=[fitParams]
    if (fitMids is None):
        fitMids=[ None for ii in range(len(hist)) ]
    elif (type(fitMids[0][0]) is not list):
        fitMids=[fitMids]
    if (fitVals is None):
        fitVals=[ None for ii in range(len(hist)) ]
    elif (type(fitVals[0][0]) is not list):
        fitVals=[fitVals]

    if (axs is None):
        nCols=1; nRows=1
        fig, axs = plt.subplots(nRows,nCols,figsize=(10*nCols,5*nRows))
    else:
        fig = None
    for myMids,myHist,myHistLabel,myFitModels,myFitParams,myFitMids,myFitVals in \
        zip(mids,hist,histLabels,fitModels,fitParams,fitMids,fitVals):
        axs=showHist(myMids,myHist,myHistLabel,myFitModels,myFitParams,myFitMids,myFitVals,myUnit,ax)
    axs.set(xlabel="%s [%s]"%(myWhat,myUnit),ylabel="counts")
    return axs, fig

def showHist(mids,hist,histLabel,fitModels,fitParams,fitMids,fitVals,myUnit,ax,myC=None):
    '''
    * histograms:
      - mids,hist: 1D-array;
      - histLabel: string
    * fits:
      - fitModels: 1D-array of string;
      - fitParams: 1D-array of 1D-array (float);
      - fitMids,fitVals: 1D-array of 1D-array (float);
    * myUnit: string;
    * ax: matplotlib axis;
    '''
    if (myC is None):
        ax.step(mids,hist,where='mid',label=histLabel)
    else:
        ax.step(mids,hist,myC,where='mid',label=histLabel)
    if (fitModels is not None):
        for iFit,fitModel in enumerate(fitModels):
            if (fitModel.upper().startswith("GAUSS")):
                myLab=r"GAUSS fit($\mu$=%g %s,$\sigma$=%g %s)"%(
                          fitParams[iFit][2],myUnit,
                          fitParams[iFit][0],myUnit)
            elif (fitModel.upper().startswith("LGAUSS")):
                myLab=r"lGAUSS fit($\sigma$=%g %s)"%(
                          fitParams[iFit][0],myUnit)
            elif (fitModel.upper().startswith("RUTH")):
                myLab=r"RUTH fit"
            elif (fitModel.upper().startswith("LANGAU")):
                myLab=r"LANGAU fit($\mu$=%g %s,$\sigma$=%g %s)"%(
                          fitParams[iFit][0],myUnit,
                          fitParams[iFit][2],myUnit)
            elif (fitModel.upper().startswith("LANDAU")):
                myLab=r"LANDAU fit($\mu$=%g %s,$\eta$=%g)"%(
                          fitParams[iFit][0],myUnit,fitParams[iFit][1])
            ax.plot(fitMids[iFit],fitVals[iFit],'--',label=myLab)
    return ax
                    
if ( __name__ == "__main__" ):
    # myCase="test_CMax"
    # Zs=[0.0,1.50E-4,1.55E-4,3.05E-4,3.10E-4,4.60E-4,4.65E-4,6.15E-4,6.20E-4] # [cm]
    # myCase="test_CMax_LIGHTER-GOLD"
    # Zs=[0.0,1.50E-4,2.50E-4,4.00E-4,5.00E-4,6.50E-4,7.50E-4,9.00E-4,1.00E-3] # [cm]
    # myCase="test_CMax_onlyKAPTON_MCS"
    # myCase="test_CMax_onlyKAPTON_SS"
    # myCase="test_CMax_onlyLGOLD_MCS"
    # myCase="test_CMax_onlyLGOLD_SS"
    # myCase="test_CMax_onlyGOLD_SS"
    # Zs=[0.0, 0.25E-4, 0.75E-4, 1.75E-4, 3.75E-4, 7.75E-4,14.75E-4,23.75E-4,35.75E-4] # [cm]
    # Zs=[0.0,  12.5E-7,  37.5E-7,  87.5E-7, 187.5E-7, 387.5E-7, 737.5E-7,1187.5E-7,1787.5E-7] # [cm]
    # myTitle="only GOLD - SS"
    # # myCase="test_CMax_onlyGOLD_SS"
    # # Zs=[0.0,  12.5E-7,  37.5E-7,  87.5E-7, 187.5E-7, 387.5E-7, 737.5E-7,1187.5E-7,1787.5E-7] # [cm]
    # # myCase="test_CMax_incrGOLD"
    # # Zs=[0.0,1.50E-4,1.55E-4,3.05E-4,3.15E-4,4.65E-4,4.80E-4,6.30E-4,6.50E-4] # [cm]
    # myPath="./%s"%(myCase)
    # fort23FileName="%s/run_00002/pepites_mcs001_fort.23.gz"%(myPath)
    # m0, A, Z, betaRel, gammaRel, pc0, Ek0 = SetBeamPart("Carbon",Eku=398.84) # [MeV/u]
    # X0=28.57  # [cm] KAPTON 
    # X0=6.461  # [cm] LIGHT GOLD
    # X0=0.3344 # [cm] GOLD
    # Zmat=79
    # Amat=197
    # rho=19.32
    # 
    # LabZs=["at entrance",
    #        "1st kapton layer",
    #        "1st gold layer",
    #        "2nd kapton layer",
    #        "2nd gold layer",
    #        "3rd kapton layer",
    #        "3rd gold layer",
    #        "4th kapton layer",
    #        "4th gold layer"    ]
    # LabZs=["at entrance",
    #        "1st layer",
    #        "2nd layer",
    #        "3rd layer",
    #        "4th layer",
    #        "5th layer",
    #        "6th layer",
    #        "7th layer",
    #        "8th layer"
    #        ]
    # figName="%s/pepites_mcs"%(myPath)
    # picName="%s/pepites_mcs"%(myPath)

    myPath="./pMax_onlyGOLD_SS"
    figName=None # figName="%s/pepites_mcs"%(myPath)
    picName=None # picName="%s/pepites_mcs"%(myPath)
    
    userInput=myPath.split("/")[-1].split("_")
    if (userInput[0]=="pMax"):
        BeamPart=BeamPart("Proton",228.57)
    elif (userInput[0]=="pMin"):
        BeamPart=BeamPart("Proton",62.73)
    elif (userInput[0]=="CMax"):
        BeamPart=BeamPart("Carbon",398.84,lPerNucleon=True)
    elif (userInput[0]=="CMin"):
        BeamPart=BeamPart("Carbon",115.23,lPerNucleon=True)
    else:
        print("...unknown element %s in case %s!"%(userInput[0],userInput))
        exit(1)
    if (userInput[1]=="onlyGOLD"):
        X0=0.3344 # Radiation Length [cm]
        Zs=[0.0,  12.5E-7,  37.5E-7,  87.5E-7, 187.5E-7, 387.5E-7, 737.5E-7,1187.5E-7,1787.5E-7] # [cm]
    elif (userInput[1]=="onlyKAPTON"):
        X0=28.57  # Radiation Length [cm]
        # Zs=[0.0,  0.25E-4,  0.75E-4,  1.75E-4,  3.75E-4,  7.75E-4, 14.75E-4, 23.75E-4, 35.75E-4] # [cm]
        Zs=[0.0,  8.46E-7, 25.39E-7, 59.24E-7, 126.94E-7,  262.35E-7,   499.31E-7,   803.98E-7,  1210.20E-7 ] # [cm]
    elif (userInput[1]=="onlyLightGOLD"):
        X0=6.461  # Radiation Length [cm] lightGOLD
        Zs=[0.0,  0.25E-4,  0.75E-4,  1.75E-4,  3.75E-4,  7.75E-4, 14.75E-4, 23.75E-4, 35.75E-4] # [cm]
    else:
        print("...unknown element %s in case %s!"%(userInput[1],userInput))
        exit(1)
    ZsDiffs=np.zeros(len(Zs))
    ZsDiffs[1:]=np.diff(Zs)
    myTitle=""
    for ii in range(len(userInput)):
        myTitle=myTitle+"{0[%d]}"%(ii)
        if (len(userInput)>1 and ii<len(userInput)-1):
            myTitle=myTitle+" - "
    myTitle=myTitle.format(userInput)
    LabZs=["at entrance",
           "end of 1st layer",
           "end of 2nd layer",
           "end of 3rd layer",
           "end of 4th layer",
           "end of 5th layer",
           "end of 6th layer",
           "end of 7th layer",
           "end of 8th layer"
           ]
    for ii in range(1,len(LabZs)):
        LabZs[ii]="%s - Dz=%.3g um"%(LabZs[ii],ZsDiffs[ii]*1E4)
    fort23FileName="%s/run_00001/pepites_mcs001_fort.23.gz"%(myPath)
    
    # get data set
    DataSet=Fort23Data.FromFile(fort23FileName)
    DataSet.RoundZ(Zs)
    DataSet.RemoveDieing()
    DataSet.GetAngles()
    DataSet.GetDeltaEnergy(BeamPart.Ek,BeamPart.mm)
    DiffDataSet=DataSet.GetDifferentialDataSet()
    
    # StandardAnalysis(fort23FileName,Zs,LabZs,BeamPart,X0,figName=figName,picName=picName,myTitle=myTitle)

    # single layer analysis
    hist={}; mids={}; fitParams={}; fitVals={}; fitMids={}; myStats={}
    myWhats={}; myCuts={}; myFitModels={}; nBins={}
    iZs=[]

    iZ=1; iZs.append(iZ)
    myWhats[LabZs[iZ]]=["TH","DE"]
    myCuts[LabZs[iZ]]={ "TH":[[0,50],[70,600]], "DE": [[-0.06,-0.006]] }
    myFitModels[LabZs[iZ]]={ "TH":["lGauss","Ruth"], "DE": ["landau"] }
    nBins[LabZs[iZ]]={"TH":400,"DE":600}
    iZ=2; iZs.append(iZ)
    myWhats[LabZs[iZ]]=["TH","DE"]
    myCuts[LabZs[iZ]]={ "TH":[[0,50],[70,600]], "DE": [[-0.2,0]] }
    myFitModels[LabZs[iZ]]={ "TH":["lGauss","Ruth"], "DE": ["landau"] }
    nBins[LabZs[iZ]]={"TH":400,"DE":600}
    iZ=3; iZs.append(iZ)
    myWhats[LabZs[iZ]]=["TH","DE"]
    myCuts[LabZs[iZ]]={ "TH":[[0,50],[70,600]], "DE": [[-0.4,0],[-0.85,-0.71]] }
    myFitModels[LabZs[iZ]]={ "TH":["lGauss","Ruth"], "DE": ["landau","langau"] }
    nBins[LabZs[iZ]]={"TH":400,"DE":600}
    iZ=4; iZs.append(iZ)
    myWhats[LabZs[iZ]]=["TH","DE"]
    myCuts[LabZs[iZ]]={ "TH":[[0,60],[100,600]], "DE": [[-0.5,0],[-0.9,-0.70]] }
    myFitModels[LabZs[iZ]]={ "TH":["lGauss","Ruth"], "DE": ["landau","langau"] }
    nBins[LabZs[iZ]]={"TH":400,"DE":600}
    iZ=5; iZs.append(iZ)
    myWhats[LabZs[iZ]]=["TH","DE"]
    myCuts[LabZs[iZ]]={ "TH":[[0,110],[150,800]], "DE": [[-0.7,0],[-1.15,-0.82]] }
    myFitModels[LabZs[iZ]]={ "TH":["lGauss","Ruth"], "DE": ["landau","langau"] }
    nBins[LabZs[iZ]]={"TH":400,"DE":600}
    iZ=6; iZs.append(iZ)
    myWhats[LabZs[iZ]]=["TH","DE"]
    myCuts[LabZs[iZ]]={ "TH":[[0,160],[250,900]], "DE": [[-0.8,0],[-1.45,-0.9]] }
    myFitModels[LabZs[iZ]]={ "TH":["lGauss","Ruth"], "DE": ["landau","langau"] }
    nBins[LabZs[iZ]]={"TH":400,"DE":600}

    for iZ in iZs:
        hist[LabZs[iZ]], mids[LabZs[iZ]], fitParams[LabZs[iZ]], fitVals[LabZs[iZ]], fitMids[LabZs[iZ]], myStats[LabZs[iZ]] = \
            singleLayerAnalysis(DiffDataSet,Zs[iZ],LabZs[iZ],figName=figName,picName=picName,myTitle=myTitle, \
             whats=myWhats[LabZs[iZ]],cuts=myCuts[LabZs[iZ]],fitModels=myFitModels[LabZs[iZ]],nBins=nBins[LabZs[iZ]])

    nCols=len(["TH","DE"]); nRows=1
    fig, axs = plt.subplots(nRows,nCols,figsize=(10*nCols,5*nRows))
    for jj,myWhat in enumerate(["TH","DE"]):
        if (myWhat in ["TH","XP","YP"]):
            myFact=1E6; myUnit=r"$\mu$rad"
        elif (myWhat in ["DE"]):
            myFact=1E3; myUnit="keV"
        for iZ in iZs:
            axs[jj]=showHist(mids[LabZs[iZ]][myWhat],hist[LabZs[iZ]][myWhat],LabZs[iZ],
                             myFitModels[LabZs[iZ]][myWhat],fitParams[LabZs[iZ]][myWhat],
                             fitMids[LabZs[iZ]][myWhat],fitVals[LabZs[iZ]][myWhat],
                             myUnit,axs[jj])
        axs[jj].grid(); axs[jj].legend()
        axs[jj].set(xlabel="%s [%s]"%(myWhat,myUnit),ylabel="counts")
        plt.tight_layout()
    
        
    plt.show()

import struct
import numpy as np

class ScoringMesh():
    def __init__(self,lDebug=True):
        if (lDebug): print("new mesh...")
        return

    def parseBin_Header(self,ff,lDebug=True):
        'usbrea.f naming convention'
        print("parsing BIN header...")
        self.MB=struct.unpack("I",ff.read(4))[0]
        self.TITUSB=str(ff.read(10),"utf-8")
        self.ITUSBN=struct.unpack("I",ff.read(4))[0]
        self.IDUSBN=struct.unpack("I",ff.read(4))[0]
        #
        self.XLOW  =struct.unpack("f",ff.read(4))[0]
        self.XHIGH =struct.unpack("f",ff.read(4))[0]
        self.NXBIN =struct.unpack("I",ff.read(4))[0]
        self.DXUSBN=struct.unpack("f",ff.read(4))[0]
        #
        self.YLOW  =struct.unpack("f",ff.read(4))[0]
        self.YHIGH =struct.unpack("f",ff.read(4))[0]
        self.NYBIN =struct.unpack("I",ff.read(4))[0]
        self.DYUSBN=struct.unpack("f",ff.read(4))[0]
        #
        self.ZLOW  =struct.unpack("f",ff.read(4))[0]
        self.ZHIGH =struct.unpack("f",ff.read(4))[0]
        self.NZBIN =struct.unpack("I",ff.read(4))[0]
        self.DZUSBN=struct.unpack("f",ff.read(4))[0]
        #
        self.LNTZER=struct.unpack("?",ff.read(1))[0]
        self.BKUSBN=struct.unpack("f",ff.read(4))[0]
        self.B2USBN=struct.unpack("f",ff.read(4))[0]
        self.TCUSBN=struct.unpack("f",ff.read(4))[0]
        if (lDebug):
            for key,val in self.__dict__.items():
                print("-->",key,":",val)

    def saveBin_Header(self,ff,lDebug=True):
        print("saving BIN header...")
        ff.write(struct.pack("I",self.MB))
        ff.write(struct.pack("10s",bytes(self.TITUSB,"utf-8")))
        ff.write(struct.pack("I",self.ITUSBN))
        ff.write(struct.pack("I",self.IDUSBN))
        #
        ff.write(struct.pack("f",self.XLOW  ))
        ff.write(struct.pack("f",self.XHIGH ))
        ff.write(struct.pack("I",self.NXBIN ))
        ff.write(struct.pack("f",self.DXUSBN))
        #
        ff.write(struct.pack("f",self.YLOW  ))
        ff.write(struct.pack("f",self.YHIGH ))
        ff.write(struct.pack("I",self.NYBIN ))
        ff.write(struct.pack("f",self.DYUSBN))
        #
        ff.write(struct.pack("f",self.ZLOW  ))
        ff.write(struct.pack("f",self.ZHIGH ))
        ff.write(struct.pack("I",self.NZBIN ))
        ff.write(struct.pack("f",self.DZUSBN))
        #
        ff.write(struct.pack("?",self.LNTZER))
        ff.write(struct.pack("f",self.BKUSBN))
        ff.write(struct.pack("f",self.B2USBN))
        ff.write(struct.pack("f",self.TCUSBN))

    def copyHeader(self,otherBin,copyKeys=[\
                                           "MB", "TITUSB", "ITUSBN", "IDUSBN", \
                                           "XLOW","XHIGH","NXBIN","DXUSBN", \
                                           "YLOW","YHIGH","NYBIN","DYUSBN", \
                                           "ZLOW","ZHIGH","NZBIN","DZUSBN", \
                                           "LNTZER", "BKUSBN", "B2USBN", "TCUSBN" ]):
        for myKey in copyKeys:
            self.__dict__[myKey]=otherBin.__dict__[myKey]

    def parseBin_Junk(self,ff,where,lDebug=True):
        if (where==0):
            # junk=ff.read(11) # b'\xf2IqV\x00\x00\x00\xe8\x03\x00\x00'
            # if (lDebug): print(junk)
            junk=ff.read(7) #
            if (lDebug): print("--> junk:",junk)
            self.junk=[]
            self.junk.append(ff.read(2))
            if (lDebug): print("--> self.junk[%i]:"%(where),self.junk[where])
            junk=ff.read(2)
            if (lDebug): print("--> junk:",junk)
        elif (where==1):
            # junk=ff.read(8) # b'\xe8\x03\x00\x00\x0e\x00\x00\x00'
            # if (lDebug): print(junk)
            self.junk.append(ff.read(2))
            if (lDebug): print("--> self.junk[%i]:"%(where),self.junk[where])
            junk=ff.read(6)
            if (lDebug): print("--> junk:",junk)
        elif (where==2):
            # junk=ff.read(8) # b'\x0e\x00\x00\x00\xe8\x03\x00\x00'
            # if (lDebug): print(junk)
            junk=ff.read(4) #
            if (lDebug): print("--> junk:",junk)
            self.junk.append(ff.read(2))
            if (lDebug): print("--> self.junk[%i]:"%(where),self.junk[where])
            junk=ff.read(2)
            if (lDebug): print("--> junk:",junk)
        elif (where==3):
            # junk=ff.read() # b'\xe8\x03\x00\x00'
            # if (lDebug): print("--> junk :",junk)
            self.junk.append(ff.read(2))
            if (lDebug): print("--> self.junk[%i]:"%(where),self.junk[where])
            junk=ff.read(2)
            if (lDebug): print("--> junk:",junk)
    
    def saveBin_Junk(self,ff,where,lDebug=True):
        if (where==0):
            # ff.write(b'\xf2IqV\x00\x00\x00\xe8\x03\x00\x00')
            ff.write(b'\xf2IqV\x00\x00\x00')
            ff.write(self.junk[where])
            ff.write(b'\x00\x00')
        elif (where==1):
            # ff.write(b'\xe8\x03\x00\x00\x0e\x00\x00\x00')
            ff.write(self.junk[where])
            ff.write(b'\x00\x00\x0e\x00\x00\x00')
        elif (where==2):
            # ff.write(b'\x0e\x00\x00\x00\xe8\x03\x00\x00')
            ff.write(b'\x0e\x00\x00\x00')
            ff.write(self.junk[where])
            ff.write(b'\x00\x00')
        elif (where==3):
            # ff.write(b'\xe8\x03\x00\x00')
            ff.write(self.junk[where])
            ff.write(b'\x00\x00')
    
    def parseBin_Data(self,ff,lDebug=True,nShow=10):
        print("parsing BIN data...")
        iterator=struct.iter_unpack("f",ff.read(self.NXBIN*self.NYBIN*self.NZBIN*4))
        self.GMSTOR=np.array([ val[0] for val in iterator ])
        if (lDebug):
            print("...acquired %d values!"%(len(self.GMSTOR)))
            print("   first %d:"%(nShow),self.GMSTOR[0:nShow])
            print("   last  %d:"%(nShow),self.GMSTOR[-nShow:])

    def saveBin_Data(self,ff,lDebug=True):
        print("saving BIN data...")
        ff.write(struct.pack("%sf"%len(self.GMSTOR),*self.GMSTOR))

    def parseBin_Stats(self,ff,lDebug=True,nShow=10):
        print("parsing BIN stats...")
        iterator=struct.iter_unpack("f",ff.read(self.NXBIN*self.NYBIN*self.NZBIN*4))
        self.GBSTOR=np.array([ val[0] for val in iterator ])
        if (lDebug):
            print("...acquired %d values!"%(len(self.GBSTOR)))
            print("   first %d [%%]:"%(nShow),self.GBSTOR[0:nShow]*100)
            print("   last  %d [%%]:"%(nShow),self.GBSTOR[-nShow:]*100)

    def saveBin_Stats(self,ff,lDebug=True):
        print("saving BIN stats...")
        ff.write(struct.pack("%sf"%len(self.GBSTOR),*self.GBSTOR))

    def reshape(self,lDebug=True):
        print("reshaping data...")
        if ("GMSTOR" not in self.__dict__.keys()):
            print("...no data loaded yet!")
            exit(1)
        if (lDebug): print("...original GMSTOR shap:",self.GMSTOR.shape)
        self.GMSTOR=np.reshape(self.GMSTOR,(self.NXBIN,self.NYBIN,self.NZBIN),order="F")
        if (lDebug): print("...GMSTOR reshaped to:",self.GMSTOR.shape)
        if ("GBSTOR" not in self.__dict__.keys()):
            print("...no stats to reshape!!")
            exit(1)
        if (lDebug): print("...original GBSTOR shap:",self.GBSTOR.shape)
        self.GBSTOR=np.reshape(self.GBSTOR,(self.NXBIN,self.NYBIN,self.NZBIN),order="F")
        if (lDebug): print("...GBSTOR reshaped to:",self.GBSTOR.shape)
        print("...done;")

    def cumsum(self,iAxis=2,lDebug=True):
        print("cumulative sums along axis=%d (python numbering)..."%(iAxis))
        if (self.ITUSBN!=10):
            print("sorry: only cartesian scorings for the time being!")
            exit(1)
        myCumSum=np.cumsum(self.GMSTOR,axis=iAxis)
        if (lDebug): print("...shape of cumsum:",myCumSum.shape)
        print("...done;")
        return myCumSum

    def sum(self,iAxis=2,lDebug=True):
        print("summing values along axis=%d (python numbering)..."%(iAxis))
        if (self.ITUSBN!=10):
            print("sorry: only cartesian scorings for the time being!")
            exit(1)
        myVals=np.sum(self.GMSTOR,axis=iAxis)
        if (lDebug): print("...shape of sum:",myVals.shape)
        print("...done;")
        return myVals

    def retMesh(self,iAxis=-2,lDebug=True):
        if (iAxis==0):
            return np.linspace(self.XLOW,self.XHIGH,self.NXBIN+1)
        elif (iAxis==1):
            return np.linspace(self.YLOW,self.YHIGH,self.NYBIN+1)
        elif (iAxis==2):
            return np.linspace(self.ZLOW,self.ZHIGH,self.NZBIN+1)
        else:
            print("please choose an axis with python indexing!")
            exit(1)

    @staticmethod
    def parseBin_Single(ff,lDebug=True):
        new=ScoringMesh(lDebug=lDebug)
        new.parseBin_Header(ff,lDebug=lDebug)
        new.parseBin_Junk(ff,where=0,lDebug=lDebug)
        new.parseBin_Data(ff,lDebug=lDebug)
        new.parseBin_Junk(ff,where=1,lDebug=lDebug)
        junk=ff.read(10) # STATISTICS
        if (lDebug): print(junk)
        junk=struct.unpack("I",ff.read(4))[0] # IB
        if (lDebug): print(junk)
        new.parseBin_Junk(ff,where=2,lDebug=lDebug)
        new.parseBin_Stats(ff,lDebug=lDebug)
        new.parseBin_Junk(ff,where=3,lDebug=lDebug)
        print("...acquired;")
        return new

    def saveBin_Single(self,ff,lDebug=True):
        print("saving mesh values...")
        self.saveBin_Header(ff,lDebug=lDebug)
        self.saveBin_Junk(ff,0,lDebug=lDebug)
        self.saveBin_Data(ff,lDebug=lDebug)
        self.saveBin_Junk(ff,1,lDebug=lDebug)
        ff.write(struct.pack("10s",bytes("STATISTICS","utf-8")))
        ff.write(struct.pack("I",1 ))
        self.saveBin_Junk(ff,2,lDebug=lDebug)
        self.saveBin_Stats(ff,lDebug=lDebug)
        self.saveBin_Junk(ff,3,lDebug=lDebug)
        print("...saved;")

    @staticmethod
    def Bin_sum(aBin,bBin,aa=1,bb=1,MB=1,TITUSB="sum",ITUSBN=1,IDUSBN=0,lDebug=True):
        if (not aBin.Bin_is_compatible(bBin,lDebug=lDebug)):
            print("NOT compatible bins!")
            exit(1)
        new=ScoringMesh(lDebug=lDebug)
        new.copyHeader(aBin)
        new.IDUSBN=IDUSBN
        new.junk=aBin.junk
        new.GMSTOR=np.add(aa*aBin.GMSTOR,bb*bBin.GMSTOR)
        # sig_sum=sqrt(aa**2*sig_a**2+bb**2*sig_b**2)
        new.GBSTOR=np.divide(np.sqrt(np.add((aa*np.multiply(aBin.GMSTOR,aBin.GBSTOR))**2,(bb*np.multiply(bBin.GMSTOR,bBin.GBSTOR))**2)),new.GMSTOR)
        return new

    @staticmethod
    def Bin_dif(aBin,bBin,aa=1,bb=1,MB=1,TITUSB="dif",ITUSBN=1,IDUSBN=0,lDebug=True):
        return ScoringMesh.Bin_sum(aBin,bBin,aa=aa,bb=-bb,MB=MB,TITUSB=TITUSB,ITUSBN=ITUSBN,IDUSBN=IDUSBN,lDebug=lDebug)

    @staticmethod
    def Bin_mul(aBin,bBin,aa=1,bb=1,MB=1,TITUSB="mul",ITUSBN=1,IDUSBN=0,lDebug=True):
        if (not aBin.Bin_is_compatible(bBin,lDebug=lDebug)):
            print("NOT compatible bins!")
            exit(1)
        new=ScoringMesh(lDebug=lDebug)
        new.copyHeader(aBin)
        new.IDUSBN=IDUSBN
        new.junk=aBin.junk
        new.GMSTOR=np.multiply(aa*aBin.GMSTOR,bb*bBin.GMSTOR)
        # sig_sum/sum=sqrt((sig_a/AA)**2+(sig_b/BB)**2)
        new.GBSTOR=np.sqrt(np.add(aBin.GBSTOR**2,bBin.GBSTOR**2))
        return new

    @staticmethod
    def Bin_div(aBin,bBin,aa=1,bb=1,MB=1,TITUSB="div",ITUSBN=1,IDUSBN=202,lDebug=True):
        '''element-wise division: aa*aBin/(bb*bBin)'''
        if (not aBin.Bin_is_compatible(bBin,lDebug=lDebug)):
            print("NOT compatible bins!")
            exit(1)
        new=ScoringMesh(lDebug=lDebug)
        new.copyHeader(aBin)
        new.IDUSBN=IDUSBN
        new.junk=aBin.junk
        doseMax=np.max(bBin.GMSTOR)
        new.GMSTOR=np.divide(aa*aBin.GMSTOR,bb*bBin.GMSTOR,out=np.zeros(np.shape(bBin.GMSTOR)),where=bBin.GMSTOR>=0.01*doseMax)
        # new.GMSTOR=np.divide(aa*aBin.GMSTOR,bb*bBin.GMSTOR,out=aa*aBin.GMSTOR,where=bBin.GMSTOR!=0.0)
        # sig_sum/sum=sqrt((sig_a/AA)**2+(sig_b/BB)**2)
        new.GBSTOR=np.sqrt(np.add(aBin.GBSTOR**2,bBin.GBSTOR**2))
        return new

    def Bin_is_compatible(self,newBin,lDebug=True):
        compatibleBins=True
        checkKeys=["XLOW","XHIGH","NXBIN","DXUSBN", \
                   "YLOW","YHIGH","NYBIN","DYUSBN", \
                   "ZLOW","ZHIGH","NZBIN","DZUSBN", \
                   "BKUSBN","B2USBN" \
        ]
        for myKey in checkKeys:
            lEq=(self.__dict__[myKey]==newBin.__dict__[myKey])
            if (lDebug and not lEq): print("self.%s!=newBin.%s"%(myKey,myKey))
            compatibleBins=compatibleBins and lEq
        return compatibleBins
        
class BINOUTPUT():

    def __init__(self):
        'usbrea.f naming convention'
        self.RUNTIT=""
        self.RUNTIM=""
        self.WEIPRI=0.0
        self.NCASE =0
        self.MCASE =0
        self.NBATCH=0
        self.meshes=[]
        return

    def parseHeader(self,ff,lDebug=True):
        print("parsing header...")
        self.RUNTIT=str(ff.read(80),"utf-8")
        self.RUNTIM=str(ff.read(32),"utf-8")
        self.WEIPRI=struct.unpack("f",ff.read(4))[0]
        self.NCASE =struct.unpack("I",ff.read(4))[0]
        self.MCASE =struct.unpack("I",ff.read(4))[0]
        self.NBATCH=struct.unpack("I",ff.read(4))[0]
        if (lDebug):
            for key,val in self.__dict__.items():
                print("-->",key,":",val)

    def saveHeader(self,ff,lDebug=True):
        print("saving header...")
        ff.write(struct.pack("80p",bytes(self.RUNTIT,"utf-8")))
        ff.write(struct.pack("32p",bytes(self.RUNTIM,"utf-8")))
        ff.write(struct.pack("f",self.WEIPRI))
        ff.write(struct.pack("I",self.NCASE))
        ff.write(struct.pack("I",self.MCASE))
        ff.write(struct.pack("I",self.NBATCH))

    @staticmethod
    def parseBin(iFileName,lDebug=True):
        print("parsing file %s..."%(iFileName))
        new=BINOUTPUT()
        ff=open(iFileName,"rb")
        junk=ff.read(4) # b'\x80\x00\x00\x00'
        if (lDebug): print("--> junk :",junk)
        new.parseHeader(ff,lDebug=lDebug)
        junk=ff.read(8) # b'\x80\x00\x00\x00V\x00\x00\x00'
        if (lDebug): print("--> junk :",junk)
        new.meshes.append(ScoringMesh.parseBin_Single(ff,lDebug=lDebug))
        ff.close()
        print("...done;")
        return new

    def save(self,oFileName,lDebug=True):
        print("saving to file %s..."%(oFileName))
        ff=open(oFileName,"wb")
        ff.write(b'\x80\x00\x00\x00')
        self.saveHeader(ff,lDebug=lDebug)
        ff.write(b'\x80\x00\x00\x00V\x00\x00\x00')
        self.meshes[0].saveBin_Single(ff,lDebug=lDebug)
        ff.close()
        print("...done;")

def testParsing(iFileName,oFileName,lDebug=True):
    BNN=BINOUTPUT.parseBin(iFileName,lDebug=lDebug)
    BNN.save(oFileName,lDebug=lDebug)

def save2Dmap(oFileName,XXX,YYY,myMap,lDebug=True):
     print("saving 2D map in file %s (gnuplot format)..."%(oFileName))
     if (lDebug):
         print("...shape of map:",myMap.shape)
         print("...shape of X and Y arrays:",XXX.shape,YYY.shape)
     with open(oFileName,"w") as oFile:
         for jj in range(len(YYY)-1):
             for ii in range(len(XXX)-1):
                 oFile.write(" % 17.10E % 17.10E % 12.5E\n"%(YYY[jj],XXX[ii],myMap[ii,jj]))
             oFile.write(" % 17.10E % 17.10E % 12.5E\n"%(YYY[jj],XXX[ii+1],myMap[ii,jj]))
             oFile.write("\n")
         for ii in range(len(XXX)-1):
             oFile.write(" % 17.10E % 17.10E % 12.5E\n"%(YYY[jj+1],XXX[ii],myMap[ii,jj]))
         oFile.write(" % 17.10E % 17.10E % 12.5E\n"%(YYY[jj+1],XXX[ii+1],myMap[ii,jj]))
     print("...done;")

if (__name__=="__main__"):
    from copy import deepcopy
    from scipy.interpolate import CubicSpline
    import glob
    import os
    lDebug=True
    myCase="GaussBeam_FWHM20p0_Parallel_He320p0u_RoundedPhantom"
    
    noObstacleBNN=BINOUTPUT.parseBin("./%s/GaussBeam_FWHM20p0_Parallel_He320p0u_NoTgt/FLUKA_HeCT_59.bnn"%(myCase),lDebug=lDebug)
    noObstacleBNN.meshes[0].reshape(lDebug=lDebug)
    myCumSum=noObstacleBNN.meshes[0].cumsum(lDebug=lDebug)
    XXX=noObstacleBNN.meshes[0].retMesh(iAxis=0)
    YYY=noObstacleBNN.meshes[0].retMesh(iAxis=1)
    if (lDebug):
        save2Dmap("noPhantom.dat",XXX,YYY,myCumSum[:,:,-1]/noObstacleBNN.meshes[0].NZBIN) # save 2D image of case without target
    # get inerpolation coefficients
    # ZZZ,cs=noObstacleBNN.meshes[0].getInterpoly()

    # processing each cumulative curve (a curve from each CALO channel)
    # - head zero
    # - get complement to tot energy
    iAxis=2
    myCumSum=np.insert(myCumSum,0,0.0,axis=iAxis)
    for ii in range(myCumSum.shape[0]):
        for jj in range(myCumSum.shape[1]):
            myCumSum[ii,jj,:]=myCumSum[ii,jj,-1]-myCumSum[ii,jj,:]

    for iFileName in sorted(glob.glob("./%s/%s_ang*/FLUKA_HeCT_59.bnn"%(myCase,myCase))):
        if (os.path.isfile(iFileName.replace("FLUKA_HeCT_59.bnn","mmH2O_59.dat"))):
            # skip cases already analysed
            print("... %s already processed! skipping..."%(iFileName))
            continue
        imageBNN=BINOUTPUT.parseBin(iFileName,lDebug=lDebug)
        imageBNN.meshes[0].reshape(lDebug=lDebug)
        myImage=imageBNN.meshes[0].sum(lDebug=lDebug) # from 3D mesh to 2D image (sum over long dim)

        # WET found as y=0 of shifted cumulative curves
        cumSubtracted=deepcopy(myCumSum)
        for ii in range(cumSubtracted.shape[0]):
            for jj in range(cumSubtracted.shape[1]):
                cumSubtracted[ii,jj,:]-=myImage[ii,jj]
        ZZZ=noObstacleBNN.meshes[0].retMesh(iAxis=iAxis)*10 # [cm] --> [mm]
        ZZZ-=ZZZ[0]  # USRBIN starts at z=20cm in geometry
        cs=CubicSpline(ZZZ,cumSubtracted,axis=iAxis)
        myImageMM=cs.roots(extrapolate=False)
        print("...image shape:",myImageMM.shape)
        for ii in range(myImageMM.shape[0]):
            for jj in range(myImageMM.shape[1]):
                myImageMM[ii,jj]=np.delete(myImageMM[ii,jj],myImageMM[ii,jj]<ZZZ[0])
                myImageMM[ii,jj]=np.delete(myImageMM[ii,jj],myImageMM[ii,jj]>ZZZ[-1])
                if (len(myImageMM[ii,jj])>=1):
                    myImageMM[ii,jj]=myImageMM[ii,jj][0]
                else:
                    myImageMM[ii,jj]=ZZZ[0]
        print("...cleaned image shape:",myImageMM.shape)
        save2Dmap(iFileName.replace("FLUKA_HeCT_59.bnn","mmH2O_59.dat"),XXX,YYY,myImageMM) # save 2D map
    exit()
    
    Particles=[ "ALL","PROTON","ELECTRON","4-HELIUM","12-CARBON" ]
    Particles=Particles+[ "DEUTERON","TRITON","3-HELIUM","6-LITIUM","7-LITIUM","9-BERYLLIU","10-BORON","11-BORON","14-NITROGEN","15-NITROGEN","16-OXYGEN","17-OXYGEN","18-OXYGEN","BEAMPART" ]
    Units=[ ii for ii in range(52,52+len(Particles)+1)]
    for particle,unit in zip(Particles,Units):
        print(particle,unit,"./VIF_exp_%d.bnn"%(unit),"./VIF_exp_DOSAVLET_fluscw_%s.bnn"%(particle))
    # exit(0)
   
    # testParsing("./VIF_exp_50.bnn","./test.bnn",lDebug=lDebug)

    doseBinFile="./VIF_exp_51.bnn"
    doseBNN=BINOUTPUT.parseBin(doseBinFile,lDebug=lDebug)
    for particle,unit in zip(Particles,Units):
        aletBNN=BINOUTPUT.parseBin("./VIF_exp_%d.bnn"%(unit),lDebug=lDebug)
        aletBNN.meshes[0]=ScoringMesh.Bin_div(aletBNN.meshes[0],doseBNN.meshes[0],lDebug=lDebug)
        aletBNN.save("./VIF_exp_DOSAVLET_fluscw_%s.bnn"%(particle),lDebug=lDebug)
    
    # # all particles
    # aletBNN=BINOUTPUT.parseBin("./VIF_exp_52.bnn",lDebug=lDebug)
    # aletBNN.meshes[0]=ScoringMesh.Bin_div(aletBNN.meshes[0],doseBNN.meshes[0],lDebug=lDebug)
    # aletBNN.save("./VIF_exp_DOSAVLET_fluscw_ALL.bnn",lDebug=lDebug)
    # # proton only
    # aletBNN=BINOUTPUT.parseBin("./VIF_exp_53.bnn",lDebug=lDebug)
    # aletBNN.meshes[0]=ScoringMesh.Bin_div(aletBNN.meshes[0],doseBNN.meshes[0],lDebug=lDebug)
    # aletBNN.save("./VIF_exp_DOSAVLET_fluscw_PROTON.bnn",lDebug=lDebug)
    # # electron only
    # aletBNN=BINOUTPUT.parseBin("./VIF_exp_54.bnn",lDebug=lDebug)
    # aletBNN.meshes[0]=ScoringMesh.Bin_div(aletBNN.meshes[0],doseBNN.meshes[0],lDebug=lDebug)
    # aletBNN.save("./VIF_exp_DOSAVLET_fluscw_ELECTRON.bnn",lDebug=lDebug)
    # # He-4 only
    # aletBNN=BINOUTPUT.parseBin("./VIF_exp_55.bnn",lDebug=lDebug)
    # aletBNN.meshes[0]=ScoringMesh.Bin_div(aletBNN.meshes[0],doseBNN.meshes[0],lDebug=lDebug)
    # aletBNN.save("./VIF_exp_DOSAVLET_fluscw_4-HELIUM.bnn",lDebug=lDebug)
    # # C-12 only
    # aletBNN=BINOUTPUT.parseBin("./VIF_exp_56.bnn",lDebug=lDebug)
    # aletBNN.meshes[0]=ScoringMesh.Bin_div(aletBNN.meshes[0],doseBNN.meshes[0],lDebug=lDebug)
    # aletBNN.save("./VIF_exp_DOSAVLET_fluscw_12-CARBON.bnn",lDebug=lDebug)
    # # DEUTERON only
    # aletBNN=BINOUTPUT.parseBin("./VIF_exp_57.bnn",lDebug=lDebug)
    # aletBNN.meshes[0]=ScoringMesh.Bin_div(aletBNN.meshes[0],doseBNN.meshes[0],lDebug=lDebug)
    # aletBNN.save("./VIF_exp_DOSAVLET_fluscw_DEUTERON.bnn",lDebug=lDebug)

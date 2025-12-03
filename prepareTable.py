import numpy as np

eleNames=[
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr"
]

class Detector():
    def __init__(self):
        self.A=[]
        self.Z=[]
        self.m=[]
        self.Act=[]
        self.Err=[]
        self.ID=0
        self.name=""
    def fromTabLis_Isotope(self,tmpLine=None,tmpData=None):
        if (tmpLine is not None):
            dataIN=tmpLine.strip().split()
        elif (tmpData is not None):
            dataIN=tmpData
        else:
            print("No data to parse when acquiring detector data!")
            exit(1)
        if (float(dataIN[2])!=0.0):
            self.A.append(int(float(dataIN[0])+1.0E-14))
            self.Z.append(int(float(dataIN[1])+1.0E-14))
            self.Act.append(float(dataIN[2]))
            self.Err.append(float(dataIN[3]))
            self.m.append(0)
    def fromTabLis_Isomer(self,tmpLine=None,tmpData=None):
        if (tmpLine is not None):
            dataIN=tmpLine.strip().split()
        elif (tmpData is not None):
            dataIN=tmpData
        else:
            print("No data to parse when acquiring detector data!")
            exit(1)
        if (float(dataIN[3])!=0.0):
            self.A.append(int(float(dataIN[0])+1.0E-14))
            self.Z.append(int(float(dataIN[1])+1.0E-14))
            self.m.append(int(float(dataIN[2])+1.0E-14))
            self.Act.append(float(dataIN[3]))
            self.Err.append(float(dataIN[4]))
    def setID(self,tmpLine):
        data=tmpLine.split()
        self.ID=int(float(data[-1])+1.0E-14)
    def setName(self,tmpLine):
        data=tmpLine.split()
        self.name=data[0]
    def echoSummary(self):
        return "--> detector name: %10s - ID: % 3d; - % 4d isotopes, % 4d isomers;"%(
            self.name,self.ID,
            sum([ tmpM==0 for tmpM in self.m ]),
            sum([ tmpM!=0 for tmpM in self.m ]))
    def echoIsotopeFmt(self,ii):
        return "%d,%s,%d,%g,%g"%(
            self.Z[ii],eleNames[self.Z[ii]-1],self.A[ii],self.Act[ii],self.Err[ii])
    def echoIsomerFmt(self,ii):
        return "%d,%s,%d,%d,%g,%g"%(
            self.Z[ii],eleNames[self.Z[ii]-1],self.A[ii],self.m[ii],self.Act[ii],self.Err[ii])
    def saveCsv(self,oFileName,lSeparate=True):
        print("saving detector named %10s (ID=% 3d) in file %s ..."%(self.name,self.ID,oFileName))
        oFile=open(oFileName,"w")
        if (lSeparate):
            # separate isomers from isotopes
            if (sum([ tmpM is None for tmpM in self.m ])>0):
                print("...saving isotopes...")
                for ii in range(len(self.m)-1,-1,-1):
                    if (self.m[ii]==0):
                        oFile.write("%s\n"%(self.echoIsotopeFmt(ii)))
            if (sum([ tmpM is not None for tmpM in self.m ])>0):
                print("...saving isomers...")
                for ii in range(len(self.m)-1,-1,-1):
                    if (self.m[ii]!=0):
                        oFile.write("%s\n"%(self.echoIsomerFmt(ii)))
        else:
            print("...saving both isomers and isotopes mixed...")
            for ii in range(len(self.m)-1,-1,-1):
                oFile.write("%s\n"%(self.echoIsomerFmt(ii)))
        oFile.close()
        print("...done.")

def ParseTabLis(iFileName):
    detectors=[]
    print("parsing file %s ..."%(iFileName))
    iFile=open(iFileName,"r")
    iRead=-1
    for tmpLine in iFile.readlines():
        myLine=tmpLine.strip()
        if (myLine.startswith("# Detector")):
            detectors.append(Detector())
            iRead=0
            detectors[-1].setID(myLine)
            continue
        elif (myLine.startswith("# A/Z Isotopes")):
            iRead=1
            continue
        elif (myLine.startswith("# A/Z/m Isomers")):
            iRead=2
            continue
              
        if (iRead==0):
            detectors[-1].setName(myLine)
        elif (iRead==1):
            detectors[-1].fromTabLis_Isotope(tmpLine=myLine)
        elif (iRead==2):
            detectors[-1].fromTabLis_Isomer(tmpLine=myLine)
    iFile.close()
    print("...done: acquired %d detectors."%(len(detectors)))
    for ii in range(len(detectors)):
        print(detectors[ii].echoSummary())
            
    return detectors

if (__name__=="__main__"):
    myCase="C_Polyethy"
    iFileName="./%s/FLUKAbench_47_tab.lis"%(myCase)
    oFileName="./%s_7d.csv"%(myCase)
    tCools=[60.,86400.,259200.,604800.,1296000.,2592000.]
    detectors=ParseTabLis(iFileName)
    detectors[3].saveCsv(oFileName,lSeparate=False)

from ROOT import *

class JetCorrectionUncertainty:
    def __init__(self, dataFile):
        JESFiles = open(dataFile,"r")
        self.angleMinJES = []
        self.angleMaxJES = []
        self.ptMinJES = []
        self.ptMaxJES = []
        self.jesSFUp = []
        self.jesSFDown = []
        lineNumber = 0
        for line in JESFiles:
            self.jesSFUp.append([])
            self.jesSFDown.append([])
            partialLine = line.split(" ")
            self.angleMinJES.append(partialLine[0])
            self.angleMaxJES.append(partialLine[1])
            for i in range(1,len(partialLine)/3):
                ind = i*3-1
                self.ptMinJES.append(partialLine[ind])
                self.ptMaxJES.append(partialLine[ind + 3])
                self.jesSFUp[lineNumber].append(partialLine[ind+1])
                self.jesSFDown[lineNumber].append(partialLine[ind+2])
            lineNumber += 1

    def getUncertainty(self,pt,eta,jesUD):
        if jesUD == 0: return 1.
        ptBin = 0
        etaBin = 0
        for i in range(len(self.ptMinJES)):
            if pt > self.ptMinJES[i] and pt <= self.ptMaxJES[i]:
                ptBin = i
                break
        for i in range(len(self.angleMinJES)):
            if eta > self.angleMinJES[i] and eta <= self.angleMaxJES[i]:
                etaBin = i
                break
        lowFact = 0.
        highFact = 0.
        if jesUD == 1:
            lowFact=float(self.jesSFUp[etaBin][ptBin])
            highFact = float(self.jesSFUp[etaBin][ptBin+1])
        else:
            lowFact=float(self.jesSFDown[etaBin][ptBin])
            highFact = float(self.jesSFDown[etaBin][ptBin+1])
        a = (highFact - lowFact)/(float(self.ptMaxJES[ptBin])-float(self.ptMinJES[ptBin]))
        b = (lowFact*float(self.ptMaxJES[ptBin]) - highFact*float(self.ptMinJES[ptBin]))/(float(self.ptMaxJES[ptBin])-float(self.ptMinJES[ptBin]))
        return a*pt + b


    def getMetAfterJESUnc(self,metPx,metPy,tree,jesUD):
        for i in range(tree.numJetPF2PAT):
            metPx += tree.jetPF2PATPx[i]
            metPy += tree.jetPF2PATPy[i]
            uncertainty = self.getUncertainty(tree.jetPF2PATPt[i],tree.jetPF2PATEta[i],jesUD)
            if jesUD == 1:
                metPx -= (1 + uncertainty)*tree.jetPF2PATPx[i]
                metPy -= (1 + uncertainty)*tree.jetPF2PATPy[i]
            else:
                metPx -= (1 - uncertainty)*tree.jetPF2PATPx[i]
                metPy -= (1 - uncertainty)*tree.jetPF2PATPy[i]
        return (metPx,metPy)



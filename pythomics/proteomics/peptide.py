import pythomics.proteomics.config as config

class Peptide(object):
    def __init__(self, sequence):
        self.sequence = sequence
        
    def getCharge(self, pH=7.0):
        d = config.PI_RESIDUE_CHARGE
        pH = float(pH)
        s = self.sequence.upper()
        l = [(i, float(s.count(i))) for i in set(list(s))]
        #do n-term/c-term first
        charge = (1.0/(1.0+10.0**(pH-9.69)))
        charge += (-1.0/(1.0+10.0**(-1.0*(pH-2.34))))
        for aa, count in l: 
            v = d.get(aa, (0.0,0.0))
            charge += (v[1]*count)/(1.0+10.0**(v[1]*(pH-v[0])))
        return charge
        
    def getPI(self):
        pH = 6.5
        left_bound = 0.0
        right_bound = 14.0
        epsilon = 1
        charge = self.getCharge(pH=pH)
        while epsilon > 0.01:
            if charge < 0:
                #drop it by half
                new_pH = (pH-left_bound)/2.0
                right_bound = pH
            else:
                #increase it by half
                new_pH = pH+(right_bound-pH)/2.0
                left_bound = pH
            epsilon = abs(pH-new_pH)
            charge = self.getCharge(pH=new_pH)
            pH = new_pH
        return pH
        
    def getMass(self, average=False):
        if average:
            pass
        else:
            return sum([config.RESIDUE_MASSES[i][0] for i in self.sequence.upper()])+config.MODIFICATION_MASSES['h2o'][0]


    def getMSMass(self, charge=0):
        return (self.getMass()+charge*config.PROTON)/float(charge)
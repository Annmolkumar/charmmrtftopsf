#!/usr/bin/python
######################################################################
## Authors: Anmol Kumar
## created :  01/09/2019                                             #
## latest modified: 01/09/2019                                       #
## E-mail: anmol@outerbanks.umaryland.edu                            #
######################################################################

import os, logging, sys
import numpy as np
import logging


#logging.basicConfig(filename=sourcepath+'/'+'crd_top.out',level=print)

class Topology():
    '''
    Topology class that parses and returns the information about a residue or a patch.
    Requires two arguments to initialize, name of the topology file and the residue name.
    Third argument is optional,if not given entire information about the molecule is returned by default.
    Available parameters to be returned can be a list or a single parameter.
    Run example:
    topo=Topology('name.str','resname',ff='additive',returnvars='all')
    print topo.params
    '''

    def __init__(self,strfile,resi,ff=None,returnvars='all'):
        self.strfile=strfile
        self.resi=resi
        if len(self.resi) < 3:
            logging.error('Residue name is smaller than 3, are you sure there is no typo ?')
        self.returnvars=returnvars
        # This class checks for the ff type if the user did not supplied it
        if ff is not None:
            # If it is not passed readstr populates the self.ff value, not ideal but couldn't think of another way for it to work
            self.ff=ff.lower()
        if ff not in ['additive','drude',None]:
            logging.error('Only additive or drude forcefields (Charmm36) are accepted. Is this a typo: %s' % ff)
        self.param=self.readstr(self.returnvars,ff)
        if self.param is None:
           self.control=False
        else:
           self.control=True
           self.mapnum=self.createmap()
           self.mapsym=self.createmap(reverse=True)
           bondlistsym,bondlistnum = self.findbonds()
           anglelistsym,anglelistnum = self.findangles()
           dihelistsym,dihelistnum = self.finddihedrals()
           self.icnum = {"bonds":bondlistnum,"angles":anglelistnum,"dihedrals":dihelistnum}
           self.icsym = {"bonds":bondlistsym,"angles":anglelistsym,"dihedrals":dihelistsym}
#        for atom in self.param['atomnames']:
#            setattr(self, atom, self.Atominfo(atom,self.param))
        self.writetop()

    def readstr(self,returnvars,ff=None):
        '''
        Reads a topology and returns the parameters as a dictionary.
        Requires the parameter type to return. 'all' returns everything.
        '''
        anamnolp={}
        anamall={}
        atyp={} #[]
        achrg={} #[]
        aalp={}
        athl={}
        bndlist={}
        
        try:
            filein=open(self.strfile, "r")
        except IOError:
            logging.error('I cannot access the topology file')
            return #sys.exit()
        molinfo=[]
        found=False
        flag='wrng'
        for line in filein:
            #Strip the line, remove the comment lines, and split using whitespace
            field = line.strip().split("!")[0].split()
            if len(field) > 0:
                if field[0].upper() in ["RESI","PRES"] and field[1].upper() !=  self.resi.upper():
                    if flag is "crct":
                        break
                    else:
                        flag = "wrng"
                if field[0].upper() in ["RESI","PRES"] and field[1].upper() ==  self.resi.upper():
                    found=True
                    #reschrg = float(field[2])
                    flag = "crct"
                    molinfo.append(field)
                if flag == "crct":
                    molinfo.append(field)
        if found is False:
            logging.error('I could not find the residue name you have given.')
            return None
        reschrg=float(molinfo[0][2])
        ring=1
        ringatoms={}
        j = 1
        for field in molinfo:
            if "RING" in field[0]:
                ringatoms[ring]=field[1:]
                ring=ring+1
                continue
            if len(field) > 3 and field[0] == "ATOM":
                anamall[field[1]] = j
                #inamall[j] = field[1]
                j = j + 1
                atyp[field[1]]=field[2]
                achrg[field[1]]=float(field[3])
                try:
                    field[4] == "ALPHA"
                    aalp[field[1]]=float(field[5])
                except IndexError:
                    pass
                try:
                    field[6] == "THOLE"
                    athl[field[1]]=float(field[7])
                except IndexError:
                    pass
            #Bond list dictionary of each atom is populated.
            if field[0] in ["BOND","DOUB"] and (len(field)-1) % 2 == 0:
                for i in range(1,len(field),2):
                    if 'LP' not in [field[i][0:2], field[i+1][0:2]]:
                        try:
                            bndlist[field[i]].append(field[i+1])
                        except KeyError:
                            bndlist[field[i]]=[field[i+1]]
                        try:
                            bndlist[field[i+1]].append(field[i])
                        except KeyError:
                            bndlist[field[i+1]]=[field[i]]
            elif field[0] in ["BOND","DOUB"] and (len(field)-1) % 2 != 0:
                logging.error('Topology file BOND line has the wrong format. It should be multiples of two! BOND At1 At2 At2 At3 At1 At3 etc.')
                return None #sys.exit()
        # This is the check for drude and additive if user didn't supply any value

        if ff is None:
           #Checks only the length of alpha values since this is always be present if it is a drude
           self.ff='drude' if len(aalp) > 0 else 'additive'
           #Now saw the need for this. Do not remove, we use the ff object below
           ff=self.ff
        #Creates the atom names and lists. It preserves the order from tpr which can be used as numbers
        natomall=len(anamall)
        i = 1
        for atom in anamall.keys(): 
            if atom[0:2] != "LP":
               anamnolp[atom] = i  
               #inamnolp[i] = atom  
               i = i + 1
        natom=len(anamnolp)
        
        if ff == 'drude':
            i = 1
            anamdlp={}
            for atom in anamall:
                anamdlp[atom] = i
                i = i + 1
                if atom[0:1] != "H" and atom[0:2] != "LP":
                    anamdlp["D"+str(atom)] = i
                    i = i + 1
            natomdlp=len(anamdlp)
        else:
            #The only other case is additive
            #anamnolp=anamnall #Passing only the nolp names here. Does this break anything for halogens.
            natom=len(anamnolp)
            natomdlp=anamdlp=None
        #Checks if there are any atoms not bonded to anything or any bonds with no defined atoms
        for atom in anamnolp.keys():
            #Issue 1, it should have failed here. But since this only logs and doesn't exist, it progresses further. We need to catch these errors and print them.
            if atom not in list(bndlist.keys()) and atom[0:2] != "LP":
                logging.error(('Atom %s is in atom list but no connectivity is found.' %(atom)))
        for atom in list(bndlist.keys()):
            if atom not in anamnolp.keys():
                #Make this check for LPs too.
                logging.error(('A bond is defined for atom %s but it is not in atom list.' %(atom)))
        #Check for residue total charge
        if abs(sum(achrg.values())-reschrg) > 0.00001:
            logging.error('Warning: Your total charge does not match the sum of partial charges.')

        #Returns none for dlp lp parameters if it is additive
        toreturn={'atomnames':anamnolp,'atomlpnames':anamall,'atomdlpnames':anamdlp,'types':atyp,'parcharges':achrg,'alphas':aalp,'tholes':athl,'bonds':bndlist,'rescharge':reschrg,'natom':natom,'natomlp':natomall,'natomdlp':natomdlp,'totalcharge':reschrg,'resnam':self.resi, 'ringinfo':[ringatoms]}

        if returnvars is 'all':
            return(toreturn)
        elif type(returnvars) is str and returnvars in toreturn:
            return({returnvars:toreturn[returnvars]})
        elif type(returnvars) is list:
            returnjoint={}
            for var in returnvars:
                if var in toreturn:
                    returnjoint[var]=toreturn[var]
                else:
                    logging.error('Parameter type %s not found.' %s(var))
                    #sys.exit()
            return(returnjoint)
        else:
            logging.error('Parameter type(s) not understood or not found.')

    def createmap(self,reverse=False):
        mapping={}
        i = 0
        if self.ff == 'additive':
           for key in self.param['atomnames'].keys():
               i = i + 1
               mapping[i]=[key,self.param['atomlpnames'][key],self.param['atomlpnames'][key]]
        elif self.ff == 'drude':
           for key in self.param['atomnames'].keys():
              i = i + 1
              mapping[i]=[key,self.param['atomlpnames'][key],self.param['atomdlpnames'][key]]

        if reverse is True:
            newmap={}
            for atom in mapping:
                newmap[mapping[atom][0]]=[atom]+mapping[atom][1:]
            mapping=newmap
            mapping['resname']=self.resi
        return(mapping)

    def findbonds(self):
        listofbonds = self.param['bonds']
        atomlist = self.param['atomnames']
        self.bonds = []
        self.bondnum = []
        for key,value in list(listofbonds.items()):
            for val in value:
                if [val,key] not in self.bonds:
                    self.bonds.append([key,val])
                    self.bondnum.append([atomlist[key],atomlist[val]])
                    #Issue1 Requires an exception so that if an atom is not defined but a bond exists it fails. But it should have failed above already!!!
        return(self.bonds,self.bondnum)
#            for key2 in self.mapping.keys():
#                if self.mapping[key2][0] == key:
#                   ia = key2
#                   break
#            for key2 in self.mapping.keys():
#                for bondat in listofbonds[key]:
#                    if self.mapping[key2][0] == bondat:
#                       ib = key2
#                       if [self.mapping[ib][0],self.mapping[ia][0]] not in self.bonds:
#                          self.bonds.append([self.mapping[ia][0],self.mapping[ib][0]])
#                          pass

    def findangles(self):
        listofbonds = self.param['bonds']
        atomlist = self.param['atomnames']
        self.angles = []
        self.anglenum = []
        for k0,v0 in list(listofbonds.items()):
            if len(v0) > 1:
               for i in range(0,(len(v0)-1)):
                   ang2 = k0
                   ang1 = v0[i]
                   for j in range(i+1,len(v0)):
                       ang3 = v0[j]
                       self.angles.append([ang1,ang2,ang3])
                       #self.anglenum.append([atomlist.index(ang1)+1,atomlist.index(ang2)+1,atomlist.index(ang3)+1])
                       self.anglenum.append([atomlist[ang1],atomlist[ang2],atomlist[ang3]])
        return(self.angles,self.anglenum)

    def finddihedrals(self):
        listofbonds = self.param['bonds']
        atomlist = self.param['atomnames']
        self.dihedrals = []
        self.dihedralnum = []
        for k0,v0 in list(listofbonds.items()):
            dih0 = k0
            for k1 in v0:
                if len(listofbonds[k1]) != 1:
                   dih1 = k1
                   for k2 in listofbonds[k1]:
                         if len(listofbonds[k2]) != 1 and k2 != dih0:
                             dih2 = k2
                             for k3 in listofbonds[k2]:
                                 #if len(listofbonds[k3]) != 1 and k3 != dih1 and k3 != dih0:
                                 if k3 != dih1 and k3 != dih0:
                                    dih3 = k3
                                    if [dih3,dih2,dih1,dih0] not in self.dihedrals:
                                       self.dihedrals.append([dih0,dih1,dih2,dih3])
                                       #self.dihedralnum.append([atomlist.index(dih0)+1,atomlist.index(dih1)+1,atomlist.index(dih2)+1,atomlist.index(dih3)+1])
                                       self.dihedralnum.append([atomlist[dih0],atomlist[dih1],atomlist[dih2],atomlist[dih3]])
            
        return(self.dihedrals,self.dihedralnum)

    def writetop(self):
        import datetime
        from masses import masses
        from itertools import cycle 
        wdict = masses()
        now = datetime.datetime.now()
        fout = open("try.psf","w")
        fout.write("PSF EXT CMAP CHEQ XPLOR\n")
        fout.write("\n")
        fout.write('{:>10} !NTITLE\n'.format(2))
        fout.write("* Created by FFParam\n")
        fout.write("* DATE:     "+now.strftime("%Y-%m-%d %H:%M")+"      CREATED BY: Dr. Anmol Kumar\n")
        fout.write("\n")
        fout.write('{:>10} !NATOM\n'.format(self.param["natom"]))
        names = []
        types = []
        charges = []
        weights = []
         
        for key,value in list(self.param["atomnames"].items()):
            names.append(key)
            types.append(self.param["types"][key])
            charges.append(self.param["parcharges"][key])
            weights.append(wdict[self.param["types"][key]])
        
        for i in range(int(self.param["natom"])):
            fout.write("{:>10} {:^}     1        {:^}  {:<}    {:<}     {:>10.6f}    {:>10.4f}    0  0.00000   -0.301140E-02\n".format(i+1,self.param["resnam"],self.param["resnam"],names[i],types[i],charges[i],weights[i]))
        fout.write("\n")

        fout.write('{:>10} !NBOND: bonds\n'.format(len(self.icnum["bonds"])))
        self.writeic(fout,"bonds")
        fout.write('{:>10} !NANGLE: angles\n'.format(len(self.icnum["angles"])))
        self.writeic(fout,"angles")
        fout.write('{:>10} !NDIHEDRAL: dihedrals\n'.format(len(self.icnum["dihedrals"])))
        self.writeic(fout,"dihedrals")
        fout.close()
        #fout.write('{:>10} !NIMPHI: impropers\n'.format(len(self.icnum["impropers"])))
        # 5         1        12        13

        #fout.write('{:>10} !NIMPHI: impropers\n'.format(len(self.icnum["impropers"])))
        # 0 !NDON: donors

        #fout.write('{:>10} !NIMPHI: impropers\n'.format(len(self.icnum["impropers"])))

        # 0 !NACC: acceptors


        #fout.write('{:>10} !NIMPHI: impropers\n'.format(len(self.icnum["impropers"])))
        # 0 !NNB

        # 0         0         0         0         0         0         0         0
        # 0         0         0         0         0

        # 1         0 !NGRP NST2
        # 0         1         0

        # 1 !MOLNT
        # 1         1         1         1         1         1         1         1
        # 1         1         1         1         1

        # 0         0 !NUMLP NUMLPH


    def writeic(self,fout,ictype):
        toprint =""
        counter = 0
        if ictype == "bonds":batch = 4
        if ictype == "angles":batch = 3
        if ictype == "dihedrals":batch = 2
        for i in range(0,len(self.icnum[ictype])):
            counter = counter + 1
            if counter%batch == 0:  
               topre   = "      ".join(map(str,self.icnum[ictype][i]))
               toprint = "      ".join((toprint,topre))
               toprint = toprint+"\n"
            else:
               topre   = "      ".join(map(str,self.icnum[ictype][i]))
               toprint = "      ".join((toprint,topre))
        fout.write("          {:}\n".format(toprint))
        fout.write("\n")

import argparse
def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-s',"--str", help='an integer for the accumulator',required=True)
    parser.add_argument('-r',"--resi", help='an integer for the accumulator',required=True)
    parser.add_argument('-ff',"--fftype",help='additive or drude')
    args = parser.parse_args()
    Topology(args.str,args.resi,args.fftype)

if __name__=="__main__":
   main()

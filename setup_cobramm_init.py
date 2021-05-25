#!/usr/bin/env python2

#******************************************
#
#    SHARC Program Suite
#
#    Copyright (c) 2019 University of Vienna
#
#    This file is part of SHARC.
#
#    SHARC is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    SHARC is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    inside the SHARC manual.  If not, see <http://www.gnu.org/licenses/>.
#
#******************************************

#!/usr/bin/env python2

# Interactive script for the setup of initial condition excitation calculations for SHARC
#
# usage: python setup_init.py

import copy
import math
import sys
import re
import os
import stat
import shutil
import datetime
from optparse import OptionParser
import readline
import time
import ast
import pprint

# =========================================================
# compatibility stuff

if sys.version_info[0]!=2:
  print 'This is a script for Python 2!'
  sys.exit(0)

if sys.version_info[1]<5:
  def any(iterable):
    for element in iterable:
      if element:
        return True
    return False

  def all(iterable):
    for element in iterable:
      if not element:
        return False
    return True

# some constants
DEBUG = False
CM_TO_HARTREE = 1./219474.6     #4.556335252e-6 # conversion factor from cm-1 to Hartree
HARTREE_TO_EV = 27.211396132    # conversion factor from Hartree to eV
U_TO_AMU = 1./5.4857990943e-4            # conversion from g/mol to amu
BOHR_TO_ANG=0.529177211
PI = math.pi

version='2.1'
versionneeded=[0.2, 1.0, 2.0, float(version)]
versiondate=datetime.date(2019,9,1)


IToMult={
         1: 'Singlet', 
         2: 'Doublet', 
         3: 'Triplet', 
         4: 'Quartet', 
         5: 'Quintet', 
         6: 'Sextet', 
         7: 'Septet', 
         8: 'Octet', 
         'Singlet': 1, 
         'Doublet': 2, 
         'Triplet': 3, 
         'Quartet': 4, 
         'Quintet': 5, 
         'Sextet': 6, 
         'Septet': 7, 
         'Octet': 8
         }

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================


def try_read(l,index,typefunc,default):
  try:
    return typefunc(l[index])
  except IndexError:
    return typefunc(default)
  except ValueError:
    print 'Could not initialize object!'
    quit(1)

# ======================================================================================================================

class ATOM:
  def __init__(self,symb='??',num=0.,coord=[0.,0.,0.],m=0.,veloc=[0.,0.,0.]):
    self.symb  = symb
    self.num   = num
    self.coord = coord
    self.mass  = m
    self.veloc = veloc
    self.Ekin=0.5*self.mass * sum( [ self.veloc[i]**2 for i in range(3) ] )

  def init_from_str(self,initstring=''):
    f=initstring.split()
    self.symb  =   try_read(f,0,str,  '??')
    self.num   =   try_read(f,1,float,0.)
    self.coord = [ try_read(f,i,float,0.) for i in range(2,5) ]
    self.mass  =   try_read(f,5,float,0.)*U_TO_AMU
    self.veloc = [ try_read(f,i,float,0.) for i in range(6,9) ]
    self.Ekin=0.5*self.mass * sum( [ self.veloc[i]**2 for i in range(3) ] )

  def __str__(self):
    s ='%2s % 5.1f '               % (self.symb, self.num)
    s+='% 12.8f % 12.8f % 12.8f '  % tuple(self.coord)
    s+='% 12.8f '                  % (self.mass/U_TO_AMU)
    s+='% 12.8f % 12.8f % 12.8f'   % tuple(self.veloc)
    return s

  def EKIN(self):
    self.Ekin=0.5*self.mass * sum( [ self.veloc[i]**2 for i in range(3) ] )
    return self.Ekin

  def geomstring(self):
    s='  %2s % 5.1f % 12.8f % 12.8f % 12.8f % 12.8f' % (self.symb,self.num,self.coord[0],self.coord[1],self.coord[2],self.mass/U_TO_AMU)
    return s

  def velocstring(self):
    s=' '*11+'% 12.8f % 12.8f % 12.8f' % tuple(self.veloc)
    return s

# ======================================================================================================================

class STATE:
  def __init__(self,i=0,e=0.,eref=0.,dip=[0.,0.,0.]):
    self.i       = i
    self.e       = e.real
    self.eref    = eref.real
    self.dip     = dip
    self.Excited = False
    self.Eexc    = self.e-self.eref
    self.Fosc    = (2./3.*self.Eexc*sum( [i*i.conjugate() for i in self.dip] ) ).real
    if self.Eexc==0.:
      self.Prob  = 0.
    else:
      self.Prob  = self.Fosc/self.Eexc**2

  def init_from_str(self,initstring):
    f=initstring.split()
    self.i       =   try_read(f,0,int,  0 )
    self.e       =   try_read(f,1,float,0.)
    self.eref    =   try_read(f,2,float,0.)
    self.dip     = [ try_read(f,i,float,0.) for i in range(3,6) ]
    self.Excited =   try_read(f,2,bool, False)
    self.Eexc    = self.e-self.eref
    self.Fosc    = (2./3.*self.Eexc*sum( [i*i.conjugate() for i in self.dip] ) ).real
    if self.Eexc==0.:
      self.Prob  = 0.
    else:
      self.Prob  = self.Fosc/self.Eexc**2

  def __str__(self):
    s ='%03i % 18.10f % 18.10f ' % (self.i,self.e,self.eref)
    for i in range(3):
      s+='% 12.8f % 12.8f ' % (self.dip[i].real,self.dip[i].imag)
    s+='% 12.8f % 12.8f %s' % (self.Eexc*HARTREE_TO_EV,self.Fosc,self.excited)
    return s

  def Excite(self,max_Prob,erange):
    try:
      Prob=self.Prob/max_Prob
    except ZeroDivisionError:
      Prob=-1.
    if not (erange[0] <= self.Eexc <= erange[1]):
      Prob=-1.
    self.excited=(random.random() < Prob)

# ======================================================================================================================

class INITCOND:
  def __init__(self,atomlist=[],eref=0.,epot_harm=0.):
    self.atomlist=atomlist
    self.eref=eref
    self.Epot_harm=epot_harm
    self.natom=len(atomlist)
    self.Ekin=sum( [atom.Ekin for atom in self.atomlist] )
    self.statelist=[]
    self.nstate=0
    self.Epot=epot_harm

  def addstates(self,statelist):
    self.statelist=statelist
    self.nstate=len(statelist)
    self.Epot=self.statelist[0].e-self.eref

  def init_from_file(self,f,eref,index):
    while True:
      line=f.readline()
      #if 'Index     %i' % (index) in line:
      if re.search('Index\s+%i' % (index),line):
        break
      if line=='\n':
        continue
      if line=='':
        print 'Initial condition %i not found in file %s' % (index,f.name)
        quit(1)
    f.readline()        # skip one line, where "Atoms" stands
    atomlist=[]
    while True:
      line=f.readline()
      if 'States' in line:
        break
      atom=ATOM()
      atom.init_from_str(line)
      atomlist.append(atom)
    statelist=[]
    while True:
      line=f.readline()
      if 'Ekin' in line:
        break
      state=STATE()
      state.init_from_str(line)
      statelist.append(state)
    epot_harm=0.
    while not line=='\n' and not line=='':
      line=f.readline()
      if 'epot_harm' in line.lower():
        epot_harm=float(line.split()[1])
        break
    self.atomlist=atomlist
    self.eref=eref
    self.Epot_harm=epot_harm
    self.natom=len(atomlist)
    self.Ekin=sum( [atom.Ekin for atom in self.atomlist] )
    self.statelist=statelist
    self.nstate=len(statelist)
    if self.nstate>0:
      self.Epot=self.statelist[0].e-self.eref
    else:
      self.Epot=epot_harm

  def __str__(self):
    s='Atoms\n'
    for atom in self.atomlist:
      s+=str(atom)+'\n'
    s+='States\n'
    for state in self.statelist:
      s+=str(state)+'\n'
    s+='Ekin      % 16.12f a.u.\n' % (self.Ekin)
    s+='Epot_harm % 16.12f a.u.\n' % (self.Epot_harm)
    s+='Epot      % 16.12f a.u.\n' % (self.Epot)
    s+='Etot_harm % 16.12f a.u.\n' % (self.Epot_harm+self.Ekin)
    s+='Etot      % 16.12f a.u.\n' % (self.Epot+self.Ekin)
    s+='\n\n'
    return s

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def check_initcond_version(string,must_be_excited=False):
  if not 'sharc initial conditions file' in string.lower():
    return False
  f=string.split()
  for i,field in enumerate(f):
    if 'version' in field.lower():
      try:
        v=float(f[i+1])
        if not v in versionneeded:
          return False
      except IndexError:
        return False
  if must_be_excited:
    if not 'excited' in string.lower():
      return False
  return True


# ======================================================================================================================

def centerstring(string,n,pad=' '):
  l=len(string)
  if l>=n:
    return string
  else:
    return  pad*((n-l+1)/2)+string+pad*((n-l)/2)

# ======================================================================================================================

def displaywelcome():
  print 'Script for setup of initial conditions started...\n'
  string='\n'
  string+='  '+'='*80+'\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Setup initial conditions for QM/MM SHARC/COBRAMM dynamics',80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Author: Davide Avagliano',80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Version:'+version,80)+'||\n'
  string+='||'+centerstring(versiondate.strftime("%d.%m.%y"),80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  string+='''
This script automatizes the setup of excited-state calculations for initial conditions
for SHARC dynamics.
  '''
  print string

# ======================================================================================================================

def open_keystrokes():
  global KEYSTROKES
  KEYSTROKES=open('KEYSTROKES.tmp','w')

def close_keystrokes():
  KEYSTROKES.close()
  shutil.move('KEYSTROKES.tmp','KEYSTROKES.setup_cobramm_init')

# ===================================

def question(question,typefunc,default=None,autocomplete=True,ranges=False):
  if typefunc==int or typefunc==float:
    if not default==None and not isinstance(default,list):
      print 'Default to int or float question must be list!'
      quit(1)
  if typefunc==str and autocomplete:
    readline.set_completer_delims(' \t\n;')
    readline.parse_and_bind("tab: complete")    # activate autocomplete
  else:
    readline.parse_and_bind("tab: ")            # deactivate autocomplete

  while True:
    s=question
    if default!=None:
      if typefunc==bool or typefunc==str:
        s+= ' [%s]' % (str(default))
      elif typefunc==int or typefunc==float:
        s+= ' ['
        for i in default:
          s+=str(i)+' '
        s=s[:-1]+']'
    if typefunc==str and autocomplete:
      s+=' (autocomplete enabled)'
    if typefunc==int and ranges:
      s+=' (range comprehension enabled)'
    s+=' '

    line=raw_input(s)
    line=re.sub('#.*$','',line).strip()
    if not typefunc==str:
      line=line.lower()

    if line=='' or line=='\n':
      if default!=None:
        KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
        return default
      else:
        continue

    if typefunc==bool:
      posresponse=['y','yes','true', 't', 'ja',  'si','yea','yeah','aye','sure','definitely']
      negresponse=['n','no', 'false', 'f', 'nein', 'nope']
      if line in posresponse:
        KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
        return True
      elif line in negresponse:
        KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
        return False
      else:
        print 'I didn''t understand you.'
        continue

    if typefunc==str:
      KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
      return line

    if typefunc==float:
      # float will be returned as a list
      f=line.split()
      try:
        for i in range(len(f)):
          f[i]=typefunc(f[i])
        KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
        return f
      except ValueError:
        print 'Please enter floats!'
        continue

    if typefunc==int:
      # int will be returned as a list
      f=line.split()
      out=[]
      try:
        for i in f:
          if ranges and '~' in i:
            q=i.split('~')
            for j in range(int(q[0]),int(q[1])+1):
              out.append(j)
          else:
            out.append(int(i))
        KEYSTROKES.write(line+' '*(40-len(line))+' #'+s+'\n')
        return out
      except ValueError:
        if ranges:
          print 'Please enter integers or ranges of integers (e.g. "-3~-1  2  5~7")!'
        else:
          print 'Please enter integers!'
        continue



# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def get_general():
  '''This routine questions from the user some general information:
  - initconds file
  - number of states
  - number of initial conditions
  - interface to use'''
  INFOS={}

  print centerstring('Initial conditions file',60,'-')+'\n'
  # open the initconds file
  try:
    initfile='initconds'
    initf=open(initfile)
    line=initf.readline()
    if check_initcond_version(line):
      print 'Initial conditions file "initconds" detected. Do you want to use this?'
      if not question('Use file "initconds"?',bool,True):
        initf.close()
        raise IOError
    else:
      initf.close()
      raise IOError
  except IOError:
    print '\nIf you do not have an initial conditions file, prepare one with wigner.py!\n'
    print 'Please enter the filename of the initial conditions file.'
    while True:
      initfile=question('Initial conditions filename:',str,'initconds')
      initfile=os.path.expanduser(os.path.expandvars(initfile))
      if os.path.isdir(initfile):
        print 'Is a directory: %s' % (initfile)
        continue
      if not os.path.isfile(initfile):
        print 'File does not exist: %s' % (initfile)
        continue
      try:
        initf=open(initfile,'r')
      except IOError:
        print 'Could not open: %s' % (initfile)
        continue
      line=initf.readline()
      if check_initcond_version(line):
        break
      else:
        print 'File does not contain initial conditions!'
        continue
  # read the header
  ninit=int(initf.readline().split()[1])
  natom=int(initf.readline().split()[1])
  INFOS['ninit']=ninit
  INFOS['natom']=natom
  initf.seek(0)                 # rewind the initf file
  INFOS['initf']=initf
  print '\nFile "%s" contains %i initial conditions.' % (initfile,ninit)
  print 'Number of atoms is %i\n' % (natom)



  print centerstring('Range of initial conditions',60,'-')
  print '\nPlease enter the range of initial conditions for which an excited-state calculation should be performed as two integers separated by space.'
  while True:
    irange=question('Initial condition range:',int,[1,ninit])
    if len(irange)!=2:
      print 'Enter two numbers separated by spaces!'
      continue
    if irange[0]>irange[1]:
      print 'Range empty!'
      continue
    if irange[0]==irange[1]==0:
      print 'Only preparing calculation at equilibrium geometry!'
      break
    if irange[1]>ninit:
      print 'There are only %i initial conditions in file %s!' % (ninit,initfile)
      continue
    if irange[0]<=0:
      print 'Only positive indices allowed!'
      continue
    break
  print '\nScript will use initial conditions %i to %i (%i in total).\n' % (irange[0],irange[1],irange[1]-irange[0]+1)
  INFOS['irange']=irange

  
  print centerstring('Full system topology file',60,'-')
  print '\nPlease specify the real.top file (AMBER topology file):'
  while True:
    path2=question('File name:',str,'real.top')
    try:
      rf=open(path2,'r')
    except IOError:
      print 'Could not open: %s' % (path2)
      continue
    break
  INFOS['realtop_location']=path2
  #realtop_data=readfile(INFOS['realtop_location'])
  print '\nTopology is read in %s\n' % (path2)

  print centerstring('QM system topology file',60,'-')
  print '\nPlease specify the model-H.top file (AMBER topology file):'
  while True:
    path3=question('File name:', str, 'model-H.top')
    try:
      mf=open(path3,'r')
    except IOError:
      print 'Could not open: %s' % (path3)
      continue
    break
  INFOS['modeltop_location']=path3
  print '\nTopology is read in %s\n' % (path3)

  print centerstring('Layer definition file',60,'-')
  print '\nPlease specify the real_layers file:'
  while True:
    path4=question('File name:', str, 'real_layers.xyz')
    try:
      lf=open(path4,'r')
    except IOError:
      print 'Could not open: %s' % (path4)
      continue
    break
  INFOS['realayers_location']=path4
  print '\nLayer definition is read in %s\n' % (path4)

  string='\n  '+'='*80+'\n'
  string+='||'+centerstring('COBRAMM Interface setup',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  print string

  print centerstring('Path to COBRAMM',60,'-')+'\n'
  path=os.getenv('COBRAM_PATH')
  path=os.path.expanduser(os.path.expandvars(path))
  if path=='':
   path=None
  else:
    path='$COBRAM_PATH/'
  print '\nPlease specify path to COBRAMM directory (SHELL variables and ~      can be used, will be expanded when interface is started).\n'
  INFOS['cobramm']=question('Path to COBRAMM:', str,path)


  print centerstring('Path to AMBER',60,'-')+'\n'
  path=os.getenv('AMBERHOME')
  path=os.path.expanduser(os.path.expandvars(path))
  if path=='':
   path=None
  else:
    path='$AMBERHOME/'
  print '\nPlease specify path to AMBER directory (SHELL variables and ~ ca     n be used, will be expanded when interface is started).\n'
  INFOS['amber']=question('Path to AMBER:', str,path)


  print centerstring('Scratch directory',60,'-')+'\n'
  print 'Please specify an appropriate scratch directory\n. This will conta     in inside a subdirectory QM.\nThis will be used to temporally store the int     egrals. The scratch directory will be deleted after the calculation.\nRemem     ber that this script cannot check whether the path is valid, since you may      run the calculations on a different machine. The path will not be expanded      by this script.'
  INFOS['scratchdir']=question('Path to scratch directory:',str)
  print ''

  print centerstring('COBRAMM input template file',60,'-')+'\n'
  print '''Please specify the path to the COBRAMM.template file. This file      must contain the following setting:

jobtype <sp, optx ecc> 
interface <name of the interface for QM calculation>

'''

  if os.path.isfile('COBRAMM.template'):
    if checktemplate_COBRAMM('COBRAMM.template',INFOS):
     print 'Valid file "COBRAMM.template" detected.'
     usethisone=question('Use this template file?',bool,True)
     if usethisone:
       INFOS['cobramm.template']='COBRAMM.template'
  if not 'cobramm.template' in INFOS:
    while True:
      filename=question('Template filename:',str)
      if not os.path.isfile(filename):
        print 'File %s does not exist!' % (filename)
        continue
      if checktemplate_COBRAMM(filename,INFOS):
        break
    INFOS['cobramm.template']=filename
  print ''

  print centerstring('COBRAMM Ressource usage',60,'-')+'\n'
  print '''Please specify the amount of memory available to Cobramm.
'''
  INFOS['cobramm.mem']=abs(question('COBRAMM memory (MB):',int,[1000])[0])
  print '''Please specify the number of CPUs to be used by EACH trajectory.
'''
  INFOS['cobramm.ncpu']=abs(question('Number of CPUs:',int,[1])[0])
 
#   return INFOS

  print centerstring('Number of states',60,'-')
  print '\nPlease enter the number of states as a list of integers\ne.g. 3 0 3 for three singlets, zero doublets and three triplets.'
  while True:
    states=question('Number of states:',int)
    if len(states)==0:
      continue
    if any(i<0 for i in states):
      print 'Number of states must be positive!'
      continue
    break
  print ''
  nstates=0
  for mult,i in enumerate(states):
    nstates+=(mult+1)*i
  print 'Number of states: '+str(states)
  print 'Total number of states: %i\n' % (nstates)
  INFOS['states']=states
  INFOS['nstates']=nstates
  print centerstring('Choose the quantum chemistry interface',60,'-')
  print '\nPlease specify the quantum chemistry interface (enter any of the following numbers):'
  for i in Interfaces:
    print '%i\t%s' % (i, Interfaces[i]['description'])
  print ''
  while True:
    num=question('Interface number:',int)[0]
    if num in Interfaces:
      break
    else:
      print 'Please input one of the following: %s!' % ([i for i in Interfaces])
  INFOS['interface']=num

  INFOS['needed']=[]

  # Setup SOCs
  print '\n'+centerstring('Spin-orbit couplings (SOCs)',60,'-')+'\n'
  if len(states)>1:
    if 'soc' in Interfaces[num]['features']:
      print 'Do you want to compute spin-orbit couplings?\n'
      soc=question('Spin-Orbit calculation?',bool,True)
      if soc:
        print 'Will calculate spin-orbit matrix.'
    else:
      print 'Interface cannot provide SOCs: not calculating spin-orbit matrix.'
      soc=False
  else:
    print 'Only singlets specified: not calculating spin-orbit matrix.'
    soc=False
  print ''
  INFOS['soc']=soc
  if INFOS['soc']:
    INFOS['needed'].extend(Interfaces[num]['features']['soc'])


  # Setup Dyson spectra
  if 'dyson' in Interfaces[num]['features']:
    n=[0,0]
    for i,j in enumerate(states):
      n[i%2]+=j
    if n[0]>=1 and n[1]>=1:
      print '\n'+centerstring('Ionization probability by Dyson norms',60,'-')+'\n'
      print 'Do you want to compute Dyson norms between neutral and ionic states?'
      INFOS['ion']=question('Dyson norms?',bool,False)
      if INFOS['ion']:
        INFOS['needed'].extend(Interfaces[num]['features']['dyson'])


  # Setup initconds with reference overlap
  if 'overlap' in Interfaces[num]['features']:
    print '\n'+centerstring('Overlaps to reference states',60,'-')+'\n'
    print 'Do you want to compute the overlaps between the states at the equilibrium geometry and the states at the initial condition geometries?'
    INFOS['refov']=question('Reference overlaps?',bool,False)
    if INFOS['refov']:
      INFOS['needed'].extend(Interfaces[num]['features']['overlap'])


  # Setup theodore
  if 'theodore' in Interfaces[num]['features']:
    print '\n'+centerstring('TheoDORE wave function analysis',60,'-')+'\n'
    print 'Do you want to run TheoDORE to obtain one-electron descriptors for the electronic wave functions?'
    INFOS['theodore']=question('TheoDORE?',bool,False)
    if INFOS['theodore']:
      INFOS['needed'].extend(Interfaces[num]['features']['theodore'])


  INFOS['cwd']=os.getcwd()
  print ''
  return INFOS

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

Interfaces={
# 1: {'script':          'SHARC_MOLPRO.py',
#     'name':            'molpro',
#     'description':     'MOLPRO (only CASSCF)',
#     'get_routine':     'get_MOLPRO',
#     'prepare_routine': 'prepare_MOLPRO',
#     'features':        {'overlap': ['wfoverlap'],
#                         'dyson':   ['wfoverlap'],
#                         'nacdr':   ['wfoverlap'],
#                         'phases':  ['wfoverlap'],
#                         'soc':     []             },
#     'pysharc':          False
#    },
# 2: {'script':          'SHARC_COLUMBUS.py',
#     'name':            'columbus',
#     'description':     'COLUMBUS (CASSCF, RASSCF and MRCISD), using SEWARD integrals',
#     'get_routine':     'get_COLUMBUS',
#     'prepare_routine': 'prepare_COLUMBUS',
#     'features':        {'overlap': ['wfoverlap'],
#                         'dyson':   ['wfoverlap'],
#                         'phases':  ['wfoverlap'],
#                         'nacdr':   [],
#                         'soc':     []               },
#     'pysharc':          False
#    },
  1: {'script':          'SHARC_MOLCAS.py',
      'name':            'molcas',
      'description':     'MOLCAS (CASSCF, CASPT2, MS-CASPT2)',
      'get_routine':     'get_MOLCAS',
      'prepare_routine': 'prepare_MOLCAS',
      'features':        {'overlap': [],
                          'dyson':   ['wfoverlap'],
                          'dipolegrad':[],
                          'phases':  [],
                          'nacdr':   [],
                          'soc':     []             },
      'pysharc':          False
     },
  2: {'script':          'SHARC_RICC2.py',
      'name':            'ricc2',
      'description':     'TURBOMOLE (ricc2 with CC2 and ADC(2))',
      'get_routine':     'get_RICC2',
      'prepare_routine': 'prepare_RICC2',
      'features':        {'overlap': ['wfoverlap'],
                          'theodore':['theodore'],
                          'phases':  ['wfoverlap'],
                          'soc':     []                 },
      'pysharc':          False
     },
# 6: {'script':          'SHARC_GAUSSIAN.py',
#     'name':            'gaussian',
#     'description':     'GAUSSIAN (DFT, TD-DFT)',
#     'get_routine':     'get_GAUSSIAN',
#     'prepare_routine': 'prepare_GAUSSIAN',
#     'features':        {'overlap': ['wfoverlap'],
#                         'dyson':   ['wfoverlap'],
#                         'theodore':['theodore'],
#                         'phases':  ['wfoverlap']        },
#     'pysharc':          False
#    },
  3: {'script':          'SHARC_ORCA.py',
      'name':            'orca',
      'description':     'ORCA (DFT, TD-DFT, HF, CIS)',
      'get_routine':     'get_ORCA',
      'prepare_routine': 'prepare_ORCA',
      'features':        {'overlap': ['wfoverlap'],
                          'dyson':   ['wfoverlap'],
                          'theodore':['theodore'],
                          'phases':  ['wfoverlap'],
                          'soc':     []},
      'pysharc':          False
     },
  }


# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
def checktemplate_COBRAMM(filename,INFOS):
  necessary=['jobtype','interface']
  try:
    f=open(filename)
    data=f.readlines()
    f.close()
  except IOError:
    print 'Could not open template file %s' % (filename)
    return False
  valid=[]
  for i in necessary:
    for l in data:
      line=l.lower()
      if i in re.sub('#.*$','',line):
        valid.append(True)
        break
    else:
      valid.append(False)
  if not all(valid):
    print 'The template %s seems to be incomplete! It should contain: ' % (          filename) +str(necessary)
    return False
  return True


# ===================================================================================================
def prepare_COBRAMM(INFOS,iconddir):
  try:
    sh2cbm=open('%s/COBRAMM.resources' % (iconddir), 'w')
  except IOError:
    print 'IOError during prepare_COBRAMM, iconddir=%s' % (iconddir)
    quit(1)
  string='''cobrammdir %s
amberdir %s
scratchdir %s/%s
memory %i
ncpu %i
''' % (INFOS['cobramm'],
       INFOS['amber'],
       INFOS['scratchdir'],
       iconddir,
       INFOS['cobramm.mem'],
       INFOS['cobramm.ncpu'])

  sh2cbm.write(string)
  sh2cbm.close()

  # copy MOs and template
  cpfrom=INFOS['cobramm.template']
  cpto='%s/COBRAMM.template' % (iconddir)
  shutil.copy(cpfrom,cpto)

  cpfrom=INFOS['realtop_location']
  cpto='%s/real.top' % (iconddir)
  shutil.copy(cpfrom,cpto)

  cpfrom=INFOS['modeltop_location']
  cpto='%s/model-H.top' % (iconddir)
  shutil.copy(cpfrom,cpto)

  cpfrom=INFOS['realayers_location']
  cpto='%s/real_layers.xyz' % (iconddir)
  shutil.copy(cpfrom,cpto)
  return



# ===================================================================================================

def checktemplate_MOLPRO(filename):
  necessary=['cobramm','basis','closed','occ','nelec','roots']
  try:
    f=open(filename)
    data=f.readlines()
    f.close()
  except IOError:
    print 'Could not open template file %s' % (filename)
    return False
  i=0
  for l in data:
    if necessary[i] in l:
      i+=1
      if i+1==len(necessary):
        return True
  print 'The template %s seems to be incomplete! It should contain: ' % (filename) +str(necessary)
  return False

# =================================================

def get_MOLPRO(INFOS):
  '''This routine asks for all questions specific to MOLPRO:
  - path to molpro
  - scratch directory
  - MOLPRO.template
  - wf.init
  '''

  string='\n  '+'='*80+'\n'
  string+='||'+centerstring('MOLPRO Interface setup',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  print string

  print centerstring('Path to MOLPRO',60,'-')+'\n'
  path=os.getenv('MOLPRO')
  path=os.path.expanduser(os.path.expandvars(path))
  if not path=='':
    path='$MOLPRO/'
  else:
    path=None
  #if path!='':
    #print 'Environment variable $MOLPRO detected:\n$MOLPRO=%s\n' % (path)
    #if question('Do you want to use this MOLPRO installation?',bool,True):
      #INFOS['molpro']=path
  #if not 'molpro' in INFOS:
  print '\nPlease specify path to MOLPRO directory (SHELL variables and ~ can be used, will be expanded when interface is started).\n'
  INFOS['molpro']=question('Path to MOLPRO executable:',str,path)
  print ''



  print centerstring('MOLPRO input template file',60,'-')+'\n'
  print '''Please specify the path to the MOLPRO.template file. This file must be a valid MOLPRO input file for a CASSCF calculation. It should contain the following settings:
- memory settings
- Basis set (possibly also Douglas-Kroll settings etc.)
- CASSCF calculation with:
  * Number of frozen, closed and occupied orbitals
  * wf and state cards for the specification of the wavefunction
MOLPRO.template files can easily be created using molpro_input.py (Open a second shell if you need to create one now).

The MOLPRO interface will generate the remaining MOLPRO input automatically.
'''
  if os.path.isfile('MOLPRO.template'):
    if checktemplate_MOLPRO('MOLPRO.template'):
      print 'Valid file "MOLPRO.template" detected. '
      usethisone=question('Use this template file?',bool,True)
      if usethisone:
        INFOS['molpro.template']='MOLPRO.template'
  if not 'molpro.template' in INFOS:
    while True:
      filename=question('Template filename:',str)
      if not os.path.isfile(filename):
        print 'File %s does not exist!' % (filename)
        continue
      if checktemplate_MOLPRO(filename):
        break
    INFOS['molpro.template']=filename
  print ''


  print centerstring('Initial wavefunction: MO Guess',60,'-')+'\n'
  print '''Please specify the path to a MOLPRO wavefunction file containing suitable starting MOs for the CASSCF calculation. Please note that this script cannot check whether the wavefunction file and the Input template are consistent!

If you optimized your geometry with MOLPRO/CASSCF you can reuse the "wf" file from the optimization.
'''
  if question('Do you have an initial wavefunction file?',bool,True):
    while True:
      filename=question('Initial wavefunction file:',str,'wf.init')
      if os.path.isfile(filename):
        break
      else:
        print 'File not found!'
    INFOS['molpro.guess']=filename
  else:
    print 'WARNING: Remember that CASSCF calculations may run very long and/or yield wrong results without proper starting MOs.'
    time.sleep(2)
    INFOS['molpro.guess']=False


  print centerstring('MOLPRO Ressource usage',60,'-')+'\n'
  print '''Please specify the amount of memory available to MOLPRO (in MB). For calculations including moderately-sized CASSCF calculations and less than 150 basis functions, around 2000 MB should be sufficient.
'''
  INFOS['molpro.mem']=abs(question('MOLPRO memory:',int,[500])[0])
  print '''Please specify the number of CPUs to be used by EACH calculation.
'''



  if 'wfoverlap' in INFOS['needed']:
    print '\n'+centerstring('WFoverlap setup',60,'-')+'\n'
    INFOS['molpro.wfpath']=question('Path to wavefunction overlap executable:',str,'$SHARC/wfoverlap.x')
    # TODO: not asked for: numfrozcore, numocc

  return INFOS

# =================================================

def prepare_MOLPRO(INFOS,iconddir):
  # write MOLPRO.resources
  try:
    sh2pro=open('%s/MOLPRO.resources' % (iconddir), 'w')
  except IOError:
    print 'IOError during prepareMOLPRO, iconddir=%s' % (iconddir)
    quit(1)
  string='''molpro %s
scratchdir %s/%s/QM/

memory %i
ncpu %i
''' % (INFOS['molpro'],INFOS['scratchdir'],iconddir,INFOS['molpro.mem'],INFOS['cobramm.ncpu'])
  if 'wfoverlap' in INFOS['needed']:
    string+='wfoverlap %s\n' % (INFOS['molpro.wfpath'])
  else:
    string+='\nnooverlap\n'
  sh2pro.write(string)
  sh2pro.close()

  # copy MOs and template
  cpfrom=INFOS['molpro.template']
  cpto='%s/MOLPRO.template' % (iconddir)
  shutil.copy(cpfrom,cpto)
  if INFOS['molpro.guess']:
    cpfrom=INFOS['molpro.guess']
    cpto='%s/wf.init' % (iconddir)
    shutil.copy(cpfrom,cpto)

  return

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def checktemplate_COLUMBUS(TEMPLATE, mult):
  '''Checks whether TEMPLATE is a file or directory. If a file or does not exist, it quits with exit code 1, if it is a directory, it checks whether all important input files are there. Does not check for all input files, since runc does this, too.

  Arguments:
  1 string: path to TEMPLATE

  returns whether input is for isc keyword or socinr keyword
  and returns the DRT of the given multiplicity'''

  exist=os.path.exists(TEMPLATE)
  if exist:
    isfile=os.path.isfile(TEMPLATE)
    if isfile:
      #print 'TEMPLATE=%s exists and is a file!' % (TEMPLATE)
      return None,None,None
    necessary=['control.run','mcscfin','tranin','propin']
    lof=os.listdir(TEMPLATE)
    for i in necessary:
      if not i in lof:
        #print 'Did not find input file %s! Did you prepare the input according to the instructions?' % (i)
        return None,None,None
    cidrtinthere=False
    ciudginthere=False
    for i in lof:
      if 'cidrtin' in i:
        cidrtinthere=True
      if 'ciudgin' in i:
        ciudginthere=True
    if not cidrtinthere or not ciudginthere:
      #print 'Did not find input file %s.*! Did you prepare the input according to the instructions?' % (i)
      return None,None,None
  else:
    #print 'Directory %s does not exist!' % (TEMPLATE)
    return None,None,None

  # get integral program
  try:
    intprog=open(TEMPLATE+'/intprogram')
    line=intprog.readline()
    if 'hermit' in line:
      INTPROG='dalton'
    elif 'seward' in line:
      INTPROG='seward'
    else:
      return None,None,None
  except IOError:
    return None,None,None

  # check cidrtin and cidrtin* for the multiplicity
  try:
    cidrtin=open(TEMPLATE+'/cidrtin')
    line=cidrtin.readline().split()
    if line[0].lower()=='y':
      maxmult=int(cidrtin.readline().split()[0])
      cidrtin.readline()
      nelec=int(cidrtin.readline().split()[0])
      if mult<=maxmult and (mult+nelec)%2!=0:
        return 1, (mult+1)/2,INTPROG    # socinr=1, single=-1, isc=0
      else:
        return None,None,None
    else:
      mult2=int(cidrtin.readline().split()[0])
      if mult!=mult2:
        #print 'Multiplicity %i cannot be treated in directory %s (single DRT)!'  % (mult,TEMPLATE)
        return None,None,None
      return -1,1,INTPROG
  except IOError:
    # find out in which DRT the requested multiplicity is
    for i in range(1,9):        # COLUMBUS can treat at most 8 DRTs
      try:
        cidrtin=open(TEMPLATE+'/cidrtin.%i' % i)
      except IOError:
        return None,None,None
      cidrtin.readline()
      mult2=int(cidrtin.readline().split()[0])
      if mult==mult2:
        return 0,i,INTPROG
      cidrtin.close()

# =================================================

def get_COLUMBUS(INFOS):
  '''This routine asks for all questions specific to COLUMBUS:
  - path to COLUMBUS
  - scratchdir
  - path to template directory
  - mocoef
  - memory
  '''

  string='\n  '+'='*80+'\n'
  string+='||'+centerstring('COLUMBUS Interface setup',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  print string


  print centerstring('Path to COLUMBUS',60,'-')+'\n'
  path=os.getenv('COLUMBUS')
  if path=='':
    path=None
  else:
    path='$COLUMBUS/'
  #path=os.path.expanduser(os.path.expandvars(path))
  #if path!='':
    #print 'Environment variable $COLUMBUS detected:\n$COLUMBUS=%s\n' % (path)
    #if question('Do you want to use this COLUMBUS installation?',bool,True):
      #INFOS['columbus']=path
  #if not 'columbus' in INFOS:
  print '\nPlease specify path to COLUMBUS directory (SHELL variables and ~ can be used, will be expanded when interface is started).\n'
  INFOS['columbus']=question('Path to COLUMBUS:',str,path)
  print ''




  print centerstring('COLUMBUS input template directory',60,'-')+'\n'
  print '''Please specify the path to the COLUMBUS template directory.
The directory must contain subdirectories with complete COLUMBUS input file sets for the following steps:
- Integrals with SEWARD/MOLCAS
- SCF
- MCSCF
- SO-MRCI (even if no Spin-Orbit couplings will be calculated)
The COLUMBUS interface will generate the remaining COLUMBUS input automatically, depending on the number of states.

In order to setup the COLUMBUS input, use COLUMBUS' input facility colinp. For further information, see the Spin-orbit tutorial for COLUMBUS [1].

[1] http://www.univie.ac.at/columbus/docs_COL70/tutorial-SO.pdf
'''
  while True:
    path=question('Path to templates:',str)
    path=os.path.expanduser(os.path.expandvars(path))
    path=os.path.abspath(path)
    if not os.path.isdir(path):
      print 'Directory %s does not exist!' % (path)
      continue

    content=os.listdir(path)
    multmap={}
    allOK=True
    for mult in range(1,1+len(INFOS['states'])):
      if INFOS['states'][mult-1]==0:
        continue
      found=False
      for d in content:
        template=path+'/'+d
        socitype,drt,intprog=checktemplate_COLUMBUS(template,mult)
        if socitype==None:
          continue
        if not d[-1]=='/':
          d+='/'
        multmap[mult]=d
        found=True
        break
      if not found:
        print 'No input directory for multiplicity %i!' % (mult)
        allOK=False
        continue
    if allOK:
      break
  print '\nAccepted path: %s\n' % (path)

  print '''Check whether the jobs are assigned correctly to the multiplicities. Use the following commands:
  mult job        make <mult> use the input in <job>
  show            show the mapping of multiplicities to jobs
  end             confirm this mapping
'''
  for i in multmap:
    print '%i ==> %s' % (i,multmap[i])
  while True:
    line=question('Adjust job mapping:',str,'end',False)
    if 'show' in line.lower():
      for i in multmap:
        print '%i ==> %s' % (i,multmap[i])
      continue
    elif 'end' in line.lower():
      break
    else:
      f=line.split()
      try:
        m=int(f[0])
        j=f[1]
      except (ValueError,IndexError):
        continue
      if not m in multmap:
        print 'Multiplicity %i not necessary!' % (m)
        continue
      if not os.path.isdir(path+'/'+j):
        print 'No template subdirectory %s!' % (j)
        continue
      if not j[-1]=='/':
        j+='/'
      multmap[m]=j
  print ''

  mocoefmap={}
  for job in set([ multmap[i] for i in multmap]):
    mocoefmap[job]=multmap[min(multmap)]
  print '''Check whether the mocoeffiles are assigned correctly to the jobs. Use the following commands:
  job mocoefjob   make <job> use the mocoeffiles from <mocoefjob>
  show            show the mapping of multiplicities to jobs
  end             confirm this mapping
'''
  width=max([ len(i) for i in mocoefmap] )
  for i in mocoefmap:
    print '%s' % (i) +' '*(width-len(i))+ ' <== %s' % (mocoefmap[i])
  while True:
    line=question('Adjust mocoef mapping:',str,'end',False)
    if 'show' in line.lower():
      for i in mocoefmap:
        print '%s <== %s' % (i,mocoefmap[i])
      continue
    elif 'end' in line.lower():
      break
    else:
      f=line.split()
      try:
        j=f[0]
        m=f[1]
      except (ValueError,IndexError):
        continue
      if not m[-1]=='/':
        m+='/'
      if not j[-1]=='/':
        j+='/'
      mocoefmap[j]=m
  print ''

  INFOS['columbus.template']=path
  INFOS['columbus.multmap']=multmap
  INFOS['columbus.mocoefmap']=mocoefmap
  INFOS['columbus.intprog']=intprog

  INFOS['columbus.copy_template']=question('Do you want to copy the template directory to each trajectory (Otherwise it will be linked)?',bool,False)
  if INFOS['columbus.copy_template']:
    INFOS['columbus.copy_template_from']=INFOS['columbus.template']
    INFOS['columbus.template']='./COLUMBUS.template/'


  print centerstring('Initial wavefunction: MO Guess',60,'-')+'\n'
  print '''Please specify the path to a COLUMBUS mocoef file containing suitable starting MOs for the CASSCF calculation.
'''
  init=question('Do you have an initial mocoef file?',bool,True)
  if init:
    while True:
      line=question('Mocoef filename:',str,'mocoef_mc.init')
      line=os.path.expanduser(os.path.expandvars(line))
      if os.path.isfile(line):
          break
      else:
        print 'File not found!'
        continue
    INFOS['columbus.guess']=line
  else:
    print 'WARNING: Remember that CASSCF calculations may run very long and/or yield wrong results without proper starting MOs.'
    time.sleep(2)
    INFOS['columbus.guess']=False
  print ''


  print centerstring('COLUMBUS Memory usage',60,'-')+'\n'
  print '''Please specify the amount of memory available to COLUMBUS (in MB). For calculations including moderately-sized CASSCF calculations and less than 150 basis functions, around 2000 MB should be sufficient.
'''
  INFOS['columbus.mem']=abs(question('COLUMBUS memory:',int)[0])


  # Ionization
  #print '\n'+centerstring('Ionization probability by Dyson norms',60,'-')+'\n'
  #INFOS['ion']=question('Dyson norms?',bool,False)
  #if INFOS['ion']:
  if 'wfoverlap' in INFOS['needed']:
    print '\n'+centerstring('WFoverlap setup',60,'-')+'\n'
    INFOS['columbus.dysonpath']=question('Path to wavefunction overlap executable:',str,'$SHARC/wfoverlap.x')
    INFOS['columbus.ciothres']=question('Determinant screening threshold:',float,[0.97])[0]
    INFOS['columbus.numfrozcore']=question('Number of frozen core orbitals for overlaps (-1=as in template):',int,[-1])[0]
    if 'ion' in INFOS and INFOS['ion']:
      INFOS['columbus.numocc']=question('Number of doubly occupied orbitals for Dyson:',int,[0])[0]

  return INFOS

# =================================================

def prepare_COLUMBUS(INFOS,iconddir):
  # write COLUMBUS.resources
  try:
    sh2col=open('%s/COLUMBUS.resources' % (iconddir), 'w')
  except IOError:
    print 'IOError during prepareCOLUMBUS, directory=%i' % (iconddir)
    quit(1)
  string= 'columbus %s\nscratchdir %s/%s/QM/WORK\n' % (INFOS['columbus'],INFOS['scratchdir'],iconddir)
  string+='savedir %s/%s/savedir\ntemplate %s\nmemory %i\n\n' % (INFOS['scratchdir'],iconddir, INFOS['columbus.template'],INFOS['columbus.mem'])
  string+='integrals %s\n' % (INFOS['columbus.intprog'])
  for mult in INFOS['columbus.multmap']:
    string+='DIR %i %s\n' % (mult,INFOS['columbus.multmap'][mult])
  string+='\n'
  for job in INFOS['columbus.mocoefmap']:
    string+='MOCOEF %s %s\n' % (job,INFOS['columbus.mocoefmap'][job])
  if 'wfoverlap' in INFOS['needed']:
    string+='wfoverlap %s\n' % (INFOS['columbus.dysonpath'])
    string+='wfthres %s\n' % (INFOS['columbus.ciothres'])
    if INFOS['columbus.numfrozcore']>=0:
      string+='numfrozcore %i\n' % (INFOS['columbus.numfrozcore'])
    if 'columbus.numocc' in INFOS:
      string+='numocc %i\n' % (INFOS['columbus.numocc'])
  else:
    string+='\nnooverlap\n'
  sh2col.write(string)
  sh2col.close()

  # copy MOs and template
  if INFOS['columbus.guess']:
    cpfrom=INFOS['columbus.guess']
    cpto='%s/mocoef_mc.init' % (iconddir)
    shutil.copy(cpfrom,cpto)

  if INFOS['columbus.copy_template']:
    copy_from=INFOS['columbus.copy_template_from']
    copy_to=iconddir+'/COLUMBUS.template/'
    if os.path.exists(copy_to):
      shutil.rmtree(copy_to)
    shutil.copytree(copy_from,copy_to)

  return

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================


def checktemplate_MOLCAS(filename,INFOS):
  necessary=['cobramm','basis','ras2','nactel','inactive']
  try:
    f=open(filename)
    data=f.readlines()
    f.close()
  except IOError:
    print 'Could not open template file %s' % (filename)
    return False
  valid=[]
  for i in necessary:
    for l in data:
      if i in re.sub('#.*$','',l):
        valid.append(True)
        break
    else:
      valid.append(False)
  if not all(valid):
    print 'The template %s seems to be incomplete! It should contain: ' % (filename) +str(necessary)
    return False
  roots_there=False
  for l in data:
    l=re.sub('#.*$','',l).lower().split()
    if len(l)==0:
      continue
    if 'roots' in l[0]:
      roots_there=True
  if not roots_there:
    for mult,state in enumerate(INFOS['states']):
      if state<=0:
        continue
      valid=[]
      for l in data:
        if 'spin' in re.sub('#.*$','',l).lower():
          f=l.split()
          if int(f[1])==mult+1:
            valid.append(True)
            break
      else:
        valid.append(False)
  if not all(valid):
    string='The template %s seems to be incomplete! It should contain the keyword "spin" for ' % (filename)
    for mult,state in enumerate(INFOS['states']):
      if state<=0:
        continue
      string+='%s, ' % (IToMult[mult+1])
    string=string[:-2]+'!'
    print string
    return False
  return True

# =================================================

def get_MOLCAS(INFOS):
  '''This routine asks for all questions specific to MOLPRO:
  - path to molpro
  - scratch directory
  - MOLPRO.template
  - wf.init
  '''

  string='\n  '+'='*80+'\n'
  string+='||'+centerstring('MOLCAS Interface setup',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  print string

  print centerstring('Path to MOLCAS',60,'-')+'\n'
  path=os.getenv('MOLCAS')
  #path=os.path.expanduser(os.path.expandvars(path))
  if path=='':
    path=None
  else:
    path='$MOLCAS/'
      #print 'Environment variable $MOLCAS detected:\n$MOLCAS=%s\n' % (path)
      #if question('Do you want to use this MOLCAS installation?',bool,True):
        #INFOS['molcas']=path
    #if not 'molcas' in INFOS:
  print '\nPlease specify path to MOLCAS directory (SHELL variables and ~ can be used, will be expanded when interface is started).\n'
  INFOS['molcas']=question('Path to MOLCAS:',str,path)
  print ''



  print centerstring('MOLCAS input template file',60,'-')+'\n'
  print '''Please specify the path to the MOLCAS.template file. This file must contain the following settings:

basis <Basis set>
ras2 <Number of active orbitals>
nactel <Number of active electrons>
inactive <Number of doubly occupied orbitals>
roots <Number of roots for state-averaging>
cobramm 

The MOLCAS interface will generate the appropriate MOLCAS input automatically.
'''
  if os.path.isfile('MOLCAS.template'):
    if checktemplate_MOLCAS('MOLCAS.template',INFOS):
      print 'Valid file "MOLCAS.template" detected. '
      usethisone=question('Use this template file?',bool,True)
      if usethisone:
        INFOS['molcas.template']='MOLCAS.template'
  if not 'molcas.template' in INFOS:
    while True:
      filename=question('Template filename:',str)
      if not os.path.isfile(filename):
        print 'File %s does not exist!' % (filename)
        continue
      if checktemplate_MOLCAS(filename,INFOS):
        break
    INFOS['molcas.template']=filename
  print ''



  print centerstring('Initial wavefunction: MO Guess',60,'-')+'\n'
  print '''Please specify the path to a MOLCAS JobIph file containing suitable starting MOs for the CASSCF calculation. Please note that this script cannot check whether the wavefunction file and the Input template are consistent!
'''
  string='Do you have initial wavefunction files for '
  for mult,state in enumerate(INFOS['states']):
    if state<=0:
      continue
    string+='%s, ' % (IToMult[mult+1])
  string=string[:-2]+'?'
  if question(string,bool,True):
    while True:
      jobiph_or_rasorb=question('JobIph files (1) or RasOrb files (2)?',int)[0]
      if jobiph_or_rasorb in [1,2]:
        break
    INFOS['molcas.jobiph_or_rasorb']=jobiph_or_rasorb
    INFOS['molcas.guess']={}
    for mult,state in enumerate(INFOS['states']):
      if state<=0:
        continue
      while True:
        if jobiph_or_rasorb==1:
          guess_file='MOLCAS.%i.JobIph.init' % (mult+1)
        else:
          guess_file='MOLCAS.%i.RasOrb.init' % (mult+1)
        filename=question('Initial wavefunction file for %ss:' % (IToMult[mult+1]),str,guess_file)
        if os.path.isfile(filename):
          INFOS['molcas.guess'][mult+1]=filename
          break
        else:
          print 'File not found!'
  else:
    print 'WARNING: Remember that CASSCF calculations may run very long and/or yield wrong results without proper starting MOs.'
    time.sleep(2)
    INFOS['molcas.guess']={}


  print centerstring('MOLCAS Ressource usage',60,'-')+'\n'
  print '''Please specify the amount of memory available to MOLCAS (in MB). For calculations including moderately-sized CASSCF calculations and less than 150 basis functions, around 2000 MB should be sufficient.
'''
  INFOS['molcas.mem']=abs(question('MOLCAS memory:',int,[1000])[0])
  print '''Please specify the number of CPUs to be used by EACH calculation.
'''



  # Ionization
  #print '\n'+centerstring('Ionization probability by Dyson norms',60,'-')+'\n'
  #INFOS['ion']=question('Dyson norms?',bool,False)
  #if INFOS['ion']:
  if 'wfoverlap' in INFOS['needed']:
    print '\n'+centerstring('WFoverlap setup',60,'-')+'\n'
    INFOS['molcas.wfoverlap']=question('Path to wavefunction overlap executable:',str,'$SHARC/wfoverlap.x')
    # TODO not asked for: numfrozcore, numocc

  return INFOS

# =================================================

def prepare_MOLCAS(INFOS,iconddir):
  # write MOLCAS.resources
  try:
    sh2cas=open('%s/MOLCAS.resources' % (iconddir), 'w')
  except IOError:
    print 'IOError during prepareMOLCAS, iconddir=%s' % (iconddir)
    quit(1)
  project='MOLCAS'
  string='molcas %s\nscratchdir %s/%s/QMMM/QM\nmemory %i\nncpu %i\nproject %s' % (INFOS['molcas'],INFOS['scratchdir'],iconddir,INFOS['molcas.mem'],INFOS['cobramm.ncpu'],project)
  if 'wfoverlap' in INFOS['needed']:
    string+='\nwfoverlap %s\n' % INFOS['molcas.wfoverlap']
  else:
    string+='\nnooverlap\n'
  if 'tinker' in INFOS:
    string+='tinker %s' % (INFOS['tinker'])
  sh2cas.write(string)
  sh2cas.close()

  # copy MOs and template
  cpfrom=INFOS['molcas.template']
  cpto='%s/MOLCAS.template' % (iconddir)
  shutil.copy(cpfrom,cpto)
  if not INFOS['molcas.guess']=={}:
    for i in INFOS['molcas.guess']:
      if INFOS['molcas.jobiph_or_rasorb']==1:
        cpfrom=INFOS['molcas.guess'][i]
        cpto='%s/%s.%i.JobIph.init' % (iconddir,project,i)
      else:
        cpfrom=INFOS['molcas.guess'][i]
        cpto='%s/%s.%i.RasOrb.init' % (iconddir,project,i)
      shutil.copy(cpfrom,cpto)


#======================================================================================================================
#======================================================================================================================
#======================================================================================================================

def checktemplate_ADF(filename,INFOS):
  necessary=['cobramm','basis','functional','charge']
  try:
    f=open(filename)
    data=f.readlines()
    f.close()
  except IOError:
    print 'Could not open template file %s' % (filename)
    return False
  valid=[]
  for i in necessary:
    for l in data:
      line=l.lower().split()
      if len(line)==0:
        continue
      line=line[0]
      if i==re.sub('#.*$','',line):
        valid.append(True)
        break
    else:
      valid.append(False)
  if not all(valid):
    print 'The template %s seems to be incomplete! It should contain: ' % (filename) +str(necessary)
    return False
  return True

# =================================================

def get_ADF(INFOS):
  '''This routine asks for all questions specific to ADF:
  - path to ADF
  - scratch directory
  - ADF.template
  - TAPE21
  '''

  string='\n  '+'='*80+'\n'
  string+='||'+centerstring('ADF Interface setup',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  print string

  print centerstring('Path to ADF',60,'-')+'\n'
  path=os.getenv('ADFHOME')
  if path:
    path='$ADFHOME/'
  adfrc=question('Setup from adfrc.sh file?',bool,True)
  if adfrc:
    if path:
      path='$ADFHOME/adfrc.sh'
    print '\nPlease specify path to the adfrc.sh file (SHELL variables and ~ can be used, will be expanded when interface is started).\n'
    path=question('Path to adfrc.sh file:',str,path)
    INFOS['adfrc']=os.path.abspath(os.path.expanduser(os.path.expandvars(path)))
    print 'Will use adfrc= %s' % INFOS['adfrc']
    INFOS['adf']='$ADFHOME'
    INFOS['scmlicense']='$SCMLICENSE'
    print ''
  else:
    print '\nPlease specify path to ADF directory (SHELL variables and ~ can be used, will be expanded when interface is started).\n'
    INFOS['adf']=question('Path to ADF:',str,path)
    print ''
    print centerstring('Path to ADF license file',60,'-')+'\n'
    path=os.getenv('SCMLICENSE')
    #path=os.path.expanduser(os.path.expandvars(path))
    if path=='':
      path=None
    else:
      path='$SCMLICENSE'
    print'\nPlease specify path to ADF license.txt\n'
    INFOS['scmlicense']=question('Path to license:',str,path)
    print ''




  # template file
  print centerstring('ADF input template file',60,'-')+'\n'
  print '''Please specify the path to the ADF.template file. This file must contain the following keywords:

basis <basis>
functional <type> <name>
charge <x> [ <x2> [ <x3> ...] ]

The ADF interface will generate the appropriate ADF input automatically.
'''
  if os.path.isfile('ADF.template'):
    if checktemplate_ADF('ADF.template',INFOS):
      print 'Valid file "ADF.template" detected. '
      usethisone=question('Use this template file?',bool,True)
      if usethisone:
        INFOS['ADF.template']='ADF.template'
  if not 'ADF.template' in INFOS:
    while True:
      filename=question('Template filename:',str)
      if not os.path.isfile(filename):
        print 'File %s does not exist!' % (filename)
        continue
      if checktemplate_ADF(filename,INFOS):
        break
    INFOS['ADF.template']=filename
  print ''



  # initial MOs
  print centerstring('Initial restart: MO Guess',60,'-')+'\n'
  print '''Please specify the path to an ADF TAPE21 file containing suitable starting MOs for the ADF calculation. Please note that this script cannot check whether the wavefunction file and the Input template are consistent!
'''
  if question('Do you have a restart file?',bool,True):
    if True:
      while True:
        filename=question('Restart file:',str,'ADF.t21.init')
        if os.path.isfile(filename):
          INFOS['adf.guess']=filename
          break
        else:
          print 'Could not find file "%s"!' % (filename)
  else:
    INFOS['adf.guess']={}


  # Resources
  print centerstring('ADF Ressource usage',60,'-')+'\n'
  print '''Please specify the number of CPUs to be used by EACH calculation.
'''
  INFOS['adf.ncpu']=abs(question('Number of CPUs:',int)[0])

  if INFOS['adf.ncpu']>1:
    print '''Please specify how well your job will parallelize.
A value of 0 means that running in parallel will not make the calculation faster, a value of 1 means that the speedup scales perfectly with the number of cores.
Typical values for ADF are 0.90-0.98 for LDA/GGA functionals and 0.50-0.80 for hybrids (better if RIHartreeFock is used).'''
    INFOS['adf.scaling']=min(1.0,max(0.0,question('Parallel scaling:',float,[0.8])[0] ))
  else:
    INFOS['adf.scaling']=0.9


  # Ionization
  #print '\n'+centerstring('Ionization probability by Dyson norms',60,'-')+'\n'
  #INFOS['ion']=question('Dyson norms?',bool,False)
  #if INFOS['ion']:
  if 'wfoverlap' in INFOS['needed']:
    print '\n'+centerstring('WFoverlap setup',60,'-')+'\n'
    INFOS['adf.wfoverlap']=question('Path to wavefunction overlap executable:',str,'$SHARC/wfoverlap.x')
    print ''
    print 'State threshold for choosing determinants to include in the overlaps'
    print 'For hybrids (and without TDA) one should consider that the eigenvector X may have a norm larger than 1'
    INFOS['adf.ciothres']=question('Threshold:',float,[0.99])[0]
    print ''
    INFOS['adf.mem']=question('Memory for wfoverlap (MB):',int,[1000])[0]
    # TODO not asked: numfrozcore and numocc

    #print 'Please state the number of core orbitals you wish to freeze for the overlaps (recommended to use for at least the 1s orbital and a negative number uses default values)?'
    #print 'A value of -1 will use the defaults used by ADF for a small frozen core and 0 will turn off the use of frozen cores'
    #INFOS['frozcore_number']=question('How many orbital to freeze?',int,[-1])[0]


  # TheoDORE
  theodore_spelling=['Om',
                    'PRNTO',
                    'Z_HE', 'S_HE', 'RMSeh',
                    'POSi', 'POSf', 'POS',
                    'PRi', 'PRf', 'PR', 'PRh',
                    'CT', 'CT2', 'CTnt',
                    'MC', 'LC', 'MLCT', 'LMCT', 'LLCT',
                    'DEL', 'COH', 'COHh']
  #INFOS['theodore']=question('TheoDORE analysis?',bool,False)
  if 'theodore' in INFOS['needed']:
    print '\n'+centerstring('Wave function analysis by TheoDORE',60,'-')+'\n'

    INFOS['adf.theodore']=question('Path to TheoDORE directory:',str,'$THEODIR')
    print ''

    print 'Please give a list of the properties to calculate by TheoDORE.\nPossible properties:'
    string=''
    for i,p in enumerate(theodore_spelling):
      string+='%s ' % (p)
      if (i+1)%8==0:
        string+='\n'
    print string
    l=question('TheoDORE properties:',str,'Om  PRNTO  S_HE  Z_HE  RMSeh')
    if '[' in l:
      INFOS['theodore.prop']=ast.literal_eval(l)
    else:
      INFOS['theodore.prop']=l.split()
    print ''

    print 'Please give a list of the fragments used for TheoDORE analysis.'
    print 'You can use the list-of-lists from dens_ana.in'
    print 'Alternatively, enter all atom numbers for one fragment in one line. After defining all fragments, type "end".'
    if qmmm_job(INFOS['ADF.template'],INFOS):
      print 'You should only include the atom numbers of QM and link atoms.'
    INFOS['theodore.frag']=[]
    while True:
      l=question('TheoDORE fragment:',str,'end')
      if 'end' in l.lower():
        break
      if '[' in l:
        try:
          INFOS['theodore.frag']=ast.literal_eval(l)
          break
        except ValueError:
          continue
      f=[ int(i) for i in l.split() ]
      INFOS['theodore.frag'].append(f)
    INFOS['theodore.count']=len(INFOS['theodore.prop'])+len(INFOS['theodore.frag'])**2


  return INFOS

# =================================================

def prepare_ADF(INFOS,iconddir):
  # write ADF.resources
  try:
    sh2cas=open('%s/ADF.resources' % (iconddir), 'w')
  except IOError:
    print 'IOError during prepareADF, iconddir=%s' % (iconddir)
    quit(1)
#  project='ADF'
  string='adfhome %s\nscmlicense %s\nscratchdir %s/%s/QM\nncpu %i\nschedule_scaling %f\n' % (INFOS['adf'],INFOS['scmlicense'],INFOS['scratchdir'],iconddir,INFOS['adf.ncpu'],INFOS['adf.scaling'])
  if 'wfoverlap' in INFOS['needed']:
    string+='wfoverlap %s\nwfthres %f\n' % (INFOS['adf.wfoverlap'],INFOS['adf.ciothres'])
    string+='memory %i\n' % (INFOS['adf.mem'])
    #string+='numfrozcore %i\n' %(INFOS['frozcore_number'])
  else:
    string+='nooverlap\n'
  if 'theodore' in INFOS['needed']:
    string+='theodir %s\n' % (INFOS['adf.theodore'])
    string+='theodore_prop %s\n' % (INFOS['theodore.prop'])
    string+='theodore_fragment %s\n' % (INFOS['theodore.frag'])
  if 'ADF.fffile' in INFOS:
    string+='qmmm_ff_file ADF.qmmm.ff\n'
  if 'ADF.ctfile' in INFOS:
    string+='qmmm_table ADF.qmmm.table\n'
  sh2cas.write(string)
  sh2cas.close()

  # copy MOs and template
  cpfrom=INFOS['ADF.template']
  cpto='%s/ADF.template' % (iconddir)
  shutil.copy(cpfrom,cpto)

  if INFOS['adf.guess']:
    cpfrom1=INFOS['adf.guess']
    cpto1='%s/ADF.t21_init' % (iconddir)
    shutil.copy(cpfrom1,cpto1)

  if 'ADF.fffile' in INFOS:
    cpfrom1=INFOS['ADF.fffile']
    cpto1='%s/ADF.qmmm.ff' % (iconddir)
    shutil.copy(cpfrom1,cpto1)

  if 'ADF.ctfile' in INFOS:
    cpfrom1=INFOS['ADF.ctfile']
    cpto1='%s/ADF.qmmm.table' % (iconddir)
    shutil.copy(cpfrom1,cpto1)

  return

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def checktemplate_RICC2(filename,INFOS):
  necessary=['cobramm','basis']
  try:
    f=open(filename)
    data=f.readlines()
    f.close()
  except IOError:
    print 'Could not open template file %s' % (filename)
    return False
  valid=[]
  for i in necessary:
    for l in data:
      line=l.lower()
      if i in re.sub('#.*$','',line):
        valid.append(True)
        break
    else:
      valid.append(False)
  if not all(valid):
    print 'The template %s seems to be incomplete! It should contain: ' % (filename) +str(necessary)
    return False
  return True

# =================================================

def get_RICC2(INFOS):
  string='\n  '+'='*80+'\n'
  string+='||'+centerstring('Turbomole RICC2 Interface setup',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  print string

  print centerstring('Path to TURBOMOLE',60,'-')+'\n'
  path=os.getenv('TURBODIR')
  if path=='':
    path=None
  else:
    path='$TURBODIR/'
  print '\nPlease specify path to TURBOMOLE directory (SHELL variables and ~ can be used, will be expanded when interface is started).\n'
  INFOS['turbomole']=question('Path to TURBOMOLE:',str,path)
  print ''

  if INFOS['soc']:
    print centerstring('Path to ORCA',60,'-')+'\n'
    path=os.getenv('ORCADIR')
    if path=='':
      path=None
    else:
      path='$ORCADIR/'
    print '\nPlease specify path to ORCA directory (SHELL variables and ~ can be used, will be expanded when interface is started).\n\nORCA is necessary for the calculation of spin-orbit couplings with ricc2.\n'
    INFOS['orca']=question('Path to ORCA:',str,path)
    print ''


  print centerstring('RICC2 input template file',60,'-')+'\n'
  print '''Please specify the path to the RICC2.template file. This file must contain the following settings:

basis <Basis set>

In addition, it can contain the following:

auxbasis <Basis set>
charge <integer>
method <"ADC(2)" or "CC2">                      # only ADC(2) can calculate spin-orbit couplings
frozen <number of frozen core orbitals>
spin-scaling <"none", "SCS", or "SOS">
douglas-kroll                                   # DKH is only used if this keyword is given

'''
  if os.path.isfile('RICC2.template'):
    if checktemplate_RICC2('RICC2.template',INFOS):
      print 'Valid file "RICC2.template" detected. '
      usethisone=question('Use this template file?',bool,True)
      if usethisone:
        INFOS['ricc2.template']='RICC2.template'
  if not 'ricc2.template' in INFOS:
    while True:
      filename=question('Template filename:',str)
      if not os.path.isfile(filename):
        print 'File %s does not exist!' % (filename)
        continue
      if checktemplate_RICC2(filename,INFOS):
        break
    INFOS['ricc2.template']=filename
  print ''




  print centerstring('Initial wavefunction: MO Guess',60,'-')+'\n'
  print '''Please specify the path to a Turbomole "mos" file containing suitable starting MOs for the calculation. Please note that this script cannot check whether the file and the input template are consistent!
'''
  string='Do you have an initial orbitals file?'
  if question(string,bool,True):
    while True:
      guess_file='mos'
      filename=question('Initial wavefunction file:',str,guess_file)
      if os.path.isfile(filename):
        INFOS['ricc2.guess']=filename
        break
      else:
        print 'File not found!'
  else:
    INFOS['ricc2.guess']=[]


  print centerstring('RICC2 Ressource usage',60,'-')+'\n'
  print '''Please specify the amount of memory available to Turbomole (in MB).
'''
  INFOS['ricc2.mem']=abs(question('RICC2 memory:',int,[1000])[0])
  print '''Please specify the number of CPUs to be used by EACH trajectory.
'''



  if 'wfoverlap' in INFOS['needed']:
    print '\n'+centerstring('WFoverlap setup',60,'-')+'\n'
    INFOS['ricc2.wfoverlap']=question('Path to wavefunction overlap executable:',str,'$SHARC/wfoverlap.x')
    print ''
    print 'State threshold for choosing determinants to include in the overlaps'
    #print 'For hybrids (and without TDA) one should consider that the eigenvector X may have a norm larger than 1'
    INFOS['ricc2.ciothres']=question('Threshold:',float,[0.99])[0]



  # TheoDORE
  theodore_spelling=['Om',
                    'PRNTO',
                    'Z_HE', 'S_HE', 'RMSeh',
                    'POSi', 'POSf', 'POS',
                    'PRi', 'PRf', 'PR', 'PRh',
                    'CT', 'CT2', 'CTnt',
                    'MC', 'LC', 'MLCT', 'LMCT', 'LLCT',
                    'DEL', 'COH', 'COHh']
  #INFOS['theodore']=question('TheoDORE analysis?',bool,False)
  if 'theodore' in INFOS['needed']:
    print '\n'+centerstring('Wave function analysis by TheoDORE',60,'-')+'\n'

    INFOS['ricc2.theodore']=question('Path to TheoDORE directory:',str,'$THEODIR')
    print ''

    print 'Please give a list of the properties to calculate by TheoDORE.\nPossible properties:'
    string=''
    for i,p in enumerate(theodore_spelling):
      string+='%s ' % (p)
      if (i+1)%8==0:
        string+='\n'
    print string
    l=question('TheoDORE properties:',str,'Om  PRNTO  S_HE  Z_HE  RMSeh')
    if '[' in l:
      INFOS['theodore.prop']=ast.literal_eval(l)
    else:
      INFOS['theodore.prop']=l.split()
    print ''

    print 'Please give a list of the fragments used for TheoDORE analysis.'
    print 'You can use the list-of-lists from dens_ana.in'
    print 'Alternatively, enter all atom numbers for one fragment in one line. After defining all fragments, type "end".'
    INFOS['theodore.frag']=[]
    while True:
      l=question('TheoDORE fragment:',str,'end')
      if 'end' in l.lower():
        break
      if '[' in l:
        try:
          INFOS['theodore.frag']=ast.literal_eval(l)
          break
        except ValueError:
          continue
      f=[ int(i) for i in l.split() ]
      INFOS['theodore.frag'].append(f)
    INFOS['theodore.count']=len(INFOS['theodore.prop'])+len(INFOS['theodore.frag'])**2




  return INFOS

# =================================================

def prepare_RICC2(INFOS,iconddir):
  # write RICC2.resources
  try:
    sh2cas=open('%s/RICC2.resources' % (iconddir), 'w')
  except IOError:
    print 'IOError during prepare_RICC2, iconddir=%s' % (iconddir)
    quit(1)
  string='''turbodir %s
scratchdir %s/%s/QM/
memory %i
ncpu %i
dipolelevel 1
''' % (INFOS['turbomole'],
       INFOS['scratchdir'],
       iconddir,
       INFOS['ricc2.mem'],
       INFOS['cobramm.ncpu'])
  if INFOS['soc']:
    string+='orcadir %s\n' % (INFOS['orca'])
  if 'wfoverlap' in INFOS['needed']:
    string+='wfoverlap %s\nwfthres %f\n' % (INFOS['ricc2.wfoverlap'],INFOS['ricc2.ciothres'])
    #string+='numfrozcore %i\n' %(INFOS['frozcore_number'])
  else:
    string+='nooverlap\n'
  if 'theodore' in INFOS['needed']:
    string+='theodir %s\n' % (INFOS['ricc2.theodore'])
    string+='theodore_prop %s\n' % (INFOS['theodore.prop'])
    string+='theodore_fragment %s\n' % (INFOS['theodore.frag'])
  if 'tinker' in INFOS:
    string+='tinker %s\n' % (INFOS['tinker'])
  if 'RICC2.fffile' in INFOS:
    string+='qmmm_ff_file RICC2.qmmm.ff\n'
  if 'RICC2.ctfile' in INFOS:
    string+='qmmm_table RICC2.qmmm.table\n'
  sh2cas.write(string)
  sh2cas.close()

  # copy MOs and template
  cpfrom=INFOS['ricc2.template']
  cpto='%s/RICC2.template' % (iconddir)
  shutil.copy(cpfrom,cpto)
  if INFOS['ricc2.guess']:
    cpfrom1=INFOS['ricc2.guess']
    cpto1='%s/mos.init' % (iconddir)
    shutil.copy(cpfrom1,cpto1)

  if 'RICC2.fffile' in INFOS:
    cpfrom1=INFOS['RICC2.fffile']
    cpto1='%s/RICC2.qmmm.ff' % (iconddir)
    shutil.copy(cpfrom1,cpto1)

  if 'RICC2.ctfile' in INFOS:
    cpfrom1=INFOS['RICC2.ctfile']
    cpto1='%s/RICC2.qmmm.table' % (iconddir)
    shutil.copy(cpfrom1,cpto1)
  return

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def checktemplate_GAUSSIAN(filename,INFOS):
  necessary=['cobramm','basis','functional','charge']
  try:
    f=open(filename)
    data=f.readlines()
    f.close()
  except IOError:
    print 'Could not open template file %s' % (filename)
    return False
  valid=[]
  for i in necessary:
    for l in data:
      line=l.lower().split()
      if len(line)==0:
        continue
      line=line[0]
      if i==re.sub('#.*$','',line):
        valid.append(True)
        break
    else:
      valid.append(False)
  if not all(valid):
    print 'The template %s seems to be incomplete! It should contain: ' % (filename) +str(necessary)
    return False
  return True

# =================================================

def get_GAUSSIAN(INFOS):
  '''This routine asks for all questions specific to GAUSSIAN:
  - path to GAUSSIAN
  - scratch directory
  - GAUSSIAN.template
  '''

  string='\n  '+'='*80+'\n'
  string+='||'+centerstring('GAUSSIAN Interface setup',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  print string

  print centerstring('Path to GAUSSIAN',60,'-')+'\n'
  tries=['g16root','g09root','g03root']
  for i in tries:
    path=os.getenv(i)
    if path:
      path='$%s/' % i
      break
  #gaussianprofile=question('Setup from gaussian.profile file?',bool,True)
  #if gaussianprofile:
    #if path:
      #path='%s/gaussian.profile' % path
    #print '\nPlease specify path to the gaussian.profile file (SHELL variables and ~ can be used, will be expanded when interface is started).\n'
    #path=question('Path to GAUSSIAN:',str,path)
    #INFOS['gaussianprofile']=os.path.abspath(os.path.expanduser(os.path.expandvars(path)))
    #print 'Will use gaussianprofile= %s' % INFOS['gaussianprofile']
    #INFOS['gaussian']='$GAUSSIANHOME'
    #print ''
  #else:
  print '\nPlease specify path to GAUSSIAN directory (SHELL variables and ~ can be used, will be expanded when interface is started).\n'
  INFOS['groot']=question('Path to GAUSSIAN:',str,path)
  print ''


  # scratch

  # template file
  print centerstring('GAUSSIAN input template file',60,'-')+'\n'
  print '''Please specify the path to the GAUSSIAN.template file. This file must contain the following keywords:
  
basis <basis>
functional <type> <name>
charge <x> [ <x2> [ <x3> ...] ] 

The GAUSSIAN interface will generate the appropriate GAUSSIAN input automatically.
'''
  if os.path.isfile('GAUSSIAN.template'):
    if checktemplate_GAUSSIAN('GAUSSIAN.template',INFOS):
      print 'Valid file "GAUSSIAN.template" detected. '
      usethisone=question('Use this template file?',bool,True)
      if usethisone:
        INFOS['GAUSSIAN.template']='GAUSSIAN.template'
  if not 'GAUSSIAN.template' in INFOS:
    while True:
      filename=question('Template filename:',str)
      if not os.path.isfile(filename):
        print 'File %s does not exist!' % (filename)
        continue
      if checktemplate_GAUSSIAN(filename,INFOS):
        break
    INFOS['GAUSSIAN.template']=filename
  print ''



  # initial MOs
  print centerstring('Initial restart: MO Guess',60,'-')+'\n'
  print '''Please specify the path to an GAUSSIAN chk file containing suitable starting MOs for the GAUSSIAN calculation. Please note that this script cannot check whether the wavefunction file and the Input template are consistent!
'''
  if question('Do you have a restart file?',bool,True):
    if True:
      while True:
        filename=question('Restart file:',str,'GAUSSIAN.chk.init')
        if os.path.isfile(filename):
          INFOS['gaussian.guess']=filename
          break
        else:
          print 'Could not find file "%s"!' % (filename)
  else:
    INFOS['gaussian.guess']={}


  # Resources
  print centerstring('GAUSSIAN Ressource usage',60,'-')+'\n'
  print '''Please specify the number of CPUs to be used by EACH calculation.
'''
  INFOS['gaussian.ncpu']=abs(question('Number of CPUs:',int)[0])

  if INFOS['gaussian.ncpu']>1:
    print '''Please specify how well your job will parallelize.
A value of 0 means that running in parallel will not make the calculation faster, a value of 1 means that the speedup scales perfectly with the number of cores.
Typical values for GAUSSIAN are 0.90-0.98.'''
    INFOS['gaussian.scaling']=min(1.0,max(0.0,question('Parallel scaling:',float,[0.9])[0] ))
  else:
    INFOS['gaussian.scaling']=0.9

  INFOS['gaussian.mem']=question('Memory (MB):',int,[1000])[0]

  # Ionization
  #print '\n'+centerstring('Ionization probability by Dyson norms',60,'-')+'\n'
  #INFOS['ion']=question('Dyson norms?',bool,False)
  #if INFOS['ion']:
  if 'wfoverlap' in INFOS['needed']:
    print '\n'+centerstring('WFoverlap setup',60,'-')+'\n'
    INFOS['gaussian.wfoverlap']=question('Path to wavefunction overlap executable:',str,'$SHARC/wfoverlap.x')
    print ''
    print 'State threshold for choosing determinants to include in the overlaps'
    print 'For hybrids without TDA one should consider that the eigenvector X may have a norm larger than 1'
    INFOS['gaussian.ciothres']=question('Threshold:',float,[0.99])[0]
    print ''
    # TODO not asked: numfrozcore and numocc

    #print 'Please state the number of core orbitals you wish to freeze for the overlaps (recommended to use for at least the 1s orbital and a negative number uses default values)?'
    #print 'A value of -1 will use the defaults used by GAUSSIAN for a small frozen core and 0 will turn off the use of frozen cores'
    #INFOS['frozcore_number']=question('How many orbital to freeze?',int,[-1])[0]


  # TheoDORE
  theodore_spelling=['Om', 
                    'PRNTO', 
                    'Z_HE', 'S_HE', 'RMSeh',
                    'POSi', 'POSf', 'POS', 
                    'PRi', 'PRf', 'PR', 'PRh',
                    'CT', 'CT2', 'CTnt',
                    'MC', 'LC', 'MLCT', 'LMCT', 'LLCT', 
                    'DEL', 'COH', 'COHh']
  #INFOS['theodore']=question('TheoDORE analysis?',bool,False)
  if 'theodore' in INFOS['needed']:
    print '\n'+centerstring('Wave function analysis by TheoDORE',60,'-')+'\n'

    INFOS['gaussian.theodore']=question('Path to TheoDORE directory:',str,'$THEODIR')
    print ''

    print 'Please give a list of the properties to calculate by TheoDORE.\nPossible properties:'
    string=''
    for i,p in enumerate(theodore_spelling):
      string+='%s ' % (p)
      if (i+1)%8==0:
        string+='\n'
    print string
    l=question('TheoDORE properties:',str,'Om  PRNTO  S_HE  Z_HE  RMSeh')
    if '[' in l:
      INFOS['theodore.prop']=ast.literal_eval(l)
    else:
      INFOS['theodore.prop']=l.split()
    print ''

    print 'Please give a list of the fragments used for TheoDORE analysis.'
    print 'You can use the list-of-lists from dens_ana.in'
    print 'Alternatively, enter all atom numbers for one fragment in one line. After defining all fragments, type "end".'
    if qmmm_job(INFOS['GAUSSIAN.template'],INFOS):
      print 'You should only include the atom numbers of QM and link atoms.'
    INFOS['theodore.frag']=[]
    while True:
      l=question('TheoDORE fragment:',str,'end')
      if 'end' in l.lower():
        break
      if '[' in l:
        try:
          INFOS['theodore.frag']=ast.literal_eval(l)
          break
        except ValueError:
          continue
      f=[ int(i) for i in l.split() ]
      INFOS['theodore.frag'].append(f)
    INFOS['theodore.count']=len(INFOS['theodore.prop'])+len(INFOS['theodore.frag'])**2


  return INFOS

# =================================================

def prepare_GAUSSIAN(INFOS,iconddir):
  # write GAUSSIAN.resources
  try:
    sh2cas=open('%s/GAUSSIAN.resources' % (iconddir), 'w')
  except IOError:
    print 'IOError during prepareGAUSSIAN, iconddir=%s' % (iconddir)
    quit(1)
#  project='GAUSSIAN'
  string='groot %s\nscratchdir %s/%s/QM\nncpu %i\nschedule_scaling %f\n' % (INFOS['groot'],INFOS['scratchdir'],iconddir,INFOS['gaussian.ncpu'],INFOS['gaussian.scaling'])
  string+='memory %i\n' % (INFOS['gaussian.mem'])
  if 'wfoverlap' in INFOS['needed']:
    string+='wfoverlap %s\nwfthres %f\n' % (INFOS['gaussian.wfoverlap'],INFOS['gaussian.ciothres'])
    #string+='numfrozcore %i\n' %(INFOS['frozcore_number'])
  else:
    string+='nooverlap\n'
  if 'theodore' in INFOS['needed']:
    string+='theodir %s\n' % (INFOS['gaussian.theodore'])
    string+='theodore_prop %s\n' % (INFOS['theodore.prop'])
    string+='theodore_fragment %s\n' % (INFOS['theodore.frag'])
  sh2cas.write(string)
  sh2cas.close()

  # copy MOs and template
  cpfrom=INFOS['GAUSSIAN.template']
  cpto='%s/GAUSSIAN.template' % (iconddir)
  shutil.copy(cpfrom,cpto)

  if INFOS['gaussian.guess']:
    cpfrom1=INFOS['gaussian.guess']
    cpto1='%s/GAUSSIAN.chk.init' % (iconddir)
    shutil.copy(cpfrom1,cpto1)

  return


#======================================================================================================================
#======================================================================================================================
#======================================================================================================================

def checktemplate_ORCA(filename,INFOS):
  necessary=['cobramm','basis','functional','charge']
  try:
    f=open(filename)
    data=f.readlines()
    f.close()
  except IOError:
    print 'Could not open template file %s' % (filename)
    return False
  valid=[]
  for i in necessary:
    for l in data:
      line=l.lower().split()
      if len(line)==0:
        continue
      line=line[0]
      if i==re.sub('#.*$','',line):
        valid.append(True)
        break
    else:
      valid.append(False)
  if not all(valid):
    print 'The template %s seems to be incomplete! It should contain: ' % (filename) +str(necessary)
    return False
  return True

# =================================================

def get_ORCA(INFOS):
  '''This routine asks for all questions specific to ORCA:
  - path to ORCA
  - scratch directory
  - ORCA.template
  - initial gbw file
  '''

  string='\n  '+'='*80+'\n'
  string+='||'+centerstring('ORCA Interface setup',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  print string

  print centerstring('Path to ORCA',60,'-')+'\n'
  print '\nPlease specify path to ORCA directory (SHELL variables and ~ can be used, will be expanded when interface is started).\n'
  INFOS['orcadir']=question('Path to ORCA:',str,'$ORCADIR')
  print ''





  # template file
  print centerstring('ORCA input template file',60,'-')+'\n'
  print '''Please specify the path to the ORCA.template file. This file must contain the following keywords:

basis <basis>
functional <type> <name>
charge <x> [ <x2> [ <x3> ...] ]

The ORCA interface will generate the appropriate ORCA input automatically.
'''
  if os.path.isfile('ORCA.template'):
    if checktemplate_ORCA('ORCA.template',INFOS):
      print 'Valid file "ORCA.template" detected. '
      usethisone=question('Use this template file?',bool,True)
      if usethisone:
        INFOS['ORCA.template']='ORCA.template'
  if not 'ORCA.template' in INFOS:
    while True:
      filename=question('Template filename:',str)
      if not os.path.isfile(filename):
        print 'File %s does not exist!' % (filename)
        continue
      if checktemplate_ORCA(filename,INFOS):
        break
    INFOS['ORCA.template']=filename
  print ''




  # initial MOs
  print centerstring('Initial restart: MO Guess',60,'-')+'\n'
  print '''Please specify the path to an ORCA gbw file containing suitable starting MOs for the ORCA calculation. Please note that this script cannot check whether the wavefunction file and the Input template are consistent!
'''
  if question('Do you have a restart file?',bool,True):
    if True:
      while True:
        filename=question('Restart file:',str,'ORCA.gbw')
        if os.path.isfile(filename):
          INFOS['orca.guess']=filename
          break
        else:
          print 'Could not find file "%s"!' % (filename)
  else:
    INFOS['orca.guess']={}


  # Resources
  print centerstring('ORCA Ressource usage',60,'-')+'\n'
  print '''Please specify the number of CPUs to be used by EACH calculation.
'''
  INFOS['orca.ncpu']=abs(question('Number of CPUs:',int)[0])

  if INFOS['orca.ncpu']>1:
    print '''Please specify how well your job will parallelize.
A value of 0 means that running in parallel will not make the calculation faster, a value of 1 means that the speedup scales perfectly with the number of cores.'''
    INFOS['orca.scaling']=min(1.0,max(0.0,question('Parallel scaling:',float,[0.8])[0] ))
  else:
    INFOS['orca.scaling']=0.9
  INFOS['orca.mem']=question('Memory (MB):',int,[1000])[0]


  # Ionization
  #print '\n'+centerstring('Ionization probability by Dyson norms',60,'-')+'\n'
  #INFOS['ion']=question('Dyson norms?',bool,False)
  #if INFOS['ion']:
  if 'wfoverlap' in INFOS['needed']:
    print '\n'+centerstring('WFoverlap setup',60,'-')+'\n'
    INFOS['orca.wfoverlap']=question('Path to wavefunction overlap executable:',str,'$SHARC/wfoverlap.x')
    print ''
    print 'State threshold for choosing determinants to include in the overlaps'
    print 'For hybrids (and without TDA) one should consider that the eigenvector X may have a norm larger than 1'
    INFOS['orca.ciothres']=question('Threshold:',float,[0.99])[0]
    print ''



  # TheoDORE
  theodore_spelling=['Om',
                    'PRNTO',
                    'Z_HE', 'S_HE', 'RMSeh',
                    'POSi', 'POSf', 'POS',
                    'PRi', 'PRf', 'PR', 'PRh',
                    'CT', 'CT2', 'CTnt',
                    'MC', 'LC', 'MLCT', 'LMCT', 'LLCT',
                    'DEL', 'COH', 'COHh']
  #INFOS['theodore']=question('TheoDORE analysis?',bool,False)
  if 'theodore' in INFOS['needed']:
    print '\n'+centerstring('Wave function analysis by TheoDORE',60,'-')+'\n'

    INFOS['orca.theodore']=question('Path to TheoDORE directory:',str,'$THEODIR')
    print ''

    print 'Please give a list of the properties to calculate by TheoDORE.\nPossible properties:'
    string=''
    for i,p in enumerate(theodore_spelling):
      string+='%s ' % (p)
      if (i+1)%8==0:
        string+='\n'
    print string
    l=question('TheoDORE properties:',str,'Om  PRNTO  S_HE  Z_HE  RMSeh')
    if '[' in l:
      INFOS['theodore.prop']=ast.literal_eval(l)
    else:
      INFOS['theodore.prop']=l.split()
    print ''

    print 'Please give a list of the fragments used for TheoDORE analysis.'
    print 'You can use the list-of-lists from dens_ana.in'
    print 'Alternatively, enter all atom numbers for one fragment in one line. After defining all fragments, type "end".'
    if qmmm_job(INFOS['ORCA.template'],INFOS):
      print 'You should only include the atom numbers of QM and link atoms.'
    INFOS['theodore.frag']=[]
    while True:
      l=question('TheoDORE fragment:',str,'end')
      if 'end' in l.lower():
        break
      if '[' in l:
        try:
          INFOS['theodore.frag']=ast.literal_eval(l)
          break
        except ValueError:
          continue
      f=[ int(i) for i in l.split() ]
      INFOS['theodore.frag'].append(f)
    INFOS['theodore.count']=len(INFOS['theodore.prop'])+len(INFOS['theodore.frag'])**2
    if 'ORCA.ctfile' in INFOS:
        INFOS['theodore.count']+=6


  return INFOS

# =================================================

def prepare_ORCA(INFOS,iconddir):
  # write ORCA.resources
  try:
    sh2cas=open('%s/ORCA.resources' % (iconddir), 'w')
  except IOError:
    print 'IOError during prepareORCA, iconddir=%s' % (iconddir)
    quit(1)
#  project='ORCA'
  string='orcadir %s\nscratchdir %s/%s/QM\nncpu %i\nschedule_scaling %f\n' % (INFOS['orcadir'],INFOS['scratchdir'],iconddir,INFOS['orca.ncpu'],INFOS['orca.scaling'])
  string+='memory %i\n' % (INFOS['orca.mem'])
  if 'wfoverlap' in INFOS['needed']:
    string+='wfoverlap %s\nwfthres %f\n' % (INFOS['orca.wfoverlap'],INFOS['orca.ciothres'])
    #string+='numfrozcore %i\n' %(INFOS['frozcore_number'])
  else:
    string+='nooverlap\n'
  if 'theodore' in INFOS['needed']:
    string+='theodir %s\n' % (INFOS['orca.theodore'])
    string+='theodore_prop %s\n' % (INFOS['theodore.prop'])
    string+='theodore_fragment %s\n' % (INFOS['theodore.frag'])
  if 'tinker' in INFOS:
    string+='tinker %s\n' % (INFOS['tinker'])
  if 'ORCA.fffile' in INFOS:
    string+='qmmm_ff_file ORCA.qmmm.ff\n'
  if 'ORCA.ctfile' in INFOS:
    string+='qmmm_table ORCA.qmmm.table\n'
  sh2cas.write(string)
  sh2cas.close()

  # copy MOs and template
  cpfrom=INFOS['ORCA.template']
  cpto='%s/ORCA.template' % (iconddir)
  shutil.copy(cpfrom,cpto)

  if INFOS['orca.guess']:
    cpfrom1=INFOS['orca.guess']
    cpto1='%s/ORCA.gbw.init' % (iconddir)
    shutil.copy(cpfrom1,cpto1)

  if 'ORCA.fffile' in INFOS:
    cpfrom1=INFOS['ORCA.fffile']
    cpto1='%s/ORCA.qmmm.ff' % (iconddir)
    shutil.copy(cpfrom1,cpto1)

  if 'ORCA.ctfile' in INFOS:
    cpfrom1=INFOS['ORCA.ctfile']
    cpto1='%s/ORCA.qmmm.table' % (iconddir)
    shutil.copy(cpfrom1,cpto1)

  return



# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def checktemplate_BAGEL(filename,INFOS):
  necessary=['cobramm','basis','df_basis','nact','nclosed']
  try:
    f=open(filename)
    data=f.readlines()
    f.close()
  except IOError:
    print 'Could not open template file %s' % (filename)
    return False
  valid=[]
  for i in necessary:
    for l in data:
      if i in re.sub('#.*$','',l):
        valid.append(True)
        break
    else:
      valid.append(False)
  if not all(valid):
    print 'The template %s seems to be incomplete! It should contain: ' % (filename) +str(necessary)
    return False
  roots_there=False
  for l in data:
    l=re.sub('#.*$','',l).lower().split()
    if len(l)==0:
      continue
    if 'nstate' in l[0]:
      roots_there=True
  if not roots_there:
    for mult,state in enumerate(INFOS['states']):
      if state<=0:
        continue
      valid=[]
      for l in data:
        if 'spin' in re.sub('#.*$','',l).lower():
          f=l.split()
          if int(f[1])==mult+1:
            valid.append(True)
            break
      else:
        valid.append(False)
  if not all(valid):
    string='The template %s seems to be incomplete! It should contain the keyword "spin" for ' % (filename)
    for mult,state in enumerate(INFOS['states']):
      if state<=0:
        continue
      string+='%s, ' % (IToMult[mult+1])
    string=string[:-2]+'!'
    print string
    return False
  return True

# =================================================

def get_BAGEL(INFOS):
  '''This routine asks for all questions specific to BAGEL:
  - path to bagel
  - scratch directory
  - BAGEL.template
  - wf.init
  '''

  string='\n  '+'='*80+'\n'
  string+='||'+centerstring('BAGEL Interface setup',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  print string

  print centerstring('Path to BAGEL',60,'-')+'\n'
  path=os.getenv('BAGEL')
  #path=os.path.expanduser(os.path.expandvars(path))
  if path=='':
    path=None
  else:
    path='$BAGEL/'
      #print 'Environment variable $MOLCAS detected:\n$MOLCAS=%s\n' % (path)
      #if question('Do you want to use this MOLCAS installation?',bool,True):
        #INFOS['molcas']=path
    #if not 'molcas' in INFOS:
  print '\nPlease specify path to BAGEL directory (SHELL variables and ~ can be used, will be expanded when interface is started).\n'
  INFOS['bagel']=question('Path to BAGEL:',str,path)
  print ''


  print centerstring('BAGEL input template file',60,'-')+'\n'
  print '''Please specify the path to the BAGEL.template file. This file must contain the following settings:

basis <Basis set>
df_basis <Density fitting basis set>
nact <Number of active orbitals>
nclosed <Number of doubly occupied orbitals>
nstate <Number of states for state-averaging>

The BAGEL interface will generate the appropriate BAGEL input automatically.
'''
  if os.path.isfile('BAGEL.template'):
    if checktemplate_BAGEL('BAGEL.template',INFOS):
      print 'Valid file "BAGEL.template" detected. '
      usethisone=question('Use this template file?',bool,True)
      if usethisone:
        INFOS['bagel.template']='BAGEL.template'
  if not 'bagel.template' in INFOS:
    while True:
      filename=question('Template filename:',str)
      if not os.path.isfile(filename):
        print 'File %s does not exist!' % (filename)
        continue
      if checktemplate_BAGEL(filename,INFOS):
        break
    INFOS['molcas.template']=filename
  print ''

  print centerstring('Dipole level',60,'-')+'\n'
  print 'Please specify the desired amount of calculated dipole moments:\n0 -only dipole moments that are for free are calculated\n1 -calculate all transition dipole moments between the (singlet) ground state and all singlet states for absorption spectra\n2 -calculate all dipole moments'
  INFOS['dipolelevel']=question('Requested dipole level:',int,[1])[0]
  print ''





  print centerstring('Initial wavefunction: MO Guess',60,'-')+'\n'
  print '''Please specify the path to a MOLCAS JobIph file containing suitable starting MOs for the CASSCF calculation. Please note that this script cannot check whether the wavefunction file and the Input template are consistent!
'''
  INFOS['bagel.guess']={}
  string='Do you have initial wavefunction files for '
  for mult,state in enumerate(INFOS['states']):
    if state<=0:
      continue
    string+='%s, ' % (IToMult[mult+1])
  string=string[:-2]+'?'
  if question(string,bool,True):
    for mult,state in enumerate(INFOS['states']):
      if state<=0:
        continue
      while True:
        guess_file='archive.%i.init' % (mult+1)
        filename=question('Initial wavefunction file for %ss:' % (IToMult[mult+1]),str,guess_file)
        if os.path.isfile(filename):
          INFOS['bagel.guess'][mult+1]=filename
          break
        else:
          print 'File not found!'
  else:
    print 'WARNING: Remember that CASSCF calculations may run very long and/or yield wrong results without proper starting MOs.'
    time.sleep(1)

  print centerstring('BAGEL Ressource usage',60,'-')+'\n'#TODO

  print '''Please specify the number of CPUs to be used by EACH calculation.
'''
  INFOS['bagel.ncpu']=abs(question('Number of CPUs:',int,[1])[0])
  
  if INFOS['bagel.ncpu']>1:
    INFOS['bagel.mpi']=question('Use MPI mode (no=OpenMP)?',bool,False)
  else:
    INFOS['bagel.mpi']=False




  ## Ionization
  #need_wfoverlap=False
  #print centerstring('Ionization probability by Dyson norms',60,'-')+'\n'
  #INFOS['ion']=question('Dyson norms?',bool,False)
  #if 'ion' in INFOS and INFOS['ion']:
    #need_wfoverlap=True

  # wfoverlap
  if 'wfoverlap' in INFOS['needed']:
    print '\n'+centerstring('WFoverlap setup',60,'-')+'\n'
    INFOS['bagel.wfoverlap']=question('Path to wavefunction overlap executable:',str,'$SHARC/wfoverlap.x')
    # TODO not asked for: numfrozcore, numocc
    print '''Please specify the path to the PyQuante directory.
'''
    INFOS['bagel.pyq']=question('PyQuante path:',str)
    print '''Please specify the amount of memory available to wfoverlap.x (in MB). \n (Note that BAGEL's memory cannot be controlled)
'''
    INFOS['bagel.mem']=abs(question('wfoverlap.x memory:',int,[1000])[0])
  else:
    INFOS['bagel.mem']=1000
    INFOS['bagel.pyq']=''

  return INFOS

# =================================================

def prepare_BAGEL(INFOS,iconddir):
  # write BAGEL.resources
  try:
    sh2cas=open('%s/BAGEL.resources' % (iconddir), 'w')
  except IOError:
    print 'IOError during prepareBAGEL, iconddir=%s' % (iconddir)
    quit(1)
  project='BAGEL'
  string='bagel %s\npyquante %s\nscratchdir %s/%s/QM\nmemory %i\nncpu %i\ndipolelevel %i\nproject %s\n' % (INFOS['bagel'],INFOS['bagel.pyq'],INFOS['scratchdir'],iconddir,INFOS['bagel.mem'],INFOS['bagel.ncpu'],INFOS['dipolelevel'],project)

  if INFOS['bagel.mpi']:
    string+='mpi\n'
  if 'wfoverlap' in INFOS['needed']:
    string+='\nwfoverlap %s\n' % INFOS['bagel.wfoverlap']
  else:
    string+='\nnooverlap\n'
#  if 'tinker' in INFOS:
#    string+='tinker %s' % (INFOS['tinker'])
  sh2cas.write(string)
  sh2cas.close()

  # copy MOs and template
  cpfrom=INFOS['bagel.template']
  cpto='%s/BAGEL.template' % (iconddir)
  shutil.copy(cpfrom,cpto)
  if not INFOS['bagel.guess']=={}:
    for i in INFOS['bagel.guess']:
      cpfrom=INFOS['bagel.guess'][i]
      cpto='%s/%s.%i.init' % (iconddir,'archive',i)
      shutil.copy(cpfrom,cpto)



  return





# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def get_runscript_info(INFOS):
  ''''''

  string='\n  '+'='*80+'\n'
  string+='||'+centerstring('Run mode setup',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  print string

  print centerstring('Run script',60,'-')+'\n'
  print '''This script can generate the run scripts for each initial condition in two modes:

  - In mode 1, the calculation is run in subdirectories of the current directory.

  - In mode 2, the input files are transferred to another directory (e.g. a local scratch directory), the calculation is run there, results are copied back and the temporary directory is deleted. Note that this temporary directory is not the same as the "scratchdir" employed by the interfaces.

Note that in any case this script will create the input subdirectories in the current working directory. 
'''
  print 'In case of mode 1, the calculations will be run in:\n%s\n' % (INFOS['cwd'])
  here=question('Use mode 1 (i.e., calculate here)?',bool,True)
  if here:
    INFOS['here']=True
  else:
    INFOS['here']=False
    print '\nWhere do you want to perform the calculations? Note that this script cannot check whether the path is valid.'
    INFOS['copydir']=question('Run directory?',str)
  print ''

  print centerstring('Submission script',60,'-')+'\n'
  print '''During the setup, a script for running all initial conditions sequentially in batch mode is generated. Additionally, a queue submission script can be generated for all initial conditions.
'''
  qsub=question('Generate submission script?',bool,False)
  if not qsub:
    INFOS['qsub']=False
  else:
    INFOS['qsub']=True
    print '\nPlease enter a queue submission command, including possibly options to the queueing system,\ne.g. for SGE: "qsub -q queue.q -S /bin/bash -cwd" (Do not type quotes!).'
    INFOS['qsubcommand']=question('Submission command?',str,None,False)
    INFOS['proj']=question('Project Name:',str,None,False)

  print ''
  return INFOS

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def make_directory(iconddir):
  '''Creates a directory'''

  if os.path.isfile(iconddir):
    print '\nWARNING: %s is a file!' % (iconddir)
    return -1
  if os.path.isdir(iconddir):
    if len(os.listdir(iconddir))==0:
      return 0
    else:
      print '\nWARNING: %s/ is not empty!' % (iconddir)
      if not 'overwrite' in globals():
        global overwrite
        overwrite=question('Do you want to overwrite files in this and all following directories? ',bool,False)
      if overwrite:
        return 0
      else:
        return -1
  else:
    try:
      os.mkdir(iconddir)
    except OSError:
      print '\nWARNING: %s cannot be created!' % (iconddir)
      return -1
    return 0

# ======================================================================================================================

def writeQMin(INFOS,iconddir):
  icond=int(iconddir[-6:-1])
  try:
    qmin=open('%s/QM.in' % (iconddir), 'w')
  except IOError:
    print 'IOError during writeQMin, icond=%s' % (iconddir)
    quit(1)
  string='%i\nInitial condition %s\n' % (INFOS['natom'],iconddir)

  if icond>0:
    searchstring='Index\s+%i' % (icond)
  else:
    searchstring='Equilibrium'
  rewinded=False
  while True:
    try:
      line=INFOS['initf'].readline()
    except EOFError:
      if not rewinded:
        rewinded=True
        INFOS['initf'].seek(0)
      else:
        print 'Could not find Initial condition %i!' % (icond)
        quit(1)
    #if searchstring in line:
    if re.search(searchstring,line):
      break
  if icond>0:
    line=INFOS['initf'].readline()        # skip one line
  for iatom in range(INFOS['natom']):
    line=INFOS['initf'].readline()
    s=line.split()
    string+='%s %s %s %s\n' % (s[0],s[2],s[3],s[4])

  string+='unit bohr\nstates '
  for i in INFOS['states']:
    string+='%i ' % (i)
  string+='\n'
  
  string+='cobramm\n'
  if ('refov' in INFOS and INFOS['refov']):
    if icond==0:
      string+='init\nsavedir ./SAVEDIRQMMM/\n'
    else:
      currentdir=os.getcwd()
      iconddir='ICOND_%05i/' % (icond)
      string+='overlap\ncleanup\nsavedir %s\n' % (currentdir + '/' + iconddir + 'SAVEDIRQMMM')
  else:
    string+='init\ncleanup\n'

  if INFOS['soc']:
    string+='\nSOC\n'
  else:
    string+='\nH\n'
  string+='DM\n'
  if 'ion' in INFOS and INFOS['ion']:
    string+='ion\n'
  if 'theodore' in INFOS and INFOS['theodore']:
    string+='theodore\n'

  qmin.write(string)
  qmin.close()
  return

# ======================================================================================================================

def writeRunscript(INFOS,iconddir):
  '''writes the runscript in each subdirectory'''

  try:
    runscript=open('%s/run.sh' % (iconddir), 'w')
  except IOError:
    print 'IOError during writeRunscript, iconddir=%s' % (iconddir)
    quit(1)
  if 'proj' in INFOS:
    projname='%4s_%5s' % (INFOS['proj'][0:4],iconddir[-6:-1])
  else:
    projname='init_%5s' % (iconddir[-6:-1])

  # ================================
  intstring=''
  if 'adfrc' in INFOS:
    intstring='. %s\nexport PYTHONPATH=$ADFHOME/scripting:$PYTHONPATH' % (INFOS['adfrc'])

  # ================================
  if ('refov' in INFOS and INFOS['refov']) and iconddir!='ICOND_00000/':
    refstring='''
if [ -d ../ICOND_00000/SAVEDIRQMMM ];
then
  if [ -d ./SAVEDIRQMMM ];
  then
    rm -r ./SAVEDIRQMMM
  fi
  cp -r ../ICOND_00000/SAVEDIRQMMM ./SAVEDIRQMMM
else
  echo "Should do a reference overlap calculation, but the reference data in ../ICOND_00000/ seems not OK."
  exit 1
fi
'''
  else:
    refstring=''

  # generate run scripts here
  # ================================ for here mode
  if INFOS['here']:
    string='''#!/bin/bash

#$-N %s

%s

PRIMARY_DIR=%s/%s/

cd $PRIMARY_DIR
%s

$SHARC/SHARC_COBRAMM.py QM.in >> QMMM.log 2>> QMMM.err
''' % (projname,intstring,INFOS['cwd'], iconddir, refstring)
   
  # ================================ for remote mode
  else:
    string='''#!/bin/bash

#$-N %s

%s

PRIMARY_DIR=%s/%s/
COPY_DIR=%s/%s/

cd $PRIMARY_DIR
%s

mkdir -p $COPY_DIR
cp -r $PRIMARY_DIR/* $COPY_DIR
cd $COPY_DIR

$SHARC/SHARC_COBRAMM.py QM.in >> QMMM.log 2>> QMMM.err

cp -r $COPY_DIR/QM.* $COPY_DIR/SAVEDIRQMMM/ $PRIMARY_DIR
rm -r $COPY_DIR
''' % (projname,intstring,INFOS['cwd'], iconddir, INFOS['copydir'], iconddir, refstring)

  # ================================
  runscript.write(string)
  runscript.close()
  filename='%s/run.sh' % (iconddir)
  os.chmod(filename, os.stat(filename).st_mode | stat.S_IXUSR)
  return

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def setup_equilibrium(INFOS):
  #iconddir='ICOND_%05i/' % (0)
  #exists=os.path.isfile(iconddir+'/QM.out')
  exists=False
  if not exists:
    iconddir='ICOND_%05i/' % (0)
    io=make_directory(iconddir)
    if io!=0:
      print 'Skipping initial condition %s!' % (iconddir)
      return

    prepare_COBRAMM(INFOS,iconddir)
    writeQMin(INFOS,iconddir)
    globals()[Interfaces[ INFOS['interface']]['prepare_routine'] ](INFOS,iconddir)
    writeRunscript(INFOS,iconddir)
  return exists

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def setup_all(INFOS):
  '''This routine sets up the directories for the initial calculations.'''

  string='\n  '+'='*80+'\n'
  string+='||'+centerstring('Setting up directories...',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  print string

  all_run=open('all_run_init.sh','w')
  string='#/bin/bash\n\nCWD=%s\n\n' % (INFOS['cwd'])
  all_run.write(string)
  if INFOS['qsub']:
    all_qsub=open('all_qsub_init.sh','w')
    string='#/bin/bash\n\nCWD=%s\n\n' % (INFOS['cwd'])
    all_qsub.write(string)

  width=50
  ninit=INFOS['irange'][1]-INFOS['irange'][0]+1
  idone=0

  EqExists=setup_equilibrium(INFOS)
  if not EqExists:
    iconddir='ICOND_%05i/' % (0)
    string='cd $CWD/%s/\nbash run.sh\ncd $CWD\necho %s >> DONE\n' % (iconddir,iconddir)
    all_run.write(string)
    if INFOS['qsub']:
      string='cd $CWD/%s/\n%s run.sh\ncd $CWD\n' % (iconddir,INFOS['qsubcommand'])
      all_qsub.write(string)

  if INFOS['irange']!=[0,0]:
    for icond in range(INFOS['irange'][0],INFOS['irange'][1]+1):
      iconddir='ICOND_%05i/' % (icond)
      idone+=1
      done=idone*width/ninit
      sys.stdout.write('\rProgress: ['+'='*done+' '*(width-done)+'] %3i%%' % (done*100/width))
      sys.stdout.flush()

      io=make_directory(iconddir)
      if io!=0:
        print 'Skipping initial condition %s!' % (iconddir)
        continue
      prepare_COBRAMM(INFOS,iconddir)
      writeQMin(INFOS,iconddir)
      globals()[Interfaces[ INFOS['interface']]['prepare_routine'] ](INFOS,iconddir)
      writeRunscript(INFOS,iconddir)

      string='cd $CWD/%s/\nbash run.sh\ncd $CWD\necho %s >> DONE\n' % (iconddir,iconddir)
      all_run.write(string)
      if INFOS['qsub']:
        string='cd $CWD/%s/\n%s run.sh\ncd $CWD\n' % (iconddir,INFOS['qsubcommand'])
        all_qsub.write(string)

  all_run.close()
  filename='all_run_init.sh'
  os.chmod(filename, os.stat(filename).st_mode | stat.S_IXUSR)
  if INFOS['qsub']:
    all_qsub.close()
    filename='all_qsub_init.sh'
    os.chmod(filename, os.stat(filename).st_mode | stat.S_IXUSR)

  print '\n'


# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def main():
  '''Main routine'''

  usage='''
python setup_init.py

This interactive program prepares the initial excited-state calculations for SHARC.
As input it takes the initconds file, number of states and range of initconds.

Afterwards, it asks for the interface used and goes through the preparation depending on the interface.
'''

  description=''
  parser = OptionParser(usage=usage, description=description)

  displaywelcome()
  open_keystrokes()

  INFOS=get_general()
  INFOS=globals()[Interfaces[ INFOS['interface']]['get_routine'] ](INFOS)
  INFOS=get_runscript_info(INFOS)

  print '\n'+centerstring('Full input',60,'#')+'\n'
  for item in INFOS:
    print item, ' '*(25-len(item)), INFOS[item]
  print ''
  setup=question('Do you want to setup the specified calculations?',bool,True)
  print ''

  if setup:
    setup_all(INFOS)

  close_keystrokes()


# ======================================================================================================================

if __name__ == '__main__':
  try:
    main()
  except KeyboardInterrupt:
    print '\nCtrl+C makes me a sad SHARC ;-(\n'
    quit(0)

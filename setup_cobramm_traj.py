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

# Interactive script for the setup of dynamics calculations for SHARC
#
# usage: python setup_traj.py

import copy
import math
import sys
import re
import os
import stat
import shutil
import datetime
import random
from optparse import OptionParser
import readline
import time
from socket import gethostname
import ast

# =========================================================0
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

version='0.1'
versionneeded=[0.2, 1.0, 2.0, 2.1, float(version)]
versiondate=datetime.date(2020,11,11)


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

# ======================================================================= #

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
# 4: {'script':          'SHARC_ADF.py',
#     'name':            'adf',
#     'description':     'ADF (DFT, TD-DFT)',
#     'get_routine':     'get_ADF',
#     'prepare_routine': 'prepare_ADF',
#     'features':        {'overlap': ['wfoverlap'],
#                         'dyson':   ['wfoverlap'],
#                         'theodore':['theodore'],
#                         'phases':  ['wfoverlap'],
#                         'soc':     []                 },
#     'pysharc':          False
#    },
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


Couplings={
  1: {'name':        'nacdt',
      'description': 'DDT     =  < a|d/dt|b >        Hammes-Schiffer-Tully scheme   '
     },
  2: {'name':        'nacdr',
      'description': 'DDR     =  < a|d/dR|b >        Original Tully scheme          '
     },
  3: {'name':        'overlap',
      'description': 'overlap = < a(t0)|b(t) >       Local Diabatization scheme     '
     }
  }

EkinCorrect={
  1: {'name':             'none',
      'description':      'Do not conserve total energy. Hops are never frustrated.',
      'description_refl': 'Do not reflect at a frustrated hop.',
      'required':   []
     },
  2: {'name':             'parallel_vel',
      'description':      'Adjust kinetic energy by rescaling the velocity vectors. Often sufficient.',
      'description_refl': 'Reflect the full velocity vector.',
      'required':   []
     },
  3: {'name':             'parallel_nac',
      'description':      'Adjust kinetic energy only with the component of the velocity vector along the non-adiabatic coupling vector.',
      'description_refl': 'Reflect only the component of the velocity vector along the non-adiabatic coupling vector.',
      'required':   ['nacdr']
     },
  4: {'name':             'parallel_diff',
      'description':      'Adjust kinetic energy only with the component of the velocity vector along the gradient difference vector.',
      'description_refl': 'Reflect only the component of the velocity vector along the gradient difference vector.',
      'required':   []
     }
  }

Decoherences={
  1: {'name':             'none',
      'description':      'No decoherence correction.',
      'required':   [],
      'params':     ''
     },
  2: {'name':             'edc',
      'description':      'Energy-based decoherence scheme (Granucci, Persico, Zoccante).',
      'required':   [],
      'params':     '0.1'
     },
  3: {'name':             'afssh',
      'description':      'Augmented fewest-switching surface hopping (Jain, Alguire, Subotnik).',
      'required':   [],
      'params':     ''
     }
  }

HoppingSchemes={
  1: {'name':             'off',
      'description':      'Surface hops off.'
     },
  2: {'name':             'sharc',
      'description':      'Standard SHARC surface hopping probabilities (Mai, Marquetand, Gonzalez).'
     },
  3: {'name':             'gfsh',
      'description':      'Global flux surface hopping probabilities (Wang, Trivedi, Prezhdo).'
     }
  }

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def try_read(l,index,typefunc,default):
  try:
    if typefunc==bool:
      return 'True'==l[index]
    else:
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
    self.dip     = [ complex( try_read(f,i,float,0.),try_read(f,i+1,float,0.) ) for i in [3,5,7] ]
    self.Excited =   try_read(f,11,bool, False)
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
    s+='% 12.8f % 12.8f %s' % (self.Eexc*HARTREE_TO_EV,self.Fosc,self.Excited)
    return s

  def Excite(self,max_Prob,erange):
    try:
      Prob=self.Prob/max_Prob
    except ZeroDivisionError:
      Prob=-1.
    if not (erange[0] <= self.Eexc <= erange[1]):
      Prob=-1.
    self.Excited=(random.random() < Prob)

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

def displaywelcome():
  print 'Script for setup of SHARC trajectories started...\n'
  string='\n'
  string+='  '+'='*80+'\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Setup QM/MM trajectories for SHARC/COBRAMM dynamics',80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Author: Davide Avagliano',80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='||'+centerstring('Version:'+version,80)+'||\n'
  string+='||'+centerstring(versiondate.strftime("%d.%m.%y"),80)+'||\n'
  string+='||'+centerstring('',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  string+='''
This script automatizes the setup of the input files for SHARC dynamics.
  '''
  print string

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def open_keystrokes():
  global KEYSTROKES
  KEYSTROKES=open('KEYSTROKES.tmp','w')

def close_keystrokes():
  KEYSTROKES.close()
  shutil.move('KEYSTROKES.tmp','KEYSTROKES.setup_cobramm_traj')

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

def itnmstates(states):
  for i in range(len(states)):
    if states[i]<1:
      continue
    for k in range(i+1):
      for j in range(states[i]):
        yield i+1,j+1,k-i/2.
  return

# ======================================================================================================================


class init_string:
  def __init__(self):
    self.strings=[]
    self.nst=0
    self.width=100
    self.group=10
    self.groups=(self.width-1)/self.group+1
    self.nrow=1
    self.lastrow=0
  def add(self,s):
    self.strings.append(s)
    self.nst+=1
    self.nrow=(self.nst-1)/self.width+1
    self.lastrow=self.nst%self.width
    if self.lastrow==0:
      self.lastrow=self.width
  def reset(self):
    self.strings=[]
    self.nst=0
    self.nrow=1
    self.lastrow=0
  def __str__(self):
    nw=int(math.log(self.nst)/math.log(10)+1.1)
    s=' '*(nw+2)
    fs='%%%ii' % (nw)
    for i in range(self.groups):
      s+=' '*(self.group-nw+1)+fs % ((i+1)*self.group)
    s+='\n'
    s+=' '*(nw+2)
    for i in range(self.groups):
      s+=' '
      for j in range(self.group-1):
        s+=' '
      s+='|'
    s+='\n'
    index=0
    for i in range(self.nrow):
      s+=fs % (i*self.width) + ' | '
      for j in range(self.width):
        try:
          s+=self.strings[index]
        except IndexError:
          return s
        index+=1
        if (j+1)%self.group==0:
          s+=' '
      s+='\n'
    s+='\n'
    return s

# ======================================================================================================================

def analyze_initconds(initlist,INFOS):
  if INFOS['show_content']:
    print 'Contents of the initconds file:'
    print '''\nLegend:
?       Geometry and Velocity
.       not selected
#       selected
'''
  n_hasexc=[]
  n_issel=[]
  display=init_string()
  for state in range(INFOS['nstates']):
    if INFOS['show_content']:
      print 'State %i:' % (state+1)
    display.reset()
    n_hasexc.append(0)
    n_issel.append(0)
    for i in initlist:
      if len(i.statelist)<state+1:
        display.add('?')
      else:
        n_hasexc[-1]+=1
        if i.statelist[state].Excited:
          display.add('#')
          n_issel[-1]+=1
        else:
          display.add('.')
    if INFOS['show_content']:
      print display
  print 'Number of excited states and selections:'
  print   'State    #InitCalc       #Selected'
  for i in range(len(n_hasexc)):
    s= '% 5i        % 5i           % 5i' % (i+1,n_hasexc[i],n_issel[i])
    if not INFOS['isactive'][i]:
      s+='  inactive'
    print s
  return n_issel

# ======================================================================================================================

def get_initconds(INFOS):
  ''''''

  INFOS['initf'].seek(0)                 # rewind the initf file
  initlist=[]
  for icond in range(1,INFOS['ninit']+1):
    initcond=INITCOND()
    initcond.init_from_file(INFOS['initf'],INFOS['eref'],icond)
    initlist.append(initcond)
  print 'Number of initial conditions in file:       %5i' % (INFOS['ninit'])

  INFOS['initlist']=initlist
  INFOS['n_issel']=analyze_initconds(initlist,INFOS)
  return INFOS

# ======================================================================================================================

def check_laserfile(filename,nsteps,dt):
  try:
    f=open(filename)
    data=f.readlines()
    f.close()
  except IOError:
    print 'Could not open laser file %s' % (filename)
    return False
  n=0
  for line in data:
    if len(line.split())>=8:
      n+=1
    else:
      break
  if n<nsteps:
    print 'File %s has only %i timesteps, %i steps needed!' % (filename,n,nsteps)
    return False
  for i in range(int(nsteps)-1):
    t0=float(data[i].split()[0])
    t1=float(data[i+1].split()[0])
    if abs(abs(t1-t0)-dt)>1e-6:
      print 'Time step wrong in file %s at line %i.' % (filename,i+1)
      return False
  return True

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

  string='\n  '+'='*80+'\n'
  string+='||'+centerstring('Initial conditions',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  print string
  print '''\nThis script reads the initial conditions (geometries, velocities, initial excited state)
from the initconds.excited files as provided by excite.py.
'''

  # open the initconds file
  try:
    initfile='initconds.excited'
    initf=open(initfile)
    line=initf.readline()
    if check_initcond_version(line,must_be_excited=True):
      print 'Initial conditions file "initconds.excited" detected. Do you want to use this?'
      if not question('Use file "initconds.excited"?',bool,True):
        initf.close()
        raise IOError
    else:
      initf.close()
      raise IOError
  except IOError:
    print 'Please enter the filename of the initial conditions file.'
    while True:
      initfile=question('Initial conditions filename:',str,'initconds.excited')
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
      if check_initcond_version(line,must_be_excited=True):
        break
      else:
        print 'File does not contain initial conditions!'
        continue
  # read the header
  INFOS['ninit']=int(initf.readline().split()[1])
  INFOS['natom']=int(initf.readline().split()[1])
  print '\nFile %s contains %i initial conditions.' % (initfile,INFOS['ninit'])
  print 'Number of atoms is %i' % (INFOS['natom'])
  INFOS['repr']=initf.readline().split()[1]
  if INFOS['repr']=='MCH':
    INFOS['diag']=False
  else:
    INFOS['diag']=True
  INFOS['eref']=float(initf.readline().split()[1])
  INFOS['eharm']=float(initf.readline().split()[1])

  # get guess for number of states
  line=initf.readline()
  if 'states' in line.lower():
    states=[]
    l=line.split()
    for i in range(1,len(l)):
      states.append(int(l[i]))
    guessstates=states
  else:
    guessstates=None

  print 'Reference energy %16.12f a.u.' % (INFOS['eref'])
  print 'Excited states are in %s representation.\n' % (['MCH','diagonal'][INFOS['diag']])
  initf.seek(0)                 # rewind the initf file
  INFOS['initf']=initf

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

#rattle

  print '\n\n'+centerstring('Rattle file',60,'-')+'\n'
  INFOS['rattle']=question('Do you want to include rattle constraint in the simulation?',bool,False)
  if INFOS['rattle']:
    print '''Please specify the file containing the restrain information (#index1 #index2 #distance)
'''
    if os.path.isfile('rattlefile'):
      #if check_rattlefile('rattlefile'):
        print 'Valid rattle file "rattlefile" detected. '
        usethisrattle=question('Use this rattle file?',bool,True)
        if usethisrattle:
          INFOS['rattlefile']='rattlefile'
    if not 'rattlefile' in INFOS:
      while True:
        rattle_filename=question('rattle filename:',str)
	try:
	  lf=open(rattle_filename)
	except IOError: 
	  print 'Could not open: %s' % (path4)
	  continue
	break
      INFOS['rattlefile']=rattle_filename
      print 'Rattle constraints are read in %s\n' % (rattle_filename)

  # Number of states
  print '\nPlease enter the number of states as a list of integers\ne.g. 3 0 3 for three singlets, zero doublets and three triplets.'
  while True:
    states=question('Number of states:',int,guessstates)
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
  # obtain the statemap
  statemap={}
  i=1
  for imult,istate,ims in itnmstates(INFOS['states']):
    statemap[i]=[imult,istate,ims]
    i+=1
  INFOS['statemap']=statemap

  # get active states
  if question('Do you want all states to be active?',bool,True):
    INFOS['actstates']=INFOS['states']
  else:
    print '\nPlease enter the number of ACTIVE states as a list of integers\ne.g. 3 0 3 for three singlets, zero doublets and three triplets.'
    while True:
      actstates=question('Number of states:',int)
      if len(actstates)!=len(INFOS['states']):
        print 'Length of nstates and actstates must match!'
        continue
      valid=True
      for i,nst in enumerate(actstates):
        if not 0<=nst<=INFOS['states'][i]:
          print 'Number of active states of multiplicity %i must be between 0 and the number of states of this multiplicity (%i)!' % (i+1,INFOS['states'][i])
          valid=False
      if not valid:
        continue
      break
    INFOS['actstates']=actstates
  isactive=[]
  for imult in range(len(INFOS['states'])):
    for ims in range(imult+1):
      for istate in range(INFOS['states'][imult]):
        isactive.append( (istate+1<=INFOS['actstates'][imult]) )
  INFOS['isactive']=isactive
  print ''


  # ask whether initfile content is shown
  INFOS['show_content']=question('Do you want to see the content of the initconds file?',bool,False)



  # read initlist, analyze it and print content (all in get_initconds)
  INFOS['initf']=initf
  INFOS=get_initconds(INFOS)


  # Generate random example for setup-states, according to Leti's wishes
  exampleset=set()
  nactive=sum(INFOS['isactive'])
  while len(exampleset)<min(3,nactive):
    i=random.randint(1,INFOS['nstates'])
    if INFOS['isactive'][i-1]:
      exampleset.add(i)
  exampleset=list(exampleset)
  exampleset.sort()
  string1=''
  string2=''
  j=0
  for i in exampleset:
    j+=1
    if j==len(exampleset) and len(exampleset)>1:
      string1+=str(i)
      string2+='and '+str(i)
    else:
      string1+=str(i)+' '
      string2+=str(i)+', '



  # ask for states to setup
  print '\nPlease enter a list specifying for which excited states trajectories should be set-up\ne.g. %s to select states %s.' % (string1,string2)
  defsetupstates=[]
  nmax=0
  for i,active in enumerate(INFOS['isactive']):
    if active and INFOS['n_issel'][i]>0:
      defsetupstates.append(i+1)
      nmax+=INFOS['n_issel'][i]
  if nmax<=0:
    print '\nZero trajectories can be set up!'
    sys.exit(1)
  while True:
    setupstates=question('States to setup the dynamics:',int,defsetupstates,ranges=True)
    valid=True
    for i in setupstates:
      if i>INFOS['nstates']:
        print 'There are only %i states!' % (INFOS['nstates'])
        valid=False
        continue
      if i<0:
        valid=False
        continue
      if not INFOS['isactive'][i-1]:
        print 'State %i is inactive!' % (i)
        valid=False
    if not valid:
      continue
    INFOS['setupstates']=set(setupstates)
    nsetupable=sum( [ INFOS['n_issel'][i-1] for i in INFOS['setupstates'] if INFOS['isactive'][i-1] ] )
    print '\nThere can be %i trajector%s set up.\n' % (nsetupable,['y','ies'][nsetupable!=1])
    if nsetupable==0:
      continue
    break


  # select range within initconds file
  # only start index needed, end index is determined by number of trajectories
  print 'Please enter the index of the first initial condition in the initconds file to be setup.'
  while True:
    firstindex=question('Starting index:',int,[1])[0]
    if not 0<firstindex<=INFOS['ninit']:
      print 'Please enter an integer between %i and %i.' % (1,INFOS['ninit'])
      continue
    nsetupable=0
    for i,initcond in enumerate(INFOS['initlist']):
      if i+1<firstindex:
        continue
      for state in set(setupstates):
        try:
          nsetupable+=initcond.statelist[state-1].Excited
        except IndexError:
          break
    print '\nThere can be %i trajector%s set up, starting in %i states.' % (nsetupable,['y','ies'][nsetupable!=1],len(INFOS['setupstates']))
    if nsetupable==0:
      continue
    break
  INFOS['firstindex']=firstindex


  # Number of trajectories
  print '\nPlease enter the total number of trajectories to setup.'
  while True:
    ntraj=question('Number of trajectories:',int,[nsetupable])[0]
    if not 1<=ntraj<=nsetupable:
      print 'Please enter an integer between %i and %i.' % (1,nsetupable)
      continue
    break
  INFOS['ntraj']=ntraj


  # Random number seed
  print '\nPlease enter a random number generator seed (type "!" to initialize the RNG from the system time).'
  while True:
    line=question('RNG Seed: ',str,'!',False)
    if line=='!':
      random.seed()
      break
    try:
      rngseed=int(line)
      random.seed(rngseed)
    except ValueError:
      print 'Please enter an integer or "!".'
      continue
    break
  print ''

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


  # Interface
  string='\n  '+'='*80+'\n'
  string+='||'+centerstring('Choose the quantum chemistry interface',80)+'||\n'
  string+='  '+'='*80+'\n'
  print string
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



  # Dynamics options
  string='\n  '+'='*80+'\n'
  string+='||'+centerstring('Surface Hopping dynamics settings',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  print string


  # Simulation time
  print centerstring('Simulation time',60,'-')+'\n'
  print 'Please enter the total simulation time.'
  while True:
    num2=question('Simulation time (fs):',float,[1000.])[0]
    if num2<=0:
      print 'Simulation time must be positive!'
      continue
    break
  INFOS['tmax']=num2


  # Timestep
  print '\nPlease enter the simulation timestep (0.5 fs recommended).'
  while True:
    dt=question('Simulation timestep (fs):',float,[0.5])[0]
    if dt<=0:
      print 'Simulation timestep must be positive!'
      continue
    break
  INFOS['dtstep']=dt
  print '\nSimulation will have %i timesteps.' % (num2/dt+1)


  # number of substeps
  print '\nPlease enter the number of substeps for propagation (25 recommended).'
  while True:
    nsubstep=question('Nsubsteps:',int,[25])[0]
    if nsubstep<=0:
      print 'Enter a positive integer!'
      continue
    break
  INFOS['nsubstep']=nsubstep


  # whether to kill relaxed trajectories
  print '\nThe trajectories can be prematurely terminated after they run for a certain time in the lowest state. '
  INFOS['kill']=question('Do you want to prematurely terminate trajectories?',bool,False)
  if INFOS['kill']:
    while True:
      tkill=question('Kill after (fs):',float,[10.])[0]
      if tkill<=0:
        print 'Must be positive!'
        continue
      break
    INFOS['killafter']=tkill
  print ''


  print '\n'+centerstring('Dynamics settings',60,'-')


  # SHARC or FISH
  print '\nDo you want to perform the dynamics in the diagonal representation (SHARC dynamics) or in the MCH representation (regular surface hopping)?'
  surf=question('SHARC dynamics?',bool,True)
  INFOS['surf']=['mch','diagonal'][surf]

  ## SOC or not
  #recommended=True
  #if len(INFOS['states'])==1:
    #recommended=False
  #print '\nDo you want to include spin-orbit couplings in the dynamics?'
  #INFOS['soc']=question('Spin-orbit couplings?',bool,recommended)

  # Setup SOCs
  if len(states)>1:
    if 'soc' in Interfaces[INFOS['interface']]['features']:
      print 'Do you want to include spin-orbit couplings in the dynamics?\n'
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
  INFOS['states']=states
  INFOS['nstates']=nstates
  INFOS['soc']=soc
  if INFOS['soc']:
    INFOS['needed'].extend(Interfaces[INFOS['interface']]['features']['soc'])

  # Coupling
  print '\nPlease choose the quantities to describe non-adiabatic effects between the states:'
  for i in Couplings:
    print '%i\t%s%s' % (i,
                        Couplings[i]['description'],
                        ['(not available)',''][Couplings[i]['name'] in Interfaces[INFOS['interface']]['features']]
                        )
  #print ''
  while True:
    default=None
    for i in Couplings:
      if Couplings[i]['name'] in Interfaces[INFOS['interface']]['features']:
        default=[i]
    #if len(Interfaces[INFOS['interface']]['couplings'])==1:
      #default=Interfaces[INFOS['interface']]['couplings']
    #else:
      #default=None
    num=question('Coupling number:',int,default)[0]
    if num in Couplings and Couplings[num]['name'] in Interfaces[INFOS['interface']]['features']:
      break
    else:
      l=[]
      for i in Couplings:
        if Couplings[i]['name'] in Interfaces[INFOS['interface']]['features']:
          l.append(i)
      print 'Please input one of the following: %s!' % (l)
  INFOS['coupling']=num
  INFOS['needed'].extend(Interfaces[INFOS['interface']]['features'][Couplings[i]['name']])


  # Phase tracking
  INFOS['phases_from_interface']=False
  if Couplings[INFOS['coupling']]['name']!='overlap':
    if 'phases' in Interfaces[INFOS['interface']]['features']:
      INFOS['phases_from_interface']=question('Do you want to track wavefunction phases through overlaps?',bool,True)
      if INFOS['phases_from_interface']:
        INFOS['needed'].extend(Interfaces[INFOS['interface']]['features']['phases'])

  # Gradient correction (only for SHARC)
  if INFOS['surf']=='diagonal':
    possible= ('nacdr' in Interfaces[INFOS['interface']]['features'])
    recommended=Couplings[INFOS['coupling']]['name']=='nacdr'
    print '\nFor SHARC dynamics, the evaluation of the mixed gradients necessitates to calculate non-adiabatic coupling vectors %s.' % (['(Extra computational cost)',' (Recommended)'][recommended])
    if possible:
      #while True:
        INFOS['gradcorrect']=question('Include non-adiabatic couplings in the gradient transformation?',bool,recommended)
        #if INFOS['gradcorrect'] and not 'nacdr' in Interfaces[INFOS['interface']]['features']:
          #print 'Not possible with the chosen interface!'
        #else:
          #break
    else:
      print '... but interface cannot provide non-adiabatic coupling vectors, turning option off.'
      INFOS['gradcorrect']=False
  else:
    INFOS['gradcorrect']=False
  if INFOS['gradcorrect']:
    INFOS['needed'].extend(Interfaces[INFOS['interface']]['features']['nacdr'])


  # Kinetic energy modification
  print '\nDuring a surface hop, the kinetic energy has to be modified in order to conserve total energy. There are several options to that:'
  cando=[]
  for i in EkinCorrect:
    recommended=len(EkinCorrect[i]['required'])==0  or  Couplings[INFOS['coupling']]['name'] in EkinCorrect[i]['required']
    possible= all([ j in Interfaces[INFOS['interface']]['features']  for j in EkinCorrect[i]['required']])
    if possible:
      cando.append(i)
    if not possible:
      print '%i\t%s%s' % (i, EkinCorrect[i]['description'],'\n\t(not possible)' )
    else:
      print '%i\t%s%s' % (i, EkinCorrect[i]['description'],['\n\t(extra computational cost)',''][ recommended ])
  while True:
    ekinc=question('EkinCorrect:',int,[2])[0]
    if ekinc in EkinCorrect and ekinc in cando:
      break
    else:
      print 'Please input one of the following: %s!' % ([i for i in cando])
  INFOS['ekincorrect']=ekinc
  if INFOS['ekincorrect']:
    for i in EkinCorrect[INFOS['ekincorrect']]['required']:
      INFOS['needed'].extend(Interfaces[INFOS['interface']]['features'][i])


  # frustrated reflection
  print '\nIf a surface hop is refused (frustrated) due to insufficient energy, the velocity can either be left unchanged or reflected:'
  cando=[]
  for i in EkinCorrect:
    recommended=len(EkinCorrect[i]['required'])==0  or  Couplings[INFOS['coupling']]['name'] in EkinCorrect[i]['required']
    possible= all([ j in Interfaces[INFOS['interface']]['features']  for j in EkinCorrect[i]['required']])
    if possible:
      cando.append(i)
    if not possible:
      print '%i\t%s%s' % (i, EkinCorrect[i]['description_refl'],'\n\t(not possible)' )
    else:
      print '%i\t%s%s' % (i, EkinCorrect[i]['description_refl'],['\n\t(extra computational cost)',''][ recommended ])
  while True:
    reflect=question('Reflect frustrated:',int,[1])[0]
    if reflect in EkinCorrect and reflect in cando:
      break
    else:
      print 'Please input one of the following: %s!' % ([i for i in cando])
  INFOS['reflect']=reflect
  if INFOS['reflect']:
    for i in EkinCorrect[INFOS['ekincorrect']]['required']:
      INFOS['needed'].extend(Interfaces[INFOS['interface']]['features'][i])


  # decoherence
  print '\nPlease choose a decoherence correction for the %s states:' % (['MCH','diagonal'][INFOS['surf']=='diagonal'])
  cando=[]
  for i in Decoherences:
    recommended=len(Decoherences[i]['required'])==0  or  Couplings[INFOS['coupling']]['name'] in Decoherences[i]['required']
    possible= all([ j in Interfaces[INFOS['interface']]['features']  for j in Decoherences[i]['required']])
    if possible:
      cando.append(i)
    if not possible:
      print '%i\t%s%s' % (i, Decoherences[i]['description'],'\n\t(not possible)' )
    else:
      print '%i\t%s%s' % (i, Decoherences[i]['description'],['\n\t(extra computational cost)',''][ recommended ])
  while True:
    decoh=question('Decoherence scheme:',int,[2])[0]
    if decoh in Decoherences and decoh in cando:
      break
    else:
      print 'Please input one of the following: %s!' % ([i for i in cando])
  INFOS['decoherence']=[Decoherences[decoh]['name'],Decoherences[decoh]['params']]
  for i in Decoherences[decoh]['required']:
    INFOS['needed'].extend(Interfaces[INFOS['interface']]['features'][i])


  # surface hopping scheme
  print '\nPlease choose a surface hopping scheme for the %s states:' % (['MCH','diagonal'][INFOS['surf']=='diagonal'])
  cando=list(HoppingSchemes)
  for i in HoppingSchemes:
    print '%i\t%s' % (i, HoppingSchemes[i]['description'])
  while True:
    hopping=question('Hopping scheme:',int,[2])[0]
    if hopping in HoppingSchemes and hopping in cando:
      break
    else:
      print 'Please input one of the following: %s!' % ([i for i in cando])
  INFOS['hopping']=HoppingSchemes[hopping]['name']

  # Forced hops to lowest state
  print '\nDo you want to perform forced hops to the lowest state based on a energy gap criterion?'
  print '(Note that this ignores spin multiplicity)'
  INFOS['force_hops']=question('Forced hops to ground state?',bool, False)
  if INFOS['force_hops']:
    INFOS['force_hops_dE']=abs( question('Energy gap threshold for forced hops (eV):',float,[0.1])[0] )
  else:
    INFOS['force_hops_dE']=9999.

  # Scaling
  print '\nDo you want to scale the energies and gradients?'
  scal=question('Scaling?',bool,False)
  if scal:
    while True:
      fscal=question('Scaling factor (>0.0): ',float)[0]
      if fscal<=0:
        print 'Please enter a positive real number!'
        continue
      break
    INFOS['scaling']=fscal
  else:
    INFOS['scaling']=False


  # Damping
  print '\nDo you want to damp the dynamics (Kinetic energy is reduced at each timestep by a factor)?'
  damp=question('Damping?',bool,False)
  if damp:
    while True:
      fdamp=question('Scaling factor (0-1): ',float)[0]
      if not 0<=fdamp<=1:
        print 'Please enter a real number 0<=r<=1!'
        continue
      break
    INFOS['damping']=fdamp
  else:
    INFOS['damping']=False


  # atommask
  INFOS['atommaskarray']=[]
  if (INFOS['decoherence'][0]=='edc') or (INFOS['ekincorrect']==2) or (INFOS['reflect']==2):
    print '\nDo you want to use an atom mask for velocity rescaling or decoherence?'
    if question('Atom masking?',bool,True):
      print '\nPlease enter all atom indices (start counting at 1) of the atoms which should be included in the velocity rescaling or decoherence. \nRemember that you can also enter ranges (e.g., "-1~-3  5  11~21").'
      arr=question('Masked atoms:',int,ranges=True)
      for i in arr:
        if 1<=i<=INFOS['natom']:
          INFOS['atommaskarray'].append(i)

  # selection of gradients (only for SHARC) and NACs (only if NAC=ddr)
  print '\n'+centerstring('Selection of Gradients and NACs',60,'-')+'\n'
  print '''In order to speed up calculations, SHARC is able to select which gradients and NAC vectors it has to calculate at a certain timestep. The selection is based on the energy difference between the state under consideration and the classical occupied state.
'''
  if INFOS['surf']=='diagonal':
    if INFOS['soc']:
      sel_g=question('Select gradients?',bool,False)
    else:
      sel_g=True
  else:
    sel_g=False
  INFOS['sel_g']=sel_g
  if Couplings[INFOS['coupling']]['name']=='ddr' or INFOS['gradcorrect'] or EkinCorrect[INFOS['ekincorrect']]['name']=='parallel_nac':
    sel_t=question('Select non-adiabatic couplings?',bool,False)
  else:
    sel_t=False
  INFOS['sel_t']=sel_t
  if sel_g or sel_t:
    if not sel_t and not INFOS['soc']:
      INFOS['eselect']=0.001
      print '\nSHARC dynamics without SOC and NAC: setting minimal selection threshold.'
    else:
      print '\nPlease enter the energy difference threshold for the selection of gradients and non-adiabatic couplings (in eV). (0.5 eV recommended, or even larger if SOC is strong in this system.)'
      eselect=question('Selection threshold (eV):',float,[0.5])[0]
      INFOS['eselect']=abs(eselect)


  # Laser file
  print '\n\n'+centerstring('Laser file',60,'-')+'\n'
  INFOS['laser']=question('Do you want to include a laser field in the simulation?',bool,False)
  if INFOS['laser']:
    print '''Please specify the file containing the complete laser field. The timestep in the file and the length of the file must fit to the simulation time, time step and number of substeps given above.

Laser files can be created using $SHARC/laser.x
'''
    if os.path.isfile('laser'):
      if check_laserfile('laser',INFOS['tmax']/INFOS['dtstep']*INFOS['nsubstep']+1,INFOS['dtstep']/INFOS['nsubstep']):
        print 'Valid laser file "laser" detected. '
        usethisone=question('Use this laser file?',bool,True)
        if usethisone:
          INFOS['laserfile']='laser'
    if not 'laserfile' in INFOS:
      while True:
        filename=question('Laser filename:',str)
        if not os.path.isfile(filename):
          print 'File %s does not exist!' % (filename)
          continue
        if check_laserfile(filename,INFOS['tmax']/INFOS['dtstep']*INFOS['nsubstep']+1,INFOS['dtstep']/INFOS['nsubstep']):
          break
      INFOS['laserfile']=filename
    # only the analytical interface can do dipole gradients
    if 'dipolegrad' in Interfaces[INFOS['interface']]['features']:
      INFOS['dipolegrad']=question('Do you want to use dipole moment gradients?',bool,False)
    else:
      INFOS['dipolegrad']=False
    print ''
  else:
    INFOS['dipolegrad']=False
  if INFOS['dipolegrad']:
    INFOS['needed'].extend(Interfaces[INFOS['interface']]['features']['dipolegrad'])



  # Setup Dyson computation
  INFOS['ion']=False
  if 'dyson' in Interfaces[INFOS['interface']]['features']:
    n=[0,0]
    for i,j in enumerate(INFOS['states']):
      n[i%2]+=j
    if n[0]>=1 and n[1]>=1:
      print '\n'+centerstring('Ionization probability by Dyson norms',60,'-')+'\n'
      print 'Do you want to compute Dyson norms between neutral and ionic states?'
      INFOS['ion']=question('Dyson norms?',bool,False)
      if INFOS['ion']:
        INFOS['needed'].extend(Interfaces[INFOS['interface']]['features']['dyson'])


  # Setup theodore
  if 'theodore' in Interfaces[INFOS['interface']]['features']:
    print '\n'+centerstring('TheoDORE wave function analysis',60,'-')+'\n'
    print 'Do you want to run TheoDORE to obtain one-electron descriptors for the electronic wave functions?'
    INFOS['theodore']=question('TheoDORE?',bool,False)
    if INFOS['theodore']:
      INFOS['needed'].extend(Interfaces[INFOS['interface']]['features']['theodore'])






  #print INFOS['needed']


  # Interface-specific section
  INFOS=globals()[Interfaces[ INFOS['interface']]['get_routine'] ](INFOS)



  # PYSHARC
  if Interfaces[ INFOS['interface']]['pysharc']:
    string='\n  '+'='*80+'\n'
    string+='||'+centerstring('PYSHARC',80)+'||\n'
    string+='  '+'='*80+'\n'
    print string
    print '\nThe chosen interface can be run very efficiently with PYSHARC.'
    print 'PYSHARC runs the SHARC dynamics directly within Python (with C and Fortran extension)'
    print 'with minimal file I/O for maximum performance.'
    INFOS['pysharc']=question('Setup for PYSHARC?',bool,True)
  else:
    INFOS['pysharc']=False

  # Dynamics options
  string='\n  '+'='*80+'\n'
  string+='||'+centerstring('Content of output.dat files',80)+'||\n'
  string+='  '+'='*80+'\n'
  print string

  # NetCDF
  print '\nSHARC or PYSHARC can produce output in ASCII format (all features supported currently)'
  print 'or in NetCDF format (more efficient file I/O, some features currently not supported).'
  INFOS['netcdf']=question('Write output in NetCDF format?',bool,INFOS['pysharc'])


  # options for writing to output.dat
  print '\nDo you want to write the gradients to the output.dat file ?'
  write_grad=question('Write gradients?',bool,False)
  if write_grad:
    INFOS['write_grad']=True
  else:
    INFOS['write_grad']=False

  print '\nDo you want to write the non-adiabatic couplings (NACs) to the output.dat file ?'
  write_NAC=question('Write NACs?',bool,False)
  if write_NAC:
    INFOS['write_NAC']=True
  else:
    INFOS['write_NAC']=False


  print '\nDo you want to write property matrices to the output.dat file  (e.g., Dyson norms)?'
  if 'ion' in INFOS and INFOS['ion']:
    INFOS['write_property2d']=question('Write property matrices?',bool,True)
  else:
    INFOS['write_property2d']=question('Write property matrices?',bool,False)


  print '\nDo you want to write property vectors to the output.dat file  (e.g., TheoDORE results)?'
  if 'theodore' in INFOS and INFOS['theodore']:
    INFOS['write_property1d']=question('Write property vectors?',bool,True)
  else:
    INFOS['write_property1d']=question('Write property vectors?',bool,False)


  print '\nDo you want to write the overlap matrix to the output.dat file ?'
  INFOS['write_overlap']=question('Write overlap matrix?',bool, (Couplings[INFOS['coupling']]['name']=='overlap') )


  print '\nDo you want to modify the output.dat writing stride?'
  stride=question('Modify stride?',bool,False)
  if stride:
    INFOS['stride']=[]
    stride=question('Enter the  *INITIAL*   output stride (e.g., "1"=write every step)',int,[1])
    INFOS['stride'].extend(stride)
    stride=question('Enter the *SUBSEQUENT* output stride (e.g., "10 2"=write every second step starting at step 10)',int,[0,1])
    INFOS['stride'].extend(stride)
    stride=question('Enter the   *FINAL*    output stride (e.g., "100 10"=write every tenth step starting at step 100)',int,[0,1])
    INFOS['stride'].extend(stride)
  else:
    INFOS['stride']=[1]


  # Add some simple keys
  INFOS['printlevel']=2
  INFOS['cwd']=os.getcwd()
  print ''

  return INFOS

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


# =========================================================================     =================
def prepare_COBRAMM(INFOS,iconddir):
  try:
    sh2cbm=open('%s/QM/COBRAMM.resources' % (iconddir), 'w')
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
  cpto='%s/QM/COBRAMM.template' % (iconddir)
  shutil.copy(cpfrom,cpto)

  cpfrom=INFOS['realtop_location']
  cpto='%s/QM/real.top' % (iconddir)
  shutil.copy(cpfrom,cpto)

  cpfrom=INFOS['modeltop_location']
  cpto='%s/QM/model-H.top' % (iconddir)
  shutil.copy(cpfrom,cpto)

  cpfrom=INFOS['realayers_location']
  cpto='%s/QM/real_layers.xyz' % (iconddir)
  shutil.copy(cpfrom,cpto)
  
  if INFOS['rattle']: 
    cpfrom=INFOS['rattlefile']
    cpto='%s/rattle' % (iconddir)
    shutil.copy(cpfrom,cpto)
  
  runname=iconddir+'/QM/runQM.sh'
  runscript=open(runname,'w')
  s='''cd QM
$SHARC/SHARC_COBRAMM.py QM.in >> QMMM.log 2>> QMMM.err
err=$?

exit $err'''
  runscript.write(s)
  runscript.close()
  os.chmod(runname, os.stat(runname).st_mode | stat.S_IXUSR)
  return

#===================================================================================================


def checktemplate_MOLPRO(filename):
  necessary=['basis','closed','occ','nelec','roots']
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


  # MOLPRO executable
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



  # MOLPRO input template
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


  # Initial wavefunction
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
  print '''Please specify the number of CPUs to be used by EACH trajectory.
'''

  # Ionization
  #print centerstring('Ionization probability by Dyson norms',60,'-')+'\n'
  #INFOS['ion']=question('Dyson norms?',bool,False)

  # wfoverlap
  if 'wfoverlap' in INFOS['needed']:
    print '\n'+centerstring('Wfoverlap code setup',60,'-')+'\n'
    INFOS['molpro.wfpath']=question('Path to wavefunction overlap executable:',str,'$SHARC/wfoverlap.x')


  # Other settings
  INFOS['molpro.gradaccudefault']=1.e-7
  INFOS['molpro.gradaccumax']=1.e-4
  INFOS['molpro.ncore']=-1
  INFOS['molpro.ndocc']=0

  return INFOS

# =================================================

def prepare_MOLPRO(INFOS,iconddir):
  # write MOLPRO.resources
  try:
    sh2pro=open('%s/QM/MOLPRO.resources' % (iconddir), 'w')
  except IOError:
    print 'IOError during prepareMOLPRO, iconddir=%s' % (iconddir)
    quit(1)
  string='''molpro %s
scratchdir %s/%s/QM
savedir %s/%s/restart
gradaccudefault %.8f
gradaccumax %f
memory %i
ncpu %i
''' % (
       INFOS['molpro'],
       INFOS['scratchdir'],
       iconddir,
       INFOS['copydir'],
       iconddir,
       INFOS['molpro.gradaccudefault'],
       INFOS['molpro.gradaccumax'],
       INFOS['molpro.mem'],
       INFOS['cobramm.ncpu']
       )
  if 'wfoverlap' in INFOS['needed']:
    string+='wfoverlap %s\n' % (INFOS['molpro.wfpath'])
  sh2pro.write(string)
  sh2pro.close()

  # copy MOs and template
  cpfrom=INFOS['molpro.template']
  cpto='%s/QM/MOLPRO.template' % (iconddir)
  shutil.copy(cpfrom,cpto)
  if INFOS['molpro.guess']:
    cpfrom=INFOS['molpro.guess']
    cpto='%s/QM/wf.init' % (iconddir)
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


  # Path to COLUMBUS directory
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




  # COLUMBUS template directory
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
  print ''

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
    mocoefmap[job]=multmap[1]
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


  # Initial mocoef
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


  # Memory
  print centerstring('COLUMBUS Memory usage',60,'-')+'\n'
  print '''Please specify the amount of memory available to COLUMBUS (in MB). For calculations including moderately-sized CASSCF calculations and less than 150 basis functions, around 2000 MB should be sufficient.
'''
  INFOS['columbus.mem']=abs(question('COLUMBUS memory:',int)[0])


  #need_wfoverlap=False
  # Ionization
  #print centerstring('Ionization probability by Dyson norms',60,'-')+'\n'
  #INFOS['ion']=question('Dyson norms?',bool,False)
  #if 'ion' in INFOS and INFOS['ion']:
    #need_wfoverlap=True
  # cioverlaps
  #if Couplings[INFOS['coupling']]['name']=='overlap':
    #need_wfoverlap=True

  # wfoverlap
  #if need_wfoverlap:
  if 'wfoverlap' in INFOS['needed']:
    if 'ion' in INFOS and INFOS['ion']:
      print 'Dyson norms requested.'
    if Couplings[INFOS['coupling']]['name']=='overlap':
      print 'Wavefunction overlaps requested.'
    print '\n'+centerstring('Wfoverlap code setup',60,'-')+'\n'
    INFOS['columbus.wfpath']=question('Path to wavefunction overlap executable:',str,'$SHARC/wfoverlap.x')
    INFOS['columbus.wfthres']=question('Determinant screening threshold:',float,[0.97])[0]
    INFOS['columbus.numfrozcore']=question('Number of frozen core orbitals for overlaps (-1=as in template):',int,[-1])[0]
    if 'ion' in INFOS and INFOS['ion']:
      INFOS['columbus.numocc']=question('Number of doubly occupied orbitals for Dyson:',int,[0])[0]

  return INFOS

# =================================================

def prepare_COLUMBUS(INFOS,iconddir):
  # write COLUMBUS.resources
  try:
    sh2col=open('%s/QM/COLUMBUS.resources' % (iconddir), 'w')
  except IOError:
    print 'IOError during prepareCOLUMBUS, directory=%i' % (iconddir)
    quit(1)
  string='''columbus %s
scratchdir %s/%s/QM
savedir %s/%s/restart
memory %i
template %s
''' % (INFOS['columbus'], INFOS['scratchdir'], iconddir, INFOS['copydir'], iconddir, INFOS['columbus.mem'], INFOS['columbus.template'])
  string+='integrals %s\n' % (INFOS['columbus.intprog'])
  for mult in INFOS['columbus.multmap']:
    string+='DIR %i %s\n' % (mult,INFOS['columbus.multmap'][mult])
  string+='\n'
  for job in INFOS['columbus.mocoefmap']:
    string+='MOCOEF %s %s\n' % (job,INFOS['columbus.mocoefmap'][job])
  string+='\n'
  if 'wfoverlap' in INFOS['needed']:
    string+='wfoverlap %s\n' % (INFOS['columbus.wfpath'])
    string+='wfthres %f\n' % (INFOS['columbus.wfthres'])
    if INFOS['columbus.numfrozcore']>=0:
      string+='numfrozcore %i\n' % (INFOS['columbus.numfrozcore'])
    if 'columbus.numocc' in INFOS:
      string+='numocc %i\n' % (INFOS['columbus.numocc'])
  else:
    string+='nooverlap\n'
  sh2col.write(string)
  sh2col.close()

  # copy MOs and template
  if INFOS['columbus.guess']:
    cpfrom=INFOS['columbus.guess']
    cpto='%s/QM/mocoef_mc.init' % (iconddir)
    shutil.copy(cpfrom,cpto)

  if INFOS['columbus.copy_template']:
    copy_from=INFOS['columbus.copy_template_from']
    copy_to=iconddir+'/QM/COLUMBUS.template/'
    if os.path.exists(copy_to):
      shutil.rmtree(copy_to)
    shutil.copytree(copy_from,copy_to)


  return

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================


def checktemplate_MOLCAS(filename,INFOS):
  necessary=['basis','ras2','nactel','inactive','cobramm']
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
  INFOS['molcas.mem']=abs(question('MOLCAS memory:',int)[0])
  print '''Please specify the number of CPUs to be used by EACH calculation.
'''




  ## Ionization
  #need_wfoverlap=False
  #print centerstring('Ionization probability by Dyson norms',60,'-')+'\n'
  #INFOS['ion']=question('Dyson norms?',bool,False)
  #if 'ion' in INFOS and INFOS['ion']:
    #need_wfoverlap=True

  # wfoverlap
  if 'wfoverlap' in INFOS['needed']:
    print '\n'+centerstring('Wfoverlap code setup',60,'-')+'\n'
    if 'ion' in INFOS and INFOS['ion']:
      print 'Dyson norms requested.'
    INFOS['molcas.wfpath']=question('Path to wavefunction overlap executable:',str,'$SHARC/wfoverlap.x')
    # TODO not asked for: numfrozcore, numocc


  return INFOS

# =================================================

def prepare_MOLCAS(INFOS,iconddir):
  # write MOLCAS.resources
  try:
    sh2cas=open('%s/QM/MOLCAS.resources' % (iconddir), 'w')
  except IOError:
    print 'IOError during prepareMOLCAS, iconddir=%s' % (iconddir)
    quit(1)
  project='MOLCAS'
  string='''molcas %s
scratchdir %s/%s/QMMM/QM
savedir %s/%s/restart
memory %i
ncpu %i
project %s''' % (INFOS['molcas'],
                 INFOS['scratchdir'],
                 iconddir,
                 INFOS['copydir'],
                 iconddir,
                 INFOS['molcas.mem'],
                 INFOS['cobramm.ncpu'],
                 project)
  if 'wfoverlap' in INFOS['needed']:
    string+='\nwfoverlap %s\n' % INFOS['molcas.wfpath']
  sh2cas.write(string)
  sh2cas.close()

  # copy MOs and template
  cpfrom=INFOS['molcas.template']
  cpto='%s/QM/MOLCAS.template' % (iconddir)
  shutil.copy(cpfrom,cpto)
  if not INFOS['molcas.guess']=={}:
    for i in INFOS['molcas.guess']:
      if INFOS['molcas.jobiph_or_rasorb']==1:
        cpfrom=INFOS['molcas.guess'][i]
        cpto='%s/QM/%s.%i.JobIph.init' % (iconddir,project,i)
      else:
        cpfrom=INFOS['molcas.guess'][i]
        cpto='%s/QM/%s.%i.RasOrb.init' % (iconddir,project,i)
      shutil.copy(cpfrom,cpto)



  return

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

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



  # QMMM
  if qmmm_job(INFOS['ADF.template'],INFOS):
    print centerstring('ADF QM/MM setup',60,'-')+'\n'
    print 'Your template specifies a QM/MM calculation. Please give the force field and connection table files.'
    while True:
      filename=question('Force field file:',str)
      if not os.path.isfile(filename):
        print 'File %s does not exist!' % (filename)
        continue
      else:
        break
    INFOS['ADF.fffile']=filename
    while True:
      filename=question('Connection table file:',str)
      if not os.path.isfile(filename):
        print 'File %s does not exist!' % (filename)
        continue
      else:
        break
    INFOS['ADF.ctfile']=filename


  # initial MOs
  print centerstring('Initial restart: MO Guess',60,'-')+'\n'
  print '''Please specify the path to an ADF TAPE21 file containing suitable starting MOs for the ADF calculation. Please note that this script cannot check whether the wavefunction file and the Input template are consistent!
'''
  if question('Do you have a restart file?',bool,True):
     if True:
       filename=question('Restart file:',str,'ADF.t21.init')
       INFOS['adf.guess']=filename
  else:
    print 'WARNING: Remember that the calculations may take longer without an initial guess for the MOs.'
    #time.sleep(2)
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
  #need_wfoverlap=False
  #print centerstring('Ionization probability by Dyson norms',60,'-')+'\n'
  #INFOS['ion']=question('Dyson norms?',bool,False)
  #if 'ion' in INFOS and INFOS['ion']:
    #need_wfoverlap=True
  #if Couplings[INFOS['coupling']]['name']=='overlap':
    #need_wfoverlap=True


  # Overlaps
  #if need_wfoverlap:
  if 'wfoverlap' in INFOS['needed']:
    print '\n'+centerstring('Wfoverlap code setup',60,'-')+'\n'
    INFOS['adf.wfoverlap']=question('Path to wavefunction overlap executable:',str,'$SHARC/wfoverlap.x')
    print ''
    print '''State threshold for choosing determinants to include in the overlaps'''
    print '''For hybrids (and without TDA) one should consider that the eigenvector X may have a norm larger than 1'''
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
  print '\n'+centerstring('Wave function analysis by TheoDORE',60,'-')+'\n'
  #INFOS['theodore']=question('TheoDORE analysis?',bool,False)
  if 'theodore' in INFOS['needed']:

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
    if 'ADF.ctfile' in INFOS:
        INFOS['theodore.count']+=7


  return INFOS

# =================================================

def prepare_ADF(INFOS,iconddir):
  # write ADF.resources
  try:
    sh2cas=open('%s/QM/ADF.resources' % (iconddir), 'w')
  except IOError:
    print 'IOError during prepareADF, iconddir=%s' % (iconddir)
    quit(1)
#  project='ADF'
  string='adfhome %s\nscmlicense %s\nscratchdir %s/%s/QM\nsavedir %s/%s/restart\nncpu %i\nschedule_scaling %f\n' % (INFOS['adf'],INFOS['scmlicense'],INFOS['scratchdir'],iconddir,INFOS['copydir'],iconddir,INFOS['adf.ncpu'],INFOS['adf.scaling'])
  if 'wfoverlap' in INFOS['needed']:
    string+='wfoverlap %s\nwfthres %f\n' % (INFOS['adf.wfoverlap'],INFOS['adf.ciothres'])
    string+='memory %i\n' % (INFOS['adf.mem'])
    #string+='numfrozcore %i\n' %(INFOS['frozcore_number'])
  else:
    string+='nooverlap\n'
  if INFOS['theodore']:
    string+='theodir %s\n' % (INFOS['adf.theodore'])
    string+='theodore_prop %s\n' % (INFOS['theodore.prop'])
    string+='theodore_fragment %s\n' % (INFOS['theodore.frag'])
  sh2cas.write(string)
  sh2cas.close()

  # copy MOs and template
  cpfrom=INFOS['ADF.template']
  cpto='%s/QM/ADF.template' % (iconddir)
  shutil.copy(cpfrom,cpto)

  if INFOS['adf.guess']:
    cpfrom1=INFOS['adf.guess']
    cpto1='%s/ADF.t21_init' % (iconddir)
    shutil.copy(cpfrom1,cpto1)

  if 'ADF.fffile' in INFOS:
    cpfrom1=INFOS['ADF.fffile']
    cpto1='%s/QM/ADF.qmmm.ff' % (iconddir)
    shutil.copy(cpfrom1,cpto1)

  if 'ADF.ctfile' in INFOS:
    cpfrom1=INFOS['ADF.ctfile']
    cpto1='%s/QM/ADF.qmmm.table' % (iconddir)
    shutil.copy(cpfrom1,cpto1)



  return

# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def checktemplate_RICC2(filename,INFOS):
  necessary=['basis','cobramm']
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

  print centerstring('Path to ORCA',60,'-')+'\n'
  path=os.getenv('ORCADIR')
  if path=='':
    path=None
  else:
    path='$ORCADIR/'
  print '\nPlease specify path to ORCA directory (SHELL variables and ~ can be used, will be expanded when interface is started).\nORCA is necessary for the calculation of spin-orbit couplings with ricc2.\n'
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
  print '''Please specify the amount of memory available to Turbomole.
'''
  INFOS['ricc2.mem']=abs(question('RICC2 memory (MB):',int,[1000])[0])
  print '''Please specify the number of CPUs to be used by EACH trajectory.
'''

  if INFOS['laser']:
    guess=2
  else:
    guess=1
  a=['','(recommended)']
  print 'For response-based methods like CC2 and ADC(2), dipole moments and transition dipole moments carry a significant computational cost. In order to speed up calculations, the interface can restrict the calculation of these properties.'
  print '''Choose one of the following dipolelevels:
0       only calculate dipole moments which are for free                                %s
1       additionally, calculate transition dipole moments involving the ground state    %s
2       calculate all elements possible with the method                                 %s
''' % (a[guess==0],a[guess==1],a[guess==2])
  INFOS['ricc2.dipolelevel']=question('Dipole level:',int,[guess])[0]


  if 'wfoverlap' in INFOS['needed']:
    print 'Wavefunction overlaps requested.'
    INFOS['ricc2.wfpath']=question('Path to wfoverlap executable:',str,'$SHARC/wfoverlap.x')
    print ''
    print '''Give threshold for choosing determinants to include in the overlaps'''
    INFOS['ricc2.wfthres']=question('Threshold:',float,[0.99])[0]



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
    sh2cc2=open('%s/QM/RICC2.resources' % (iconddir), 'w')
  except IOError:
    print 'IOError during prepare_RICC2, iconddir=%s' % (iconddir)
    quit(1)
  string='''turbodir %s
orcadir %s
scratchdir %s/%s/QM
memory %i
ncpu %i
dipolelevel %i
''' % (INFOS['turbomole'],
       INFOS['orca'],
       INFOS['scratchdir'],
       iconddir,
       INFOS['ricc2.mem'],
       INFOS['cobramm.ncpu'],
       INFOS['ricc2.dipolelevel'])
  if 'wfoverlap' in INFOS['needed']:
    string+='wfoverlap %s\n' % (INFOS['ricc2.wfpath'])
    string+='wfthres %f\n' %(INFOS['ricc2.wfthres'])
  else:
    string+='nooverlap\n'
  if 'theodore' in INFOS['needed']:
    string+='theodir %s\n' % (INFOS['ricc2.theodore'])
    string+='theodore_prop %s\n' % (INFOS['theodore.prop'])
    string+='theodore_fragment %s\n' % (INFOS['theodore.frag'])

  sh2cc2.write(string)
  sh2cc2.close()

  # copy MOs and template
  cpfrom=INFOS['ricc2.template']
  cpto='%s/QM/RICC2.template' % (iconddir)
  shutil.copy(cpfrom,cpto)
  if INFOS['ricc2.guess']:
    cpfrom1=INFOS['ricc2.guess']
    cpto1='%s/QM/mos.init' % (iconddir)
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
  - TAPE21
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
    print '\n'+centerstring('Wfoverlap code setup',60,'-')+'\n'
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
    sh2cas=open('%s/QM/GAUSSIAN.resources' % (iconddir), 'w')
  except IOError:
    print 'IOError during prepareGAUSSIAN, iconddir=%s' % (iconddir)
    quit(1)
#  project='GAUSSIAN'
  string='groot %s\nscratchdir %s/%s/QM\nsavedir %s/%s/restart\nncpu %i\nschedule_scaling %f\n' % (INFOS['groot'],INFOS['scratchdir'],iconddir,INFOS['scratchdir'],iconddir,INFOS['gaussian.ncpu'],INFOS['gaussian.scaling'])
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
  cpto='%s/QM/GAUSSIAN.template' % (iconddir)
  shutil.copy(cpfrom,cpto)

  if INFOS['gaussian.guess']:
    cpfrom1=INFOS['gaussian.guess']
    cpto1='%s/QM/GAUSSIAN.chk.init' % (iconddir)
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
    sh2cas=open('%s/QM/ORCA.resources' % (iconddir), 'w')
  except IOError:
    print 'IOError during prepareORCA, iconddir=%s' % (iconddir)
    quit(1)
#  project='ORCA'
  string='orcadir %s\nscratchdir %s/%s/QM\nsavedir %s/%s/restart\nncpu %i\nschedule_scaling %f\n' % (INFOS['orcadir'],INFOS['scratchdir'],iconddir,INFOS['copydir'],iconddir,INFOS['orca.ncpu'],INFOS['orca.scaling'])
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
  sh2cas.write(string)
  sh2cas.close()

  # copy MOs and template
  cpfrom=INFOS['ORCA.template']
  cpto='%s/QM/ORCA.template' % (iconddir)
  shutil.copy(cpfrom,cpto)

  if INFOS['orca.guess']:
    cpfrom1=INFOS['orca.guess']
    cpto1='%s/QM/ORCA.gbw.init' % (iconddir)
    shutil.copy(cpfrom1,cpto1)

  if 'ORCA.fffile' in INFOS:
    cpfrom1=INFOS['ORCA.fffile']
    cpto1='%s/QM/ORCA.qmmm.ff' % (iconddir)
    shutil.copy(cpfrom1,cpto1)

  if 'ORCA.ctfile' in INFOS:
    cpfrom1=INFOS['ORCA.ctfile']
    cpto1='%s/QM/ORCA.qmmm.table' % (iconddir)
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
  INFOS['dipolelevel']=question('Requested dipole level:',int,[0])[0]
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
    sh2cas=open('%s/QM/BAGEL.resources' % (iconddir), 'w')
  except IOError:
    print 'IOError during prepareBAGEL, iconddir=%s' % (iconddir)
    quit(1)
  project='BAGEL'
  string='bagel %s\npyquante %s\nscratchdir %s/%s/QM\nmemory %i\nncpu %i\ndipolelevel %i\nproject %s' % (INFOS['bagel'],INFOS['bagel.pyq'],INFOS['scratchdir'],iconddir,INFOS['bagel.mem'],INFOS['bagel.ncpu'],INFOS['dipolelevel'],project)
  
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
  cpto='%s/QM/BAGEL.template' % (iconddir)
  shutil.copy(cpfrom,cpto)
  if not INFOS['bagel.guess']=={}:
    for i in INFOS['bagel.guess']:
      cpfrom=INFOS['bagel.guess'][i]
      cpto='%s/QM/%s.%i.init' % (iconddir,'archive',i)
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
    INFOS['copydir']=INFOS['cwd']
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

def writeSHARCinput(INFOS,initobject,iconddir,istate):

  inputfname=iconddir+'/input'
  try:
    inputf=open(inputfname, 'w')
  except IOError:
    print 'IOError during writeSHARCinput, iconddir=%s\n%s' % (iconddir,inputfname)
    quit(1)

  s='printlevel 2\n\ngeomfile "geom"\nveloc external\nvelocfile "veloc"\n\n'
  s+='nstates '
  for nst in INFOS['states']:
    s+='%i ' % nst
  s+='\nactstates '
  for nst in INFOS['actstates']:
    s+='%i ' % nst
  s+='\nstate %i %s\n' % (istate,['mch','diag'][INFOS['diag']])
  s+='coeff auto\n'
  s+='rngseed %i\n\n' % (random.randint(-32768,32767))
  s+='ezero %18.10f\n' % (INFOS['eref'])

  s+='tmax %f\nstepsize %f\nnsubsteps %i\n' % (INFOS['tmax'],INFOS['dtstep'],INFOS['nsubstep'])
  if INFOS['kill']:
    s+='killafter %f\n' % (INFOS['killafter'])
  s+='\n'

  if INFOS['atommaskarray']:
    s+='atommask external\natommaskfile "atommask"\n\n'
  
  if INFOS['rattle']:
    s+='rattle\nrattlefile "%s"\n' % (INFOS['rattlefile'])

  s+='surf %s\n' % (INFOS['surf'])
  s+='coupling %s\n' % (Couplings[INFOS['coupling']]['name'])
  s+='%sgradcorrect\n' % (['no',''][INFOS['gradcorrect']])
  s+='ekincorrect %s\n' % (EkinCorrect[INFOS['ekincorrect']]['name'])
  s+='reflect_frustrated %s\n' % (EkinCorrect[INFOS['reflect']]['name'])
  s+='decoherence_scheme %s\n' % (INFOS['decoherence'][0])
  if INFOS['decoherence'][1]:
    s+='decoherence_param %s\n' % (INFOS['decoherence'][1])
  s+='hopping_procedure %s\n' % (INFOS['hopping'])
  if INFOS['force_hops']:
    s+='force_hop_to_gs %f\n' % (INFOS['force_hops_dE'])
  if INFOS['scaling']:
    s+='scaling %f\n' % (INFOS['scaling'])
  if INFOS['damping']:
    s+='dampeddyn %f\n' % (INFOS['damping'])
  if INFOS['phases_from_interface']:
    s+='phases_from_interface\n'
  if INFOS['pysharc']:
    s+='notrack_phase\n'

  if INFOS['sel_g']:
    s+='grad_select\n'
  else:
    s+='grad_all\n'
  if INFOS['sel_t']:
    s+='nac_select\n'
  else:
    if Couplings[INFOS['coupling']]['name']=='ddr' or INFOS['gradcorrect'] or EkinCorrect[INFOS['ekincorrect']]['name']=='parallel_nac':
      s+='nac_all\n'
  if 'eselect' in INFOS:
    s+='eselect %f\n' % (INFOS['eselect'])
#  if Interfaces[INFOS['interface']]['script']=='SHARC_COLUMBUS.py':
#    s+='select_directly\n'
#  if Interfaces[INFOS['interface']]['script']=='SHARC_ADF.py':
#    s+='select_directly\n'
#  if Interfaces[INFOS['interface']]['script']=='SHARC_GAUSSIAN.py':
#    s+='select_directly\n'
#  if Interfaces[INFOS['interface']]['script']=='SHARC_RICC2.py':
#    s+='select_directly\n'
#  if Interfaces[INFOS['interface']]['script']=='SHARC_MOLPRO.py':
#    s+='select_directly\n'
#  if Interfaces[INFOS['interface']]['script']=='SHARC_MOLCAS.py':
#    s+='select_directly\n'
#  if Interfaces[INFOS['interface']]['script']=='SHARC_ORCA.py':
#    s+='select_directly\n'
  # TODO: maybe could be deactivated in the future for some interfaces
  s+='select_directly\n'

  if not INFOS['soc']:
    s+='nospinorbit\n'

  if INFOS['write_grad']:
    s+='write_grad\n'
  if INFOS['write_NAC']:
    s+='write_nacdr\n'
  if INFOS['write_overlap']:
    s+='write_overlap\n'
  if INFOS['write_property1d']:
    s+='write_property1d\n'
    if 'theodore.count' in INFOS:
      s+='n_property1d %i\n' % (INFOS['theodore.count'])
    else:
      s+='n_property1d %i\n' % (1)
  if INFOS['write_property2d']:
    s+='write_property2d\n'
    s+='n_property2d %i\n' % (1)

  # NetCDF or ASCII
  if INFOS['netcdf']:
    out='netcdf'
  else:
    out='ascii'
  s+='output_format %s\n' % out

  # stride
  if 'stride' in INFOS:
    s+='output_dat_steps'
    for i in  INFOS['stride']:
      s+=' %i' % i
    s+='\n'

  # laser
  if INFOS['laser']:
    s+='laser external\n'
    s+='laserfile "laser"\n'
    if INFOS['dipolegrad']:
      s+='dipole_gradient'

  if 'ion' in INFOS and INFOS['ion']:
    s+='ionization\n'
    s+='ionization_step 1\n'

  if 'theodore' in INFOS and INFOS['theodore']:
    s+='theodore\n'
    s+='theodore_step 1\n'

  inputf.write(s)
  inputf.close()

  # geometry file
  geomfname=iconddir+'/geom'
  geomf=open(geomfname,'w')
  for atom in initobject.atomlist:
    geomf.write(atom.geomstring()+'\n')
  geomf.close()

  # velocity file
  velocfname=iconddir+'/veloc'
  velocf=open(velocfname,'w')
  for atom in initobject.atomlist:
    velocf.write(atom.velocstring()+'\n')
  velocf.close()

  # laser file
  if INFOS['laser']:
    laserfname=iconddir+'/laser'
    shutil.copy(INFOS['laserfile'],laserfname)

  # atommask file
  if INFOS['atommaskarray']:
    atommfname=iconddir+'/atommask'
    atommf=open(atommfname,'w')
    for i,atom in enumerate(initobject.atomlist):
      if i+1 in INFOS['atommaskarray']:
        atommf.write('T\n')
      else:
        atommf.write('F\n')
    atommf.close()

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
    projname='traj_%5s' % (iconddir[-6:-1])

  # ================================
  intstring=''
  if 'adfrc' in INFOS:
    intstring='. %s\nexport PYTHONPATH=$ADFHOME/scripting:$PYTHONPATH' % (INFOS['adfrc'])

  # ================================
  if INFOS['pysharc']:
    driver=Interfaces[ INFOS['interface'] ]['pysharc_driver']
    exestring='. $SHARC/sharcvars.sh\n$ANACONDA/bin/python2 $SHARC/%s input' % driver
  else:
    exestring='$SHARC/sharc.x input'


  # ================================ for here mode
  if INFOS['here']:
    string='''#!/bin/bash

#$-N %s

%s

PRIMARY_DIR=%s/%s

cd $PRIMARY_DIR

%s
''' % (projname,intstring,INFOS['cwd'],iconddir,exestring)
  #
  # ================================ for remote mode
  else:
    string='''#!/bin/bash

#$-N %s
''' % (projname)
    if INFOS['qsub']:
      string+='#$ -v USER_EPILOG=%s/epilog.sh' % (iconddir)

    string+='''
%s

PRIMARY_DIR=%s/%s
COPY_DIR=%s/%s

mkdir -p $COPY_DIR
cp -r $PRIMARY_DIR/* $COPY_DIR
cd $COPY_DIR
echo $HOSTNAME > $PRIMARY_DIR/host_info
echo $(pwd) >> $PRIMARY_DIR/host_info
echo $(date) >> $PRIMARY_DIR/host_info

%s
err=$?

cp -r $COPY_DIR/output.* $COPY_DIR/restart.* $COPY_DIR/restart/ $PRIMARY_DIR

if [ $err == 0 ];
then
  rm -r $COPY_DIR
else
  echo "The calculation crashed at
date = $(date)
with error code $err.
Please inspect the trajectory on
host = $HOSTNAME
in
dir  = $(pwd)
" > $PRIMARY_DIR/README
fi
''' % (intstring,INFOS['cwd'], iconddir, INFOS['copydir'], iconddir,exestring)

  runscript.write(string)
  runscript.close()
  filename=iconddir+'/run.sh'
  os.chmod(filename, os.stat(filename).st_mode | stat.S_IXUSR)

  # also write an epilog script
  if not INFOS['here'] and INFOS['qsub']:
    try:
      episcript=open(iconddir+'/epilog.sh','w')
      string='''#/bin/bash

PRIMARY_DIR=%s/%s
COPY_DIR=%s/%s

cp $COPY_DIR/output.* $COPY_DIR/restart.* $PRIMARY_DIR
rm -r $COPY_DIR
''' % (INFOS['cwd'], iconddir, INFOS['copydir'], iconddir)
      episcript.write(string)
      episcript.close()
    except IOError:
      print 'Could not write epilog script for %s.' % (iconddir)
  return


# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def get_iconddir(istate,INFOS):
  if INFOS['diag']:
    dirname='State_%i' % (istate)
  else:
    mult,state,ms=INFOS['statemap'][istate]
    dirname=IToMult[mult]+'_%i' % (state-(mult==1 or mult==2))
  return dirname

# ====================================

def setup_all(INFOS):
  '''This routine sets up the directories for the initial calculations.'''

  string='\n  '+'='*80+'\n'
  string+='||'+centerstring('Setting up directories...',80)+'||\n'
  string+='  '+'='*80+'\n\n'
  print string

  all_run=open('all_run_traj.sh','w')
  string='#/bin/bash\n\nCWD=%s\n\n' % (INFOS['cwd'])
  all_run.write(string)
  if INFOS['qsub']:
    all_qsub=open('all_qsub_traj.sh','w')
    string='#/bin/bash\n\nCWD=%s\n\n' % (INFOS['cwd'])
    all_qsub.write(string)

  for istate in INFOS['setupstates']:
    dirname=get_iconddir(istate,INFOS)
    io=make_directory(dirname)
    if io!=0:
      print 'Could not make directory %s' % (dirname)
      quit(1)

  width=50
  ntraj=INFOS['ntraj']
  idone=0
  finished=False

  initlist=INFOS['initlist']

  for icond in range(INFOS['firstindex'],INFOS['ninit']+1):

    for istate in INFOS['setupstates']:

      if len(initlist[icond-1].statelist)<istate:
        continue
      if not initlist[icond-1].statelist[istate-1].Excited:
        continue

      idone+=1

      done=idone*width/ntraj
      sys.stdout.write('\rProgress: ['+'='*done+' '*(width-done)+'] %3i%%' % (done*100/width))

      dirname=get_iconddir(istate,INFOS)+'/TRAJ_%05i/' % (icond)
      io=make_directory(dirname)
      if io!=0:
        print 'Skipping initial condition %i %i!' % (istate, icond)
        continue

      writeSHARCinput(INFOS,initlist[icond-1],dirname,istate)
      io=make_directory(dirname+'/QM')
      io+=make_directory(dirname+'/restart')
      if io!=0:
        print 'Could not make QMMM or restart directory!'
        continue
      globals()[Interfaces[ INFOS['interface']]['prepare_routine'] ](INFOS,dirname)

      prepare_COBRAMM(INFOS,dirname)
      writeRunscript(INFOS,dirname)

      string='cd $CWD/%s/\nbash run.sh\ncd $CWD\necho %s >> DONE\n' % (dirname,dirname)
      all_run.write(string)
      if INFOS['qsub']:
        string='cd $CWD/%s/\n%s run.sh\ncd $CWD\n' % (dirname,INFOS['qsubcommand'])
        all_qsub.write(string)

      if idone==ntraj:
        finished=True
        break
    if finished:
      print '\n\n%i trajectories setup, last initial condition was %i in state %i.\n' % (ntraj,icond,istate)
      setup_stat=open('setup_traj.status','a+')
      string='''*** %s %s %s
  First index:          %i
  Last index:           %i
  Trajectories:         %i
  State of last traj.:  %i

''' % (datetime.datetime.now(),
       gethostname(),
       os.getcwd(),
       INFOS['firstindex'],
       icond,
       ntraj,
       istate)
      setup_stat.write(string)
      setup_stat.close()
      break

  all_run.close()
  filename='all_run_traj.sh'
  os.chmod(filename, os.stat(filename).st_mode | stat.S_IXUSR)
  if INFOS['qsub']:
    all_qsub.close()
    filename='all_qsub_traj.sh'
    os.chmod(filename, os.stat(filename).st_mode | stat.S_IXUSR)

  print '\n'


# ======================================================================================================================
# ======================================================================================================================
# ======================================================================================================================

def main():
  '''Main routine'''

  usage='''
python setup_traj.py

This interactive program prepares SHARC dynamics calculations.
'''

  description=''
  parser = OptionParser(usage=usage, description=description)

  displaywelcome()
  open_keystrokes()

  INFOS=get_general()
  INFOS=get_runscript_info(INFOS)

  print '\n'+centerstring('Full input',60,'#')+'\n'
  for item in INFOS:
    if not 'initlist' in item:
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

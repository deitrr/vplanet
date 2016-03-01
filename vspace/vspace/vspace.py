from __future__ import print_function
import os 
import sys
import subprocess as sb
import numpy as np
import re
from IPython.core.debugger import Tracer
import itertools as it


if len(sys.argv) < 2:
  raise IOError("Must enter an input file name")
else:
  inputf = sys.argv[1]

try:
  f = open(inputf,'r')
except IOError:
  print ("%s is not a valid file name. Please reenter." %inputf)
  
lines = f.readlines()
f.close()

src = None
dest = None
trial = 'default'
fnum = 0   #number of files to be copied/edited
flist = [] #list of files
fline = [] #list of line numbers associated with each file
iter_var = []  #list that will contain all iteration variables.
iter_file = []  #which file each variable belongs to
iter_name = []  #the name of each variable
numtry = 1     #number of trials to be generated
numvars = 0   #number of iter_vars
for i in range(len(lines)):
  if lines[i].split() == []:
    pass  #nothing on this line
  elif lines[i].split()[0] == 'srcfolder':
    src = lines[i].split()[1]
  elif lines[i].split()[0] == 'destfolder':
    dest = lines[i].split()[1]
  elif lines[i].split()[0] == 'trialname':
    trial = lines[i].split()[1]
  elif lines[i].split()[0] == 'file':
    flist.append(lines[i].split()[1])
    fline.append(i)
    fnum += 1
    
  if re.search('\[',lines[i]) != None:
    spl = re.split('[\[\]]',lines[i])
    name = lines[i].split()[0]
    values = spl[1].split(',')
    if len(values) != 3:
      raise IOError("Attempt to iterate over '%s' for '%s', but incorrect number of values provided. Syntax should be [<low>, <high>, <spacing>], [<low>, <high>, n<number of points>], or [<low>, <high>, l<number of points>] (log spacing) "%(name,flist[fnum-1]))
    
    for j in range(len(values)):
      values[j] = values[j].strip()  #remove any leading white spaces

    if values[2][0] == 'n':
      values[2] = values[2][1:]
      if np.float(values[2]).is_integer():
        array = np.linspace(np.float(values[0]),np.float(values[1]),np.float(values[2]))
      else:
        raise IOError("Attempt to iterate over '%s' for '%s', but number of points provided not an integer value"%(name,flist[fnum-1]))
    elif values[2][0] == 'l':
      values[2] = values[2][1:]
      if np.float(values[0]) < 0:
        if np.float(values[2]).is_integer():
          array = -np.logspace(np.log10(-np.float(values[0])),np.log10(-np.float(values[1])),np.float(values[2]))
        else:
          raise IOError("Attempt to iterate over '%s' for '%s', but number of points provided not an integer value"%(name,flist[fnum-1]))
      else:
        if np.float(values[2]).is_integer():
          array = np.logspace(np.log10(np.float(values[0])),np.log10(np.float(values[1])),np.float(values[2]))
        else:
          raise IOError("Attempt to iterate over '%s' for '%s', but number of points provided not an integer value"%(name,flist[fnum-1]))
      
        
    else:
      if np.float(values[0]) > np.float(values[1]) and np.float(values[2]) > 0:
        raise IOError("Attempt to iterate over '%s' for '%s', but start value > end value and spacing is positive"%(name,flist[fnum-1]))
      else:
        array = np.arange(np.float(values[0]),np.float(values[1])+np.float(values[2]),np.float(values[2]))
    iter_var.append(array)
    iter_file.append(fnum-1)
    iter_name.append(name)
    numtry *= len(array)
    numvars += 1
      
fline.append(i+1)

#error handling
if src == None:
  raise IOError("Name of source folder not provided in file '%s'. Use syntax 'srcfolder <foldername>'"%inputf)
if dest == None:
  raise IOError("Name of destination folder not provided in file '%s'. Use syntax 'destfolder <foldername>'"%inputf)
if flist == []:
  raise IOError("No files-to-be-copied provided in file '%s'. Use syntax 'file <filename>'")
  
if not os.path.exists(src):
  raise IOError("Source folder '%s' does not exist"%src)
  
if re.search('\/',dest) != None:
  dest_parent = '/'.join(dest.split("/")[0:-1])
  if not os.path.exists(dest_parent):
    raise IOError("Destination's parent folder '%s' does not exist"%dest_parent)
if not os.path.exists(dest):
  os.system('mkdir '+dest)

if numvars == 0:  
  for i in range(fnum):
    #read in file to be copied
    try:
      fIn = open(src+'/'+flist[i],'r')
    except IOError:
      print "%s is not a valid file name. Please reenter."%(src+'/'+flist[i])
  
    #find the lines in 'inputf' that correspond to this file
    slines = lines[fline[i]+1:fline[i+1]]
    slines = [slines[j] for j in range(len(slines)) if slines[j].split() != []]
    sflag = np.zeros(len(slines)) - 1  #bookkeeping (has option been found in file?)
    spref = [slines[j].split()[0] for j in range(len(slines))] #option name
    dlines = fIn.readlines()
    fIn.close()

    #create file out
    fOut = open(dest+'/'+flist[i],'w')
  
    for j in range(len(dlines)):
      for k in range(len(spref)):
        if dlines[j].split()[0] == spref[k]:
          sflag[k] = 1   #option in file to be copied matched with option from inputf
          dlines[j] = slines[k]
          if dlines[j][-1] != '\n' and j < (len(dlines)-1):
           dlines[j] = dlines[j]+'\n'  #add a newline, just in case
        elif dlines[j].split()[0] == 'rm':
          #remove an option by placing a comment!
          if dlines[j].split()[1] == spref[k]:
            dlines[j] = '#'+dlines[j]
            sflag[k] = 1
  
      fOut.write(dlines[j])  #write to the copied file
    
    for k in range(len(spref)):
      #check if any options were not already present in the copied file, then write them
      if sflag[k] < 0:
        fOut.write('\n'+slines[k])
  
    fOut.close()

elif numvars >= 1:
  count = 0  #suffix of directory to be generated
  iter_file = np.array(iter_file)  #convert to numpy array for useful methods
  iter_name = np.array(iter_name)
  iterables0 = [x for x in iter_var]  #create list of iterables

  #iterate over all possible combinations
  for tup in it.product(*iterables0):
    destfull = dest+'/'+trial+'%d'%count  #create directory for this combination
    if not os.path.exists(destfull):
       os.system('mkdir '+destfull)
    
    for i in range(fnum):
      #read in file to be copied
      try:
        fIn = open(src+'/'+flist[i],'r')
      except IOError:
        print "%s is not a valid file name. Please reenter."%(src+'/'+flist[i])
      
      #find the lines in 'inputf' that correspond to this file
      slines = lines[fline[i]+1:fline[i+1]]
      slines = [slines[j] for j in range(len(slines)) if slines[j].split() != []]
      sflag = np.zeros(len(slines)) - 1   #bookkeeping (has option been found in file?)
      spref = [slines[j].split()[0] for j in range(len(slines))]   #option name
      dlines = fIn.readlines()
      fIn.close()
      
      #create file out
      fOut = open(destfull+'/'+flist[i],'w')
          
      for j in range(len(dlines)):
        for k in range(len(spref)):
          if dlines[j].split()[0] == spref[k]:
            sflag[k] = 1   #option in file to be copied matched with option from inputf
            dlines[j] = slines[k]
            for m in range(len(iter_file)):
              if iter_file[m] == i and iter_name[m] == dlines[j].split()[0]:
                dlines[j] = dlines[j].split()[0] + ' ' +str(tup[m])
            if dlines[j][-1] != '\n' and j < (len(dlines)-1):
              dlines[j] = dlines[j]+'\n'  #add a newline, just in case
          elif slines[k].split()[0] == 'rm':
            #remove an option by placing a comment!
            if dlines[j].split()[0] == slines[k].split()[1]:
              dlines[j] = '#'+dlines[j]
              sflag[k] = 1
              
        fOut.write(dlines[j])  #write to the copied file
      
      for k in range(len(spref)):
        #check if any were not already present in the copied file, then write them
        if sflag[k] < 0:
          fOut.write('\n'+slines[k])
      
    fOut.close()        
    count += 1    #move to next combination

  
"""
Batch compute an arbitrary filter  

@author: Nick Gibbons
"""

from paraview.simple import *
from glob import glob
from sys import argv

scriptname = argv[1]
pattern = argv[2]
files = glob(pattern)
files.sort()

with open(scriptname) as fp:
     script = fp.read()

for i,filename in enumerate(files):
    print "Computing file: ", filename, scriptname

    soln = OpenDataFile(filename)
    prgfil = ProgrammableFilter(soln)
    my_script = script.replace('REPLACE',str(i).zfill(3))
    prgfil.Script = my_script

    Show(prgfil)
    Delete(prgfil)
    Delete(soln)

print "... Done"


import csv
import numpy
from matplotlib import pyplot
import sys,os
files = os.listdir(sys.argv[1])
data = []
for file in files:
   if file.endswith('.dat'):
      reader = csv.reader(open(sys.argv[1]+'/'+file, "rb"), delimiter=' ',skipinitialspace=True, quoting=csv.QUOTE_NONE)
      run = []
      for row in reader:
          run.append([float(x) for x in row[1:6]])
      data.append(run[-300:])
#print data
data = numpy.array(data)
avg = data.mean(axis=0)
std = data.std(axis=0)

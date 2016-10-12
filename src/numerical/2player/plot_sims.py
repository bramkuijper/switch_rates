#!/usr/bin/env python3

import sys, re, os.path, itertools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np

from pylab import *

params={'axes.linewidth' : .5}
rcParams.update(params)

# get the filename from the command line
filename = sys.argv[1]

# first figure out last line of dataset
f = open(filename)
fl = f.readlines()
f.close()

nrow = len(fl)-8

def find_parline(filename, searchvar):

    f = open(filename)
    fl = f.readlines()
    f.close()

    rx = "^" + searchvar + ".*"

    for row, line in enumerate(fl):

        if re.match("^" + searchvar + ".*", line) is not None:
            return(row)

    return(-1)

parline = find_parline(filename, "mum")

if parline == -1:
    print("failed to get params")
    sys.exit(1)

# read the parameters
parameters = pd.read_csv(filename, 
                            sep=";",
                            skiprows=parline-1,
                            header=None,
                            names=["var","val"])

s1 = float(parameters[parameters["var"] == "s1"]["val"])
s2 = float(parameters[parameters["var"] == "s2"]["val"])
q = float(parameters[parameters["var"] == "q"]["val"])

print(s1)
print(s2)


data = pd.read_csv(filename, sep=";", nrows=parline-3)


# initialize and specify size 
fig = plt.figure(figsize=(10,20))

num_rows = 9

# names of columns in dataframe
colnames = list(data.columns.values)

# add first subplot depicting switch rate modulation
plt.subplot(num_rows,1,1)

plt.plot(data["time"], np.exp(q * data["zs1"])*s1,'b',
        data["time"], np.exp(q * data["zs2"])*s2,'r',
        linewidth=1)

plt.legend((r'$\hat{s}_{1}$',r'$\hat{s}_{2}$'))
plt.ylabel(r'Switch rate modulation')

# add subplot depicting zs values
plt.subplot(num_rows,1,2)

plt.plot(data["time"],data["zs1"],'b',
        data["time"],data["zs2"],'r',
        linewidth=1)

plt.legend((r'$z_{s_{1}}$',r'$z_{s_{2}}$'))
plt.ylabel(r'Switch rate modulation')

# add 2nd subplot depicting p1, p2
plt.subplot(num_rows,1,3)

plt.plot(data["time"],data["p1"],'b',
        data["time"],data["p2"],'r',
        linewidth=1)
plt.legend((r'$p_{1}$',r'$p_{2}$'))
plt.ylabel(r'% $z_{1}$ offspring')
plt.ylim(-0.05,1.05)

# add 3rd subplot depicting patch frequencies in e1
plt.subplot(num_rows,1,4)

plt.plot(data["time"],data["f1aa"],'c',
        data["time"],data["f1am"],'m',
        data["time"],data["f1mm"],'y',
        linewidth=1)
plt.legend((r'$f_{1aa}$',r'$f_{1am}$',r'$f_{1mm}$'))
plt.ylabel(r'Patch frequencies $e_{1}$')
plt.ylim(-0.05,1.05)

# add subplot depicting patch frequencies in e2
plt.subplot(num_rows,1,5)

plt.plot(data["time"],data["f2aa"],'c',
        data["time"],data["f2am"],'m',
        data["time"],data["f2mm"],'y',
        linewidth=1)
plt.legend((r'$f_{2aa}$',r'$f_{2am}$',r'$f_{2mm}$'))
plt.ylabel(r'Patch frequencies $e_{2}$')
plt.ylim(-0.05,1.05)

# add  subplot depicting reproductive values in e1
plt.subplot(num_rows,1,6)

plt.plot(data["time"],data["v1aa"],'c',
        data["time"],data["v1am"],'m',
        data["time"],data["v1ma"],'y',
        data["time"],data["v1mm"],'k',
        linewidth=1)
plt.legend((r'$v_{1aa}$',r'$v_{1am}$',r'$v_{1ma}$',r'$v_{1mm}$'))
plt.ylabel(r'Reproductive value $e_{1}$')

# add  subplot depicting reproductive values in e2
plt.subplot(num_rows,1,7)

plt.plot(data["time"],data["v2aa"],'c',
        data["time"],data["v2am"],'m',
        data["time"],data["v2ma"],'y',
        data["time"],data["v2mm"],'k',
        linewidth=1)
plt.legend((r'$v_{2aa}$',r'$v_{2am}$',r'$v_{2ma}$',r'$v_{2mm}$'))
plt.ylabel(r'Reproductive value $e_{2}$')


# add relatedness in in e1
plt.subplot(num_rows,1,8)

plt.plot(data["time"],data["raa1"],'c',
        data["time"],data["ram1"],'m',
        data["time"],data["rmm1"],'y',
        linewidth=1)
plt.legend((r'$r_{aa,1}$',r'$r_{am,1}$',r'$r_{mm,1}$'))
plt.ylabel(r'Relatedness $e_{1}$')

# add relatedness in in e2
plt.subplot(num_rows,1,9)

plt.plot(data["time"],data["raa2"],'c',
        data["time"],data["ram2"],'m',
        data["time"],data["rmm2"],'y',
        linewidth=1)
plt.legend((r'$r_{aa,2}$',r'$r_{am,2}$',r'$r_{mm,2}$'))
plt.ylabel(r'Relatedness $e_{2}$')

graphname = os.path.dirname(filename)
if graphname != '':
    graphname += "/"
graphname += "graph_" + os.path.basename(filename) + ".pdf"

plt.savefig(graphname,format="pdf")

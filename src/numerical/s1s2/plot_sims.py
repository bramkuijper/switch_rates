#!/usr/bin/env python3

import sys, re, os.path, itertools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pandas as pd

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

print(s1)
print(s2)


data = pd.read_csv(filename, sep=";", nrows=parline-3)


# initialize and specify size 
fig = plt.figure(figsize=(10,10))

num_rows = 4

# names of columns in dataframe
colnames = list(data.columns.values)

# add first subplot depicting switch rate modulation
plt.subplot(num_rows,1,1)

plt.plot(data["time"],data["zs1"]*s1,'b',
        data["time"],data["zs2"]*s2,'r',
        linewidth=1)

plt.legend((r'$x_{s_{1}}$',r'$x_{s_{2}}$'))
plt.ylabel(r'Switch rate modulation')

# add 2nd subplot depicting p1, p2
plt.subplot(num_rows,1,2)

plt.plot(data["time"],data["p1"],'b',
        data["time"],data["p2"],'r',
        linewidth=1)
plt.legend((r'$p_{1}$',r'$p_{2}$'))
plt.ylabel(r'% $z_{1}$ offspring')
plt.ylim(-0.05,1.05)

# add 3rd subplot depicting patch frequencies
plt.subplot(num_rows,1,3)

plt.plot(data["time"],data["f1a"],'c',
        data["time"],data["f1m"],'m',
        data["time"],data["f2a"],'y',
        data["time"],data["f2m"],'k',
        linewidth=1)
plt.legend((r'$f_{1,a}$',r'$f_{1,m}$',r'$f_{2,a}$',r'$f_{2,m}$'))
plt.ylabel(r'Patch frequencies')
plt.ylim(-0.05,1.05)

# add 4th subplot depicting reproductive values
plt.subplot(num_rows,1,4)

plt.plot(data["time"],data["v1a"],'c',
        data["time"],data["v1m"],'m',
        data["time"],data["v2a"],'y',
        data["time"],data["v2m"],'k',
        linewidth=1)
plt.legend((r'$v_{1,a}$',r'$v_{1,m}$',r'$v_{2,a}$',r'$v_{2,m}$'))
plt.ylabel(r'Reproductive value')

graphname = os.path.dirname(filename)
if graphname != '':
    graphname += "/"
graphname += "graph_" + os.path.basename(filename) + ".pdf"

plt.savefig(graphname,format="pdf")

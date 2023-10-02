#!/bin/python3
from  globaldata import GlobalData
from matplotlib import pyplot as pp

def graph2sp3():
    x = [x for x in GlobalData.mRadii[10]]
    y = [GlobalData.mRadii[10][z] for z in x ]
    fig, ax = pp.subplots()
    ax.plot(x,y)
    ax.scatter(x,y)
    ax.set_xlabel('Z')
    ax.set_ylabel('Radii (Angstrom)')
    ax.set_title('2SP3 Radius vs. Atomic Number')
    #Label
    for i, j in zip(x, y):
        ax.annotate(str(i), (i, j), textcoords="offset points", xytext=(0, 10), ha='center')
    #Save
    fig.set_figwidth(15)
    fig.set_figheight(5) 
    pp.savefig('testfig', dpi=400)

def graph2sp3():
    x = [x for x in GlobalData.mRadii[10]]
    y = [GlobalData.mRadii[10][z] for z in x ]
    fig, ax = pp.subplots()
    ax.plot(x,y)
    ax.scatter(x,y)
    ax.set_xlabel('Z')
    ax.set_ylabel('Radii (Angstrom)')
    ax.set_title('2SP3 Radius vs. Atomic Number')
    #Label
    for i, j in zip(x, y):
        ax.annotate(str(i), (i, j), textcoords="offset points", xytext=(0, 10), ha='center')
    #Save
    fig.set_figwidth(15)
    fig.set_figheight(5) 
    pp.savefig('testfig', dpi=400)

graph2sp3()
import globaldata as gd

data = gd.GlobalData()

def ListFullElectrons():
    for name in data.mElemNames:
        number = data.GetElementAtt(name, 'AtomicNumber')
        g = data.GetElementAtt(name, 'Group')[0]
        p = data.GetElementAtt(name, 'Period')[0]
        if g != '':
            print(f"{number}: {data.GetFullElecCount(int(g),int(p))}")
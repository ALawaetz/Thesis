import os
import pandas as pd
import sys
import glob

thefolder = sys.argv[1]

for file in glob.glob(thefolder + '/*.txt'):
    df = pd.read_csv(file, sep = '\t', header = None)
    names = []
    for item in df[8]:
        alist = item.split(';')
        z = 0
        for a in alist:
            if a.find('Name=') == 0:
                names.append(a[5:])
                z += 1
                break
            else:
                pass
        if z == 0:
            for a in alist:
                if a.find('name=') == 0:
                    names.append(a[5:])
                    z += 1
                    break
                else:
                    pass
        if z == 0:
            for a in alist:
                if a.find('ID=') == 0:
                    names.append(a[3:])
                    z += 1
                    break
                else:
                    pass
        if z == 0:
            for a in alist:
                if a.find('symbol=') == 0:
                    names.append(a[7:])
                    z += 1
                    break
                else:
                    pass
        if z == 0:
            for a in alist:
                if a.find('Description=') == 0:
                    names.append(a[12:])
                    z += 1
                    break
                else:
                    pass
        if z == 0:
            for a in alist:
                if a.find('product=') == 0:
                    names.append(a[8:])
                    z += 1
                    break
                else:
                    pass
        if z == 0:
            names.append('No_name')
        else:
            pass

    newdf = pd.DataFrame()
    newdf['names'] = names
    newdf['value'] = df[9]
    newdf.to_csv(file, sep = '\t', index = False, header = None)

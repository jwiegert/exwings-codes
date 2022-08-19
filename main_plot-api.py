# A local API for more easily explore the data within r3dresults.
#
# 20220817: 
# As of now, my results are only from various tests. But in the future there
# will be a lot of data, a lot of simulated images and SEDs, all in various
# viewing angles, resolutions etc etc. For this reason I need a much more
# convenient way of looking at these. As such I designed this simple API that I
# only run locally, to much faster look at results and inter-compare results.
#



#    API't ska automatiskt lista alla kataloger i r3dresults, sen alla i den man väljer där, och sen alla sed-filer, spectrum*.out och alla bildfiler, image*.out
#    Patherna innehåller info om modell, fas, inclination och våglängder
#    Plotta alla bilder med samma inklination, plotta alla bilder med samma våglängd
#    Plotta alla SEDer med samma inklination i samma figur
#        os.listdir(path='../r3dresults/st28gm06n056/140/')    
#    eller
#        import glob
#        glob.glob('../r3dresults/st28gm06n056/140/image*.out')
#    och för kataloger
#        glob.glob('../r3dresults/st28gm06n056/*/')
#    Ger då listor som kan vara inputs till mina menyer


import os

path='../r3dresults/'


# List all folders in r3dresults:
modelnames = [modelname for modelname in os.listdir(path=path) if os.path.isdir(path+modelname) == True]
modelnames.sort()

# List all folders in chosen folder
# TODO chose folder, done in API...
#      this is a test choice
modelname = modelnames[-1]+'/'

phases = [phase for phase in os.listdir(path=path+modelname) if os.path.isdir(path+modelname+phase) == True]


#for modelname in modelnames:
#    phases.append(os.listdir(path=f'../r3dresults/{modelname}'))

print(modelnames)
print()
print(phases)



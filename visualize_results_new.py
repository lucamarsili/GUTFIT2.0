from gutfit.experimentalneutrinomassmatrix import ExperimentalNeutrinoMassMatrix
from gutfit.type1and2seesaw_v4 import Type1And2SeeSaw_v4, Predictions
import examples.multinest_v4
from gutfit import parameterlist
from pymultinest import watch, analyse, plot


a = Type1And2SeeSaw_v4()

class Visualize(watch.ProgressWatcher):
    import numpy as np
    import shutil,os

    def load_parameters():
        x = np.loadtxt(self.live) #look at end of multinest and see if you can make it more low level
        return x

    def evaluate():
        return P

#build OBS as a daughter class of Type... and then use the trick with sample

import optparse, os, sys
import numpy as np

op = optparse.OptionParser(usage=__doc__)
opts, args = op.parse_args() ###look at it!!!!


def sample(x,E,S,PL,N,pnames):
    D = []
    T = []
    trials=0
    for _ in range(N):
        p=PL()
#        from IPython import embed
#        embed()
#        exit(0)
        for num, pn in enumerate(pnames):
            p[pn] = x[num]
        data = E(p)
#        foundgoodpoint=False
#        while not foundgoodpoint:
#            if trials>1000000:
#                return D,T,False
#            p=PL()
##            print(pnames, len(pnames))
##            from IPython import embed
##            embed()
##            exit(1)
#            for num, pn in enumerate(pnames):
#                p[pn] = x[num]
#            if S.constraint(p):
#                foundgoodpoint = True
#            trials+=1
        theo = S(p)
       # print(theo)
        D.append(data)
        T.append(theo)
    return D, T, True



def load_parameter2():
    points = np.loadtxt("%s/GUTFIT.txt"%args[0])
   # N = len(points[0,:]);
    likelihood = points[:,1]
    Parameters = points[:,2:]
    #print(x[-1])
    return Parameters, likelihood, points

if __name__ == "__main__":
    PL = parameterlist.ParameterList.fromConfigFile(args[1])#"examples/param_card.dat")
    Parameters, likelihood, points = load_parameter2()
   # x = Visualize.evaluate(parameters)

    O = Predictions()
    E = ExperimentalNeutrinoMassMatrix()

    usethese = []
    Observations = []
    bounds, pnames = PL.getBox(usethese)
    for i in range(0,len(likelihood)):
        PP = [Parameters[i,j] for j in range(len(pnames))]
        #print(PP)
        D, T, isgood = sample(PP,E, O, PL, 1, pnames)
        Observations.append(T[0])
        #print(T[0])

    Observations = np.asarray(Observations)
    #print(Observations)[:,-1]
   # print(Observations)
    np.savetxt('%s/GUTFITpredictions.txt'%args[0], Observations)
    np.savetxt('%s/GUTFITpara.txt'%args[0],points)
    
    #plot
    import matplotlib
    import matplotlib.pyplot as plt
    x = Parameters[:,3]
    #print(x)
    y = Parameters[:,4]
    z = -likelihood
    
    fig, ax = plt.subplots(1,1,figsize=(10,10))
    normalize = matplotlib.colors.Normalize(vmin=z.min(), vmax= z.max())
    plot = ax.scatter(x,y,c =z,s = 10, cmap = plt.cm.PuBu, norm = normalize, marker= "o") 
    ax.set(xlabel = r"$r1$", ylabel= r"$r_2$")

    cax1 = fig.add_axes([0.355, .03, 0.3, 0.05])
    cbar1= fig.colorbar(plot,cax1, orientation = 'horizontal')
    #cax2 = fig.add_axes([0.590, .48, 0.3, 0.02])
   # cbar2= fig.colorbar(plot2,cax2, orientation = 'horizontal')
   
    #cbar1.set_label(label = r"$\chi^2$", fontsize=30, weight='bold') 
    cbar1.ax.tick_params(labelsize='large')
    ax.set_xlim([-0.03, 0.03])
    ax.set_ylim([-4,-2.5])

    plt.savefig("r1r2paperrange500.pdf")
    print("done")

    
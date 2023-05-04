from gutfit.experimentalneutrinomassmatrix import ExperimentalNeutrinoMassMatrix
from gutfit.type1and2seesaw_v4_sign_lin import Type1And2SeeSaw_v4
from gutfit import parameterlist
#Try to use the battacharruia measure in order to confront two distribution

#Parameters from 1octant non susy, different from param_card3sigma which are second octant

'''
Ye = 2.7e-6
Ymu = 5.71e-4
Ytau = 9.7e-3
sigma_e = 0.023e-6
sigma_mu = 0.00557e-4
sigma_tau = 0.0727e-3
theta12 = 0.583516667
theta13 = 0.150371111
theta23 = 0.734411111
dm21 = 7.42e-5
dm31 = 2.514e-3
sigmatheta12 = 0.013432222
sigmatheta13 = 0.002093333
sigmatheta23 = 0.019188889
sigmadm21 = 0.21e-5
sigmadm31 = 0.028e-3

'''

#Data = [theta13, theta12, theta23, dm21, dm31,Ye, Ymu, Ytau] 
#Uncertainties = [sigmatheta13, sigmatheta12,sigmatheta23,sigmadm21,sigmadm31,sigma_e,sigma_mu,sigma_tau]
#Best fit values from mathematica



#Current measure: similar to chi squared, fix uncertainties in param_card3sigma_v4.dat

#mR logarithmic

#ce,cnu,r1,r2 linear, check always the consistency with parameter in param_card3sigma_v4.dat before running



def plotting(D, T, measure, fout): #Plot the observable distribution for the best fit point
    from matplotlib import pyplot as plt
    import matplotlib as mpl
    plt.style.use('ggplot')

    plt.clf()


    fig, axes = plt.subplots(figsize=(80, 10), sharex=False, sharey=False, ncols=8, nrows=1)
    plt.title("Measure: {}".format(measure))
    for i in range(8):
                values = [d[i] for d in D]
                theos  = [t[i] for t in T]
                axes[i].hist(theos, bins=1)
                axes[i].hist(values, bins=1)
      
    axes[0].set_title(r"$\theta_{13}$")
    axes[1].set_title(r"$\theta_{12}$")
    axes[2].set_title(r"$\theta_{23}$")
    axes[3].set_title(r"$\Delta m^2_{21}$")
    axes[4].set_title(r"$\Delta m^2_{31}$")
    axes[5].set_title(r"$m_e$")
   # axes[5].set_xlim([0.0])
    axes[6].set_title(r"$m_\mu$")
    axes[7].set_title(r"$m_\tau$")
    print(measure)
        
    plt.savefig(fout)
    
def readparameter(fname):  #reads the parameter card 
    d = {}
    with open(fname) as f:
        for line in f:
            l = line.strip()
            if len(l) == 0 or l.startswith("#"):
                continue
                print(l)

            k,v  = l.split(" ")
            d[k] = float(v)
    return d
    

def overlap(A,B): #Unused method at the moment
    A.sort()
    B.sort()

    minA = min(A)
    minB = min(B)
    maxA = max(A)
    maxB = max(B)

    Ni = 0
    if   maxA < minB: return 1000
    elif maxB < minA: return 1000
    elif minA < minB and maxA > maxB: return 0
    elif minA > minB and maxA < maxB: return 0


    elif minA < minB and maxA < maxB:
        i,j = 0,0
        while(A[i]<minB): i+=1
        while(B[j]<maxA): j+=1
        return (len(A)-i + j)/float((len(A)+len(B)))

    elif minA > minB and maxA > maxB:
        i,j = 0,0
        while(A[i]<maxB): i+=1
        while(B[j]<minA): j+=1
        return (len(A) + i - j)/float((len(A)+len(B)))

    else:
        print("whoopsiedoodles")


def sample(x,E,S,PL,N,pnames): #Creates a distribution with N points for theory and data 
    D = []
    T = []
    trials=0
    for _ in range(N):
        p=PL()
        for num, pn in enumerate(pnames):
            p[pn] = x[num]
        data = E(p)
        theo = S(p)
        D.append(data)
        T.append(theo)
    return D, T, True

#Current measure use only DS and DT, similar to chi square
def measure(D, T): 
  
    DS =[   abs(np.mean([t[0] for t in T]) - [d[0] for d in D][0]),
	    abs(np.mean([t[1] for t in T]) - [d[1] for d in D][0]),
    	    abs(np.mean([t[2] for t in T]) - [d[2] for d in D][0]),
    	    abs(np.mean([t[3] for t in T]) - [d[3] for d in D][0]),
    	    abs(np.mean([t[4] for t in T]) - [d[4] for d in D][0]),
    	    abs(np.mean([t[5] for t in T]) - [d[5] for d in D][0]),
    	    abs(np.mean([t[6] for t in T]) - [d[6] for d in D][0]),
    	    abs(np.mean([t[7] for t in T]) - [d[7] for d in D][0])
    ]

    #FIND a WAY TO IMPLEMENT UNCERTAINTIES DEBUG

    UNC = [np.mean([d[8] for d in D]),
           np.mean([d[9] for d in D]),
           np.mean([d[10] for d in D]),
           np.mean([d[11] for d in D]),
           np.mean([d[12] for d in D]),
           np.mean([d[13] for d in D]),
           np.mean([d[14] for d in D]),
           np.mean([d[15] for d in D]),

    ]

    from functools import reduce
    return sum([(d*d)/(i*i) for d,i in zip(DS, UNC)]) # works technically


#discrete points just activate again the rephasing, multimodal mode of multinest should distinguish the different regions if we are interested in the signs

#https://stackoverflow.com/questions/13377046/scipy-fill-a-histogram-reading-from-a-db-event-by-event-in-a-loop
import numpy as np

if __name__=="__main__":

    import optparse, os, sys
    op = optparse.OptionParser(usage=__doc__)
    op.add_option("-o", "--output",    dest="OUTPUT",      default="nestout", type=str, help="Prefix for outputs (default: %default)")
    op.add_option("-v", "--debug",     dest="DEBUG",       default=False, action="store_true", help="Turn on some debug messages")
    op.add_option("-q", "--quiet",     dest="QUIET",       default=False, action="store_true", help="Turn off messages")
    op.add_option("--mn-seed",         dest="SEED",        default=-1, type=int,              help="Multinest seed (default: %default)")
    op.add_option("--mn-resume",       dest="RESUME",      default=False, action='store_true', help="Resume on previous run.")
    op.add_option("--mn-multi-mod",    dest="MULTIMODE",   default=False, action='store_true', help="Set multimodal to true.")
    op.add_option("--mn-update",       dest="UPDATE",      default=1000, type=int, help="Update inteval (default: %default iterations)")
    op.add_option("--mn-tol",          dest="TOL",         default=0.5, type=float, help="Evidence tolerance (default: %default)")
    op.add_option("--mn-eff",          dest="EFF",         default=0.8, type=float, help="Sampling efficiency (default: %default)")
    op.add_option("--mn-points",       dest="POINTS",      default=400, type=int,              help="Number of live points in PyMultinest (default: %default)")
    op.add_option("--mn-imax",         dest="ITMAX",       default=0, type=int, help="Max number of iterations PyMultinest, 0 is infinite (default: %default)")
    op.add_option("--mn-multimodal",   dest="MULTIMODAL",  default=False, action='store_true', help="Run in multimodal mode.")
    op.add_option("--mn-no-importance",dest="NOIMPSAMP",   default=False, action='store_true', help="Turn off importance sampling.")

    opts, args = op.parse_args()

    
    # Output directory
    rank=0
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.Get_rank()
    except Exception as e:
        print("MPI (mpi4py) is not available (try e.g. pip install mpi4py), proceeding serially:", e)
        size = 1

    try:
        os.makedirs(opts.OUTPUT)
    except:
        pass
    




    import sys
    E  = ExperimentalNeutrinoMassMatrix()
    S  = Type1And2SeeSaw_v4()
    PL = parameterlist.ParameterList.fromConfigFile(args[0])
    N  = 1 # number of samples of PL
    
    if len(args) > 1:
        pdict        =   readparameter(args[1])
        pnames       = [ k for k,v in pdict.items() ]
        D, T, DUMMY  = sample([pdict[pn] for pn in pnames] ,E,S,PL,10*N, pnames)
        valuemeasure = measure(D,T)
        print(valuemeasure)
        #plotting(D, T, valuemeasure, args[2])
    else:
        
    
        

        usethese = [] #it select free parameters, look at parameterlist.py for understand better what happens 

        bounds, pnames = PL.getBox(usethese)

        PMIN   = [b[0] for b in bounds]
        PMAX   = [b[1] for b in bounds]
        PLEN   = [PMAX[i] - PMIN[i] for i in range(len(pnames))]

        def scaleParam(p, idx):
            return PMIN[idx] + p * PLEN[idx]
    # here the constraint should be applied, i.e. rather than sampling a box in the theory space, a "line/hypersurface" which satisfies the constraint should be sampled from
        def myprior(cube, ndim, nparams): #prior definition, after a rough general scan one can updates the prior
            for i in range(ndim):
                cube[i] = scaleParam(cube[i], i)

        def loglike(cube, ndim, nparams): #propose a gaussian
            PP=[cube[j] for j in range(ndim)]
            
            DATA, THEORY, goodpoint = sample(PP, E, S, PL, N, pnames)
            if not goodpoint:
                print("give up")
                return -1e101
            val = measure(DATA, THEORY)
            print(measure(DATA, THEORY))
            if val == 0:
                return -1e101


            loglikelihood = -val # Ad-hoc
            return loglikelihood


        import pymultinest

        pymultinest.run(loglike, myprior, len(pnames),
                importance_nested_sampling = not opts.NOIMPSAMP,
                verbose = False if opts.QUIET else True,
                multimodal=opts.MULTIMODAL,
                resume=opts.RESUME,
                n_iter_before_update=opts.UPDATE,
                evidence_tolerance=opts.TOL,
                sampling_efficiency = opts.EFF,
                init_MPI=False,
                n_live_points = opts.POINTS,
                max_iter=opts.ITMAX,
                seed=opts.SEED,
                outputfiles_basename='%s/GUTFIT'%opts.OUTPUT)

        if (rank ==0):
        
                import json
                json.dump(pnames, open('%s/GUTFITparams.json'%opts.OUTPUT, 'w'))
                json.dump(pnames, open('%s/params.json'%opts.OUTPUT, 'w'))

                NP = len(pnames)
                print("Now analyzing output from {}/GUTFIT.txt".format(opts.OUTPUT))
        
        
                a = pymultinest.Analyzer(n_params = NP, outputfiles_basename='%s/GUTFIT'%opts.OUTPUT)
                a.get_data()
                try:
                    s = a.get_stats()
                except:
                    print("There was an error accumulating statistics. Try increasing the number of iterations, e.g. --mn-iterations -1")
                    sys.exit(1)

                from collections import OrderedDict
                resraw = a.get_best_fit()["parameters"]
                D, T, DUMMY  = sample(resraw ,E,S,PL,10*N, pnames)
                bestval=measure(D,T)
                PP=OrderedDict.fromkeys(pnames)
                for num, pname in enumerate(pnames): PP[pname] = resraw[num]
                out="# Best fit point (measure: {}):\n".format(bestval)
                for k in PP: out+= "%s %.16f\n"%(k,PP[k])
                with open("%sconfig.best"%a.outputfiles_basename, "w") as f: f.write(out)
                print(out)

                print("Measure at best fit point is {}".format(bestval))
        
                plotting(D, T, bestval, "Bestfitoverlap.pdf")
        else: 
            print("Altri rank")
        


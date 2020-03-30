import sys
import numpy as np
import inifile
import mddm
import emcee
from getdist import plots,MCSamples

class PostProcess():

    def __init__(self,chains,chainlabels,paramnames,paramlabels,dparamnames,dparamlabels):
        self.BG = mddm.Background(30,10,1e-8,verbose=0)
        self.chains = chains
        self.chainlabels = chainlabels
        self.paramnames = paramnames
        self.paramlabels = paramlabels
        self.dparamnames = dparamnames
        self.dparamlabels = dparamlabels
        self.getdist_samples = []
        
    def ReadChains(self):
        
        for i in range(len(self.chains)):
            fname = self.chains[i]
            label = self.chainlabels[i]
            
            print("\n#chain name:",fname)

            Ini = inifile.IniFile(self.chains[i].replace(".h5","_params.ini"))
            section = "MCMC"
            paramranges = {}
            for key in self.paramnames:
                paramranges[key] = Ini.ReadFloatArray(section,key)[0:2]
            reader = emcee.backends.HDFBackend(fname,read_only=True)
            #tau = reader.get_autocorr_time()
            
            samples = reader.get_chain(flat=True)
            log_prob_samples = reader.get_log_prob(flat=True)
            weights = np.array(log_prob_samples>-np.inf).astype(np.float)
            blob_samples = reader.get_blobs(flat=True)
            blobnames = blob_samples.dtype.names

            flag = False
            for dname in self.dparamnames:
                if not dname in blobnames:
                    print("error: not found a derived parameter ",dname)
                    flag = True
            if(flag): sys.exit(1)
            derived_samples = np.array([blob_samples[dname] for dname in self.dparamnames]).transpose()

            a = np.nanmin(derived_samples,axis=0)
            b = np.nanmax(derived_samples,axis=0)
            dparamranges = np.concatenate([a[:,None],b[:,None]],axis=1)
            
            print(" flat chain shape: {0}".format(samples.shape))
            print(" flat log prob shape: {0}".format(log_prob_samples.shape))
            print(" flat blob shape: {0}".format(blob_samples.shape))
            print(" flat derived chain shape: {0}".format(derived_samples.shape))

            paramnames = []
            paramnames.extend(self.paramnames)
            paramnames.extend(self.dparamnames)
            paramlabels = []
            paramlabels.extend(self.paramlabels)
            paramlabels.extend(self.dparamlabels)
            for j in range(len(self.dparamnames)):
                paramranges[self.dparamnames[j]] = dparamranges[j,:]

            self.getdist_samples.append(
                MCSamples(samples=np.concatenate([samples,derived_samples],axis=1),
                          names=paramnames,loglikes=-log_prob_samples,
                          labels=paramlabels,label=label,ranges=paramranges,weights=weights)
            )

            # derived parameters?
            #print("###")
            #pidx = {}
            #for p in ("ob", "odm", "ol", "decay_rate", "mratio", "nnu", "mnu"):
            #    try:
            #        pidx[p] = self.paramnames.index(p)
            #    except ValueError:
            #        pass
            #if("odm" in self.paramnames and "ol" in self.paramnames):
            #    self.getdist_samples[i].addDerived(np.sqrt(samples[:,pidx["odm"]]+samples[:,pidx["ol"]])*100,name="h_nodecay",label="H_\emptyset",range=None)
            #if("decay_rate" in self.paramnames):
            #    self.getdist_samples[i].addDerived(samples[:,pidx["decay_rate"]],name="lt",label="\\tau {\\rm [Gyr]}",range=None)

            #self.getdist_samples[i].getConvergeTests()

            tau = np.array([self.getdist_samples[i].getCorrelationLength(j) for j in range(len(self.paramnames))])
            burnin = int(2.*np.max(tau))
            thin = max(int(0.5*np.min(tau)),1)
            print(" correlation lengthes: {}".format(tau))
            print(" burn-in: {0}".format(burnin))
            print(" thin: {0}".format(thin))
            print(" burn-in fraction: {0}".format(burnin/len(log_prob_samples)))
            print(" thinned sample number: {0}".format(int(len(log_prob_samples)/thin)))
            self.getdist_samples[i].deleteZeros()
            self.getdist_samples[i].removeBurn(burnin)
            self.getdist_samples[i].thin(thin)

            # constraints
            for j, mean in enumerate(self.getdist_samples[i].getMeans()):
                #print(" %s = %f +/- %f +/- %f" % (self.getdist_samples[i].parLabel(j),\
                #                                  mean, mean - self.getdist_samples[i].confidence(j,upper=True,limfrac=0.16),\
                #                                  mean - self.getdist_samples[i].confidence(j,upper=False,limfrac=0.16)) )
                print(" %s mean = %f" % (self.getdist_samples[i].parLabel(j),mean))
                print("  two-tail error 68 C.L. = +%f, -%f" %(self.getdist_samples[i].confidence(j,upper=True,limfrac=0.16)-mean,\
                                                                  mean-self.getdist_samples[i].confidence(j,upper=False,limfrac=0.16)))
                upper = self.getdist_samples[i].confidence(j,upper=True,limfrac=0.05)
                print("  one-tail upper limit 95 C.L. = %f" % upper)
                lower = self.getdist_samples[i].confidence(j, upper=False, limfrac=0.05)
                print("  one-tail lower limit 95 C.L. = %f" % lower)

        

def main():
    args = sys.argv
    if(len(args)<2):
        print("error: the number of input parameters is not correct; input command must be")
        print(">$ python post.py input_file_name[e.g. dist.ini]")
        sys.exit(1)

    Ini = inifile.IniFile(args[1])
    Ini.Dump()

    section = "POSTPROCESS"
    postroot = Ini.ReadString(section,"postroot")
    Ini.Dump(postroot+"_dist.ini")

    chains = Ini.ReadString(section,"chains").split(",")
    chainlabels = Ini.ReadString(section,"chainlabels").split(",")
    names = Ini.ReadString(section,"paramnames").split(",")
    labels = Ini.ReadString(section,"paramlabels").split(",")
    dnames = Ini.ReadString(section,"dparamnames").split(",")
    dlabels = Ini.ReadString(section,"dparamlabels").split(",")

    PP = PostProcess(chains,chainlabels,names,labels,dnames,dlabels)
    PP.ReadChains()

    # triangle plots
    g = plots.getSubplotPlotter()
    g.triangle_plot(PP.getdist_samples,filled=True)
    g.export(postroot+"_triangle.pdf")
    
main()

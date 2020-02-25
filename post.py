import sys
import numpy as np
import inifile
import mddm
import emcee
from getdist import plots,MCSamples

class PostProcess():

    def __init__(self,chains,chainlabels,paramnames,paramlabels):
        self.BG = mddm.Background(30,10,1e-8,False,verbose=0)
        self.chains = chains
        self.chainlabels = chainlabels
        self.paramnames = paramnames
        self.paramlabels = paramlabels
        self.getdist_samples = []
        
    def ReadChains(self):
        
        for i in range(len(self.chains)):
            fname = self.chains[i]
            label = self.chains[i]
            
            print("\n#chain name:",fname)
            reader = emcee.backends.HDFBackend(fname,read_only=True)
            
            tau = reader.get_autocorr_time()
            burnin = int(2 * np.max(tau))
            thin = int(0.5 * np.min(tau))
            samples = reader.get_chain(discard=burnin,flat=True,thin=thin)
            log_prob_samples = reader.get_log_prob(discard=burnin,flat=True,thin=thin)
            #log_prior_samples = reader.get_blobs(discard=burnin,flat=True,thin=thin)
        
            print(" burn-in: {0}".format(burnin))
            print(" thin: {0}".format(thin))
            print(" flat chain shape: {0}".format(samples.shape))
            print(" flat log prob shape: {0}".format(log_prob_samples.shape))
            #print(" flat log prior shape: {0}".format(log_prior_samples.shape))
            
            all_samples = np.concatenate((log_prob_samples[:,None],samples),axis=1)
            
            self.getdist_samples.append(
                MCSamples(samples=all_samples[:,1:],names=self.paramnames[1:],labels=self.paramlabels[1:],label=label)
            )

            # constraints
            for j, mean in enumerate(self.getdist_samples[i].getMeans()):
                print(" %s = %f +/- %f +/- %f" % (self.getdist_samples[i].parLabel(j),\
                                                  mean, mean - self.getdist_samples[i].confidence(j, limfrac=0.16),\
                                                  mean - self.getdist_samples[i].confidence(j, limfrac=0.025)) )
                upper = self.getdist_samples[i].confidence(j,upper=True,limfrac=0.05)
                print("  upper limit 95 C.L. = %f" % upper)
                lower = self.getdist_samples[i].confidence(i, upper=False, limfrac=0.05)
                print("  lower limit 95 C.L. = %f" % lower)

        

def main():
    args = sys.argv
    if(len(args)<2):
        print("error: the number of input parameters is not correct; input command must be")
        print(">$ python post.py input_file_name[e.g. dist.ini]")
        sys.exit(1)

    Ini = inifile.IniFile(args[1])
    Ini.Dump()


    section = "POSTPROCESS"
    chains = Ini.ReadString(section,"chains").split(",")
    chainlabels = Ini.ReadString(section,"chainlabels").split(",")
    names = Ini.ReadString(section,"paramnames").split(",")
    names.insert(0,"lnlike")
    labels = Ini.ReadString(section,"paramlabels").split(",")
    labels.insert(0,"\ln P")
    PP = PostProcess(chains,chainlabels,names,labels)
    PP.ReadChains()

    # triangle plots
    g = plots.getSubplotPlotter()
    g.triangle_plot(PP.getdist_samples[:],filled=True)
    g.export(Ini.ReadString(section,"postroot")+"_triangle.pdf")
    
main()

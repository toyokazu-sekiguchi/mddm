import sys
import const
import mddm
import numpy as np
import likelihoods
import inifile
import mcmc

def main():
    args = sys.argv
    if(len(args)<2):
        print("error: the number of input parameters is not correct; input command must be")
        print(">$ python inifiles.py input_file_name")
        sys.exit(1)

    Ini = inifile.IniFile(args[1])
    Ini.Dump()

    section = "OUTPUT"
    root = Ini.ReadString(section,"root")
    
    # model calculation with fiducial parameters
    section = "FIDUCIAL COSMOLOGY"
    paramsfid = np.array([Ini.ReadFloat(section,"ob"),Ini.ReadFloat(section,"odm"),Ini.ReadFloat(section,"ol"),
                          Ini.ReadFloat(section,"decay_rate"),Ini.ReadFloat(section,"mratio"),
                          Ini.ReadFloat(section,"nnu"),Ini.ReadFloat(section,"mnu")])
    section = "DDM SETUP"
    BG = mddm.Background(Ini.ReadInt(section,"num_a"),Ini.ReadInt(section,"max_it"),
                         Ini.ReadFloat(section,"tol"),verbose=1)
    BG.SetParams(paramsfid)
    BG.Iterate()
    
    # parameters for output scale factors
    #BG.Output(np.logspace(np.log10(1e-3),np.log10(1),50,base=10))
    
    # likelihood calculation
    section = "LIKELIHOODS"
    LF = likelihoods.Likelihood(Ini.ReadBoolean(section,"use_BAO"),Ini.ReadBoolean(section,"use_H0"),
                                Ini.ReadBoolean(section,"use_CMB"),verbose=1)
    lnL = LF.LnLike(BG)
    print(" ln(L)=",lnL[:])
    
    # MCMC
    #inp = input("\nDoes everything seem going fine? Then let's run MCMC [Y/n]:")
    #if(inp=='n'):
    #    sys.exit()
    section ='MCMC'
    BG.verbose  = LF.verbose  = 0
    paramrange = [Ini.ReadFloatArray(section,"ob"),Ini.ReadFloatArray(section,"odm"),Ini.ReadFloatArray(section,"ol"),
                  Ini.ReadFloatArray(section,"decay_rate"),Ini.ReadFloatArray(section,"mratio"),
                  Ini.ReadFloatArray(section,"nnu"),Ini.ReadFloatArray(section,"mnu")]
    map_varied = np.array([i for i in range(len(paramrange)) if len(paramrange[i]>0)])
    range_varied = np.array([paramrange[i] for i in map_varied])
    driver = mcmc.MCMC(paramsfid,Ini.ReadBoolean(section,"parallel"))
    driver.SetParams(map_varied,range_varied,Ini.ReadInt(section,"nwalkers"),verbose=1)
    chain = root+".h5"
    driver.Run(chain,Ini.ReadInt(section,"nsteps"),BG,LF)
    
main()

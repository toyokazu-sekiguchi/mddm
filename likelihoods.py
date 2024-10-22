import const
import mddm
import numpy as np

class DataBAO:

    def __init__(self,verbose=0):
        self.verbose = verbose
        
        # Alam et al. 2017 (arXiv:1607.03155)
        self.nz0 = 3
        self.nzobs0 = 2
        self.rdrag0 = 147.78*const.Mpc
        nobs = self.nz0*self.nzobs0
        ncorr = nobs*(nobs+1)
        self.z0 = np.array([0.38,0.51,0.61])
        self.mean0 = np.array([1518,81.5,1977,90.4,2283,97.3])
        err = np.array([22,1.9,27,1.9,32,2.1])
        corr = np.array([1.0000,0.2280,1.0000,0.4970,0.1536,1.0000,0.1117,0.4873,0.2326,
                         1.0000,0.1991,0.0984,0.5120,0.1571,1.0000,0.0520,0.2307,0.1211,0.5449,0.2408,1.0000])
        covmat = np.empty([nobs,nobs])
        for i in range(nobs):
            for j in range(i+1):
                ij = i*(i+1)//2+j
                covmat[i,j] = err[i]*err[j]*corr[ij]
                covmat[j,i] = covmat[i,j]
        self.invcov0 = np.linalg.inv(covmat)
        
        # Zarrouk et al. 2018 (arXiv:1801.03062)
        self.nz1 = 1
        self.nzobs1 = 2
        nobs = self.nz1*self.nzobs1
        ncorr = nobs*(nobs+1)
        self.z1 = np.array([1.52])
        self.mean1 = np.array([23.5e3,12.58])
        err = np.array([1.78e3,0.68])
        corr = np.array([1.0000,0.05,1.0000]) # sign of off diagonal term is flipped because Hr is inversely proportional to alpha_par
        covmat = np.empty([nobs,nobs])
        for i in range(nobs):
            for j in range(i+1):
                ij = i*(i+1)//2+j
                covmat[i,j] = err[i]*err[j]*corr[ij]
                covmat[j,i] = covmat[i,j]
        self.invcov1 = np.linalg.inv(covmat)
        
        # Bautista et al. 2017 (arXiv:1702.00176)
        self.nz2 = 1
        self.nzobs2 = 2
        nobs = self.nz1*self.nzobs1
        ncorr = nobs*(nobs+1)
        self.z2 = np.array([2.33])
        self.mean2 = np.array([9.07,37.77])
        err = np.array([0.31,2.13])
        corr = np.array([1.0000,0,1.0000]) # assumed diagonal
        covmat = np.empty([nobs,nobs])
        for i in range(nobs):
            for j in range(i+1):
                ij = i*(i+1)//2+j
                covmat[i,j] = err[i]*err[j]*corr[ij]
                covmat[j,i] = covmat[i,j]
        self.invcov2 = np.linalg.inv(covmat)

    def LnLike(self,BG):
        adrag = 1/(1+BG.zdrag)
        rdrag = BG.SoundHorizon(adrag)
        if(self.verbose>0):
            print(' z_drag:%f, r_drag[Mpc]:%f' %(BG.zdrag,rdrag/const.Mpc))
            
        lnl=0

        f0 = rdrag/self.rdrag0
        obs0 = np.empty(self.nz0*self.nzobs0)
        for i in range(self.nz0):
            obs0[2*i] = BG.DeltaTau(1/(1+self.z0[i]),1)/f0*const.c/const.Mpc # [Mpc]
            obs0[2*i+1] = (1+self.z0[i])**2/BG.dtauda(1/(1+self.z0[i]))*f0*const.Mpc/const.km # in [km/sec/Mpc]
        obs0 -= self.mean0
        lnl += np.dot(obs0,np.dot(self.invcov0,obs0))

        obs1 = np.empty(self.nz1*self.nzobs1)
        for i in range(self.nz1):
            obs1[2*i] = (1+self.z1[i])**2/BG.dtauda(1/(1+self.z1[i]))*rdrag/const.km # in [km/sec]
            obs1[2*i+1] = BG.DeltaTau(1/(1+self.z1[i]),1)/rdrag*const.c/(1+self.z1[i])
        obs1 -= self.mean1
        lnl += np.dot(obs1,np.dot(self.invcov1,obs1))

        obs2 = np.empty(self.nz2*self.nzobs2)
        for i in range(self.nz2):
            obs2[2*i] = (1+self.z2[i])**2/BG.dtauda(1/(1+self.z2[i]))*rdrag
            obs2[2*i] = const.c/obs2[2*i]
            obs2[2*i+1] = BG.DeltaTau(1/(1+self.z2[i]),1)/rdrag*const.c
        obs2 -= self.mean2
        lnl += np.dot(obs2,np.dot(self.invcov2,obs2))

        return -0.5*lnl

class DataH0:

    def __init__(self,verbose=0):
        self.verbose = verbose
        # Riess 2020 (arXiv:2001.03624)
        self.mean = 73.8
        self.error = 1.0

    def LnLike(self,BG):
        if(self.verbose>0):
            print(' H_0[km/sec/Mpc]:%f' %BG.H0())
        return -0.5*((BG.H0()-self.mean)/self.error)**2

class DataCMB:

    def __init__(self,verbose):
        self.verbose = verbose
        # model-independent prior on theta_MC according to Planck 2018 results VI (arXiv:1807.06209)
        self.mean = 0.010409
        self.error = 6e-6

    def LnLike(self,BG):
        astar = 1/(1+BG.zstar)
        rstar = BG.SoundHorizon(astar)
        DA = BG.DeltaTau(astar,1)*const.c
        theta = rstar/DA
        if(self.verbose>0):
            print(' z_*:%f, r_*[Mpc]:%f'%(BG.zstar,rstar/const.Mpc))
            print(' 100*theta_*:%e'%(100*theta))
        return -0.5*((theta-self.mean)/self.error)**2
    
class Likelihood:

    def __init__(self,useBAO,useH0,useCMB,verbose=0):
        self.useBAO = useBAO
        self.useH0 = useH0
        self.useCMB = useCMB
        self.verbose = verbose
        
        if(self.useBAO):
            self.BAO=DataBAO(self.verbose)
        if(self.useH0):
            self.H0=DataH0(self.verbose)
        if(self.useCMB):
            self.CMB=DataCMB(self.verbose)

    def LnLike(self,BG,lnPprior):
        if(self.useBAO or self.useCMB):
            BG.UpdateTherm()
        
        if(self.verbose>0):
            print("\n# likelihoods")

        lnL = [lnPprior]
        if(self.useBAO):
            self.BAO.verbose = self.verbose
            lnL_BAO = self.BAO.LnLike(BG)
            lnL.append(lnL_BAO)
            lnL[0] += lnL_BAO

        if(self.useH0):
            self.H0.verbose = self.verbose
            lnL_H0 = self.H0.LnLike(BG)
            lnL.append(lnL_H0)
            lnL[0] += lnL_H0

        if(self.useCMB):
            self.CMB.verbose = self.verbose
            lnL_CMB = self.CMB.LnLike(BG)
            lnL.append(lnL_CMB)
            lnL[0] += lnL_CMB

        return lnL


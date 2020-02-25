import sys
import numpy as np
from scipy import integrate,interpolate,optimize,special
import const

class BBN:
    
    def __init__(self):
        filename = "./bbn.dat" # based on PArthENoPE taken from CLASS 
        data = np.loadtxt(filename)
        self.spl_yp = interpolate.interp2d(data[:,0],data[:,1],data[:,2],kind="linear") # cubic doesn't work; why?
        
    def yp(self,ob,dnnu):
        return self.spl_yp(ob,dnnu)

class MassiveNu:

    def __init__(self):
        self.nmass = 3 # number of mass eigen states
        self.nnu = const.nnu_standard
        self.mass = np.zeros(self.nmass)

        self.lammax = 1e3 # max of am/T
        self.lammin = 1e-3 # min of am/T
        
        nlam = 1000
        lnlam = np.linspace(np.log(self.lammin),np.log(self.lammax),nlam)
        lnrho = np.zeros(nlam)
        for i in range(nlam):
            lam = np.exp(lnlam[i])
            lnrho[i] = integrate.quad(lambda x:x**2*np.sqrt(x**2+lam**2)/(np.exp(x)+1),0,100)[0]
        lnrho = np.log(lnrho[:]*const.nu_energy) # ratio of massive to massless
        self.spl_lnrho = interpolate.make_interp_spline(lnlam,lnrho)
        
    def Rho(self,a):
        lam = a*self.mass[:]*const.eV/const.TCnuB
        r = np.empty(self.nmass)
        for i in range(self.nmass):
            if(lam[i]>self.lammax): # non-relativistic
                r[i] = lam[i]*const.nu_number
            elif(lam[i]<self.lammin):
                r[i] = 1
            else:
                r[i] = np.exp(self.spl_lnrho(np.log(lam[i])))
        return r
        
    def SetParams(self,hierarchy,nnu,mnu): # mnu [eV]
        self.hierarchy = hierarchy
        self.nnu = nnu
        mlower = 0
        mupper = mnu
        mlight = optimize.brentq(lambda x:self.lightest2tot(x)-mnu,mlower,mupper)
        if(self.hierarchy==1): # normal
            self.mass[0] = mlight
            self.mass[1] = np.sqrt(mlight**2+const.m2nu21)
            self.mass[2] = np.sqrt(mlight**2+const.m2nu21+const.m2nu32)
        elif(self.hierarchy==-1): # inverted
            self.mass[0] = np.sqrt(mlight**2-const.m2nu21+const.m2nu32)
            self.mass[1] = np.sqrt(mlight**2+const.m2nu32)
            self.mass[2] = mlight
        else: # degenerate
            self.mass[0:2] = mlight
        
    def lightest2tot(self,mlight):
        if(self.hierarchy==1): # normal
            return mlight+np.sqrt(mlight**2+const.m2nu21)+np.sqrt(mlight**2+const.m2nu21+const.m2nu32)
        elif(self.hierarchy==-1): # inverted
            return np.sqrt(mlight**2-const.m2nu21+const.m2nu32)+np.sqrt(mlight**2+const.m2nu32)+mlight
        else: # degenerate
            return mlight*self.nmass
        
class Background:
    
    ndim = 7
    
    def __init__(self,num_a,max_it,tol,calcthrm,verbose=0):
        self.verbose = verbose
        # default cosmological parameters
        self.og = np.pi**2/15*const.TCMB**4/(const.c*const.hbar)**3/const.rhoch2
        self.nu = MassiveNu()
        self.onu_nomass = self.og*7/8*(const.TCnuB/const.TCMB)**4*self.nu.nnu
        self.ob = 0.0224
        self.odm = 0.120
        self.ol = 0.311
        self.Gamma_in_BigH = 0
        self.epsilon = 0
        self.bbn = BBN()
        self.yp = 0.24714 # based on PRIMAT assuming ob = 0.02236, neff = 3.046

        # default binning of a
        self.num_a = num_a # >=30 usually works
        self.max_it = max_it
        self.tol = tol
        
        # recombination history
        self.calcthrm = calcthrm

        # arrays
        self.arr_a = np.empty(self.num_a)
        self.arr_Gdt = np.zeros_like(self.arr_a)
        self.arr_Gt = np.zeros_like(self.arr_a)
        self.arr_og = np.empty_like(self.arr_a)
        self.arr_onu = np.empty_like(self.arr_a)
        self.arr_ob = np.empty_like(self.arr_a)
        self.arr_odm0 = np.empty_like(self.arr_a) # decaying DM
        self.arr_odm1 = np.empty_like(self.arr_a) # massless decay product 
        self.arr_odm2 = np.empty_like(self.arr_a) # massive decay product
        self.arr_ol = np.empty_like(self.arr_a)
            
    def EarlyGt(self,a): # in matter-radiation dominated universe well before decay
        aeq = (self.og+self.onu_nomass)/(self.odm+self.ob)
        y = a/aeq
        t = 2/3*(2+(y-2)*np.sqrt(1+y))
        return self.Gamma_in_BigH*aeq**2*t/np.sqrt(self.og+self.onu_nomass)
    
    def SetParams(self,params):
        #params is an array containing [ob,odm,ol,Gamma_Gyr,mratio,nnu,mnu]
        self.ob = params[0]
        self.odm = params[1]
        self.ol = params[2]
        self.Gamma_in_BigH = params[3]/const.Gyr/const.BigH
        self.epsilon = (1-params[4]**2)/2

        self.nu.SetParams(1,params[5],params[6]) # normal hierarchy
        self.onu_nomass = self.og*7/8*(const.TCnuB/const.TCMB)**4*self.nu.nnu
        self.yp = self.bbn.yp(self.ob,self.nu.nnu-const.nnu_standard)[0]
        
        if(self.verbose>0):
            print('# cosmological parameters:')
            print(' omega_g*h^2: %e'%self.og)
            print(' omega_nu*h^2 (no mass): %e'%self.onu_nomass)
            print(' omega_dm*h^2 (no decay): %e'%self.odm)
            print(' omega_lambda*h^2: %e'%self.ol)
            print(' Gamma [Gyr^{-1}]: %e'%(self.Gamma_in_BigH*const.Gyr*const.BigH))
            print(' epsilon: %e'%self.epsilon)
            print(' neutrino mass hierarchy [1:NH,-1:IH],0:NH]:',self.nu.hierarchy)
            print(' neutrino masses [eV]:',self.nu.mass[:])
            print(' neutrino NR redshifts:',const.TCnuB/const.eV/self.nu.mass[:])
            print(' yp: %e'%self.yp)
        
        self.a_ini = np.min([(self.ob+self.odm)/self.ol,np.min(const.TCnuB/const.eV/self.nu.mass[:])*0.1])
        while (self.EarlyGt(self.a_ini)>1e-3): # this works
            self.a_ini /= 1.5
        self.a_fin = 2
        
        self.arr_a = np.logspace(np.log10(self.a_ini),np.log10(self.a_fin),self.num_a,base=10)
        self.arr_lna = np.log(self.arr_a)
        
        # initialization for iterations assuming no decay
        self.arr_og[:] = self.og/self.arr_a[:]**4
        self.arr_onu[:] = self.onu_nomass/self.arr_a[:]**4
        for i in range(self.num_a):
            self.arr_onu[i] *= np.average(self.nu.Rho(self.arr_a[i]))
        self.arr_ob[:] = self.ob/self.arr_a[:]**3
        self.arr_odm0[:] = self.odm/self.arr_a[:]**3 # no decay
        self.arr_odm1[:] = 0
        self.arr_odm2[:] = 0
        self.arr_ol[:] = self.ol

        self.arr_Gt[0] = self.EarlyGt(self.a_ini)
        if(self.verbose>0):
            print('#numerical settings for background iteratinons')
            print(' num_a: %i'%self.num_a)
            print(' a_ini: %e, Gt_ini: %e'%(self.arr_a[0],self.arr_Gt[0]))
            print(' a_fin: %e'%self.a_fin)
        
    def Update(self):
        self.arr_Gdt[:] = self.arr_a[:]*np.sqrt(
            self.arr_og[:]+self.arr_onu[:]+self.arr_ob[:]+self.arr_odm0[:]+self.arr_odm1[:]+self.arr_odm2[:]+self.arr_ol[:])
        self.arr_Gdt[:] = self.Gamma_in_BigH/self.arr_Gdt[:] # Gamma*dt/da
        self.spl_lnGdt = interpolate.make_interp_spline(self.arr_lna,np.log(self.arr_Gdt))
        for i in range(1,self.num_a):
            self.arr_Gt[i] = self.arr_Gt[i-1]+integrate.quad(
                lambda x:np.exp(self.spl_lnGdt(np.log(x))),self.arr_a[i-1],self.arr_a[i])[0]
            
        # decay products
        spl_lnf = interpolate.make_interp_spline(self.arr_lna,np.log(self.arr_Gdt[:])-self.arr_Gt[:])
        for i in range(self.num_a):
            a = self.arr_a[i]
            self.arr_odm1[i] = integrate.quad(lambda x: x/a*np.exp(spl_lnf(np.log(x))),0,a)[0]
            self.arr_odm2[i] = integrate.quad(
                lambda x: np.sqrt((x/a)**2+(1-2*self.epsilon)/self.epsilon**2)*np.exp(spl_lnf(np.log(x))),0,a)[0]
            
        # decaying DM
        self.arr_odm0[:] = self.odm/self.arr_a[:]**3
        self.arr_odm1[:] *= self.arr_odm0[:]*self.epsilon
        self.arr_odm2[:] *= self.arr_odm0[:]*self.epsilon
        self.arr_odm0[:] *= np.exp(-self.arr_Gt[:])

    def Iterate(self):
        arr_conv = np.empty(self.max_it) # this stores odm2[num_a] at each iteration and is used for convergence check
        for it in range(self.max_it):
            self.Update()

            # convergence check using massless decay-product
            arr_conv[it] = self.arr_odm1[self.num_a-1]
            if(it==0):
                continue
            
            ferr = (arr_conv[it-1]-arr_conv[it])/arr_conv[it]
            if(ferr < self.tol):
                if(self.verbose>0):
                    print('iteration converged after %d times' % it)
                self.Update()
                break
            elif(it==self.max_it-1):
                print('warning: iteration does not converge')
    
    def Output(self,arr_aout):
        spl0 = interpolate.make_interp_spline(self.arr_lna,np.log(self.arr_odm0))
        spl1 = interpolate.make_interp_spline(self.arr_lna,np.log(self.arr_odm1))
        spl2 = interpolate.make_interp_spline(self.arr_lna,np.log(self.arr_odm2))
        spl3 = interpolate.make_interp_spline(self.arr_lna,np.log(self.arr_Gt))
    
        print('')
        print('# a, omega_dm0, omega_dm1, omega_dm2, t[Gyr]')
        for aout in arr_aout:
            lnaout = np.log(aout)
            odm0 = np.exp(spl0(lnaout))
            odm1 = np.exp(spl1(lnaout))
            odm2 = np.exp(spl2(lnaout))
            t_in_Gyr = np.exp(spl3(lnaout))/(self.Gamma_in_BigH*const.BigH)/const.Gyr
            print('%e %e %e %e %e' %(aout,odm0,odm1,odm2,t_in_Gyr))

        print('')
        print('# initial density parameters')
        print('  omega_dm (no decay) = %e' % self.odm)
        print('  omega_l = %e' % self.ol)
        print('  h (no decay)= %e' % np.sqrt(self.og+self.onu_nomass*np.average(self.nu.Rho(1))+self.ob+self.odm+self.ol))
        print('# final parameters')
        odm0 = np.exp(spl0(0))
        odm1 = np.exp(spl1(0))
        odm2 = np.exp(spl2(0))
        print('  omega_dm0 = %e' % odm0)
        print('  omega_dm1 = %e' % odm1)
        print('  omega_dm2 = %e' % odm2)
        print('  h = %e' % np.sqrt(self.og+self.onu_nomass*np.average(self.nu.Rho(1))+self.ob+odm0+odm1+odm2+self.ol))

    def dtauda(self,a):
        if(a<self.a_ini):
            x = 1/np.sqrt(self.og+self.onu_nomass*np.average(self.nu.Rho(a))+(self.ob+self.odm)*a+self.ol*a**4)
        else:
            x = np.exp(self.spl_lnGdt(np.log(a)))/a/self.Gamma_in_BigH
        return x/const.BigH
        
    def H0(self):
        return np.exp(-self.spl_lnGdt(0))*self.Gamma_in_BigH*100
        
    def DeltaTau(self,a1,a2):
        return integrate.quad(self.dtauda,a1,a2)[0]

    def drsda(self,a):
        cs = np.sqrt(1/3/(1+0.75*a*self.ob/self.og))*const.c
        return cs*self.dtauda(a)
    
    def SoundHorizon(self,a):
        return integrate.quad(self.drsda,0,a)[0]

    def UpdateTherm(self):
        from HyRec import pyrec
        if(not self.calcthrm):
            sys.exit()

        if(self.verbose>0):
            print("thermal history is computed by HyRec")

        ok = 0
        w0 = -1
        wa = -1
        pyrec.rec_build_history_wrap(const.TCMB/const.kB,self.ob,self.odm,ok,self.ol,w0,wa,self.yp,self.nu.nnu,tuple(self.nu.mass))
        
        # initial guess based on Hu & Sugiyama 1996
        ob = self.ob
        om = self.ob+self.odm
        g1 = 0.0783*ob**-0.238/(1+39.5*ob**0.763)
        g2 = 0.560/(1+21.1*ob**1.81)
        self.zstar = 1048*(1+0.00124*ob**-0.738)*(1+g1*om**g2)
        b1 = 0.313*om**-0.419*(1+0.607*om**0.674)
        b2 = 0.238*om**0.223
        self.zdrag = 1345*om**0.251/(1+0.659*om**0.828)*(1+b1*ob**b2)
        
        zstar1 = self.zstar*0.9
        zstar2 = self.zstar*1.1
        zdrag1 = self.zdrag*0.9
        zdrag2 = self.zdrag*1.1
        
        akthom = ob*const.rhoch2/const.c*(1-self.yp)/const.m_H*const.sigmaT
        
        nthrm = 1250
        zlower = 0.
        zupper = 1250.
        zthrm = np.linspace(zlower,zupper,nthrm)
        tau_opt = np.zeros(nthrm)
        tau_drag = np.zeros(nthrm)
        for i in range(1,nthrm):
            tau_opt[i] = tau_opt[i-1]+integrate.quad(lambda z:
                                pyrec.hyrec_xe(1/(z+1))*akthom*self.dtauda(1/(z+1)),zthrm[i-1],zthrm[i])[0]
            tau_drag[i] = tau_drag[i-1]+integrate.quad(lambda z:
                                pyrec.hyrec_xe(1/(z+1))*akthom*self.dtauda(1/(z+1))/(0.75*self.ob/self.og/(1+z)),zthrm[i-1],zthrm[i])[0]
        spl_opt = interpolate.make_interp_spline(np.log(zthrm[1:]),np.log(tau_opt[1:]))
        spl_drag = interpolate.make_interp_spline(np.log(zthrm[1:]),np.log(tau_drag[1:]))
            
        self.zstar = optimize.brentq(lambda x:spl_opt(np.log(x)),zstar1,zstar2)
        self.zdrag = optimize.brentq(lambda x:spl_drag(np.log(x)),zdrag1,zdrag2)
        

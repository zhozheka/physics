import numpy as np
import math

def rich_error(fun, t, u, t_prev, u_prev, metrics='l2'):

    if fun.dim == 1:

        err_est_C = 0
        err_est_L = 0
        
        for s in xrange(u_prev.shape[0]):
            eT = (t[s*2] - t_prev[s] ) / 15
            eU = (u[s*2] - u_prev[s] ) / 15

            e1 = abs(eU - eT*fun.Ut( t[s*2], u[s*2], fun.params)[0] )

            err_est_C = np.max([err_est_C, e1])
            err_est_L += e1**2

        err_est_L = np.sqrt(err_est_L / u_prev.shape[0])


    elif fun.dim > 1:
        err_est_C = np.zeros(fun.dim)
        err_est_L = np.zeros(fun.dim)
        
        for i in xrange(fun.dim):
            
            for s in xrange(u_prev.shape[0]):

                eT = (t[s*2] - t_prev[s] ) / 15
                eU_i = (u[s*2][i] - u_prev[s][i] ) / 15

                e1 = abs(eU_i - eT*fun.Ut( t[s*2], u[s*2], fun.params)[i] )
            
                err_est_C[i] = np.max([err_est_C[i], e1])
                err_est_L[i] += e1**2
                
            err_est_L[i] = np.sqrt(err_est_L[i] / u_prev.shape[0])
        
        #err_est_C = np.sqrt((err_est_C**2).mean())
        #err_est_L = np.sqrt((err_est_L**2).mean())
        
    
    if metrics == 'l2':
        return err_est_L
           
    elif metrics == 'c':
        return err_est_C
    
    
def dense_grid(H):
    n = H.shape[0]
    H_new = np.empty(2*n)
    
    
    H_new[0]  = H[0] * H[0]**0.5 / (H[0]**0.5 + H[1]**0.5)
    H_new[1] = H[0] * H[1]**0.5 / (H[0]**0.5 + H[1]**0.5)
    
    for s in xrange(1, n-1):
        
        H_new[2*s] = H[s] * H[s-1]**0.25 / (H[s+1]**0.25 + H[s-1]**0.25)        
        H_new[2*s+1] = H[s] * H[s+1]**0.25 / (H[s+1]**0.25 + H[s-1]**0.25)
      
    H_new[2*n-2] = H[n-1] * H[n-2]**0.5 / (H[n-2]**0.5 + H[n-1]**0.5)
    H_new[2*n-1] = H[n-1] * H[n-1]**0.5 / (H[n-2]**0.5 + H[n-1]**0.5)
    
    return H_new


def U_real(t, params):
    
    u0 = params['u0']
    a_ = params['a_']
    l0_ = params['lambda0']
    L0 = u0 / (u0**2 - a_**2) + l0_ * np.sin(t);
    return -2 * L0 * a_**2 / (1 + np.sqrt(1 + (2 *a_*L0)**2 ) );

    
def Ut_lambda(t, u, params):
    a_ = params['a_']
    lambda0 = params['lambda0']
 
    return  np.array([-lambda0*math.cos(t)*(u**2 - a_**2)**2 / (u**2 + a_**2)])



def Ut_vdp(t, u, params):
   
    w_ = params['w_']
    sigma_ = params['sigma_']
    temp = -w_**2 * u[0] - sigma_ * (u[0]**2 - 1) * u[1] 
    return np.append( u[1], temp)



def F_l(R, params, f):
    
    t = R[0]
    u = R[1:]
    
    f_ = f(t, u, params) 
    
    Tl = 1 / (1. + np.sum(f_**2))**0.5
    Ul = f_ / (1. + np.sum(f_**2))**0.5
    
    return np.append(Tl, Ul)

class adaptive_step_solver:
    
    def __init__(self, Ut, u0=0, t0=0, cappa=1, h0=0.1, L=20, scheme='rk4', **kwargs):
        
        self.u0 = self.u = np.array(u0)
        self.t0 = self.t = t0
        self.l = 0
        self.L_ = L
        self.c1 = np.zeros(3);
        self.c2 = np.zeros(3);
        self.cappa = cappa
        self.h = h0
        self.h_ = h0
        self.params = kwargs
        self.params['u0'] = u0
        self.scheme = scheme
        self.Ut = Ut
        
        if (type(u0) == int) or (type(u0) == float):
            self.dim = 1
            
        elif type(u0) == list:
            self.dim = len(u0)
    
    def F(self, R):
        #return F_lambda(R, self.params)
        return F_l(R, f=self.Ut, params=self.params)
    
    def step(self, h=None):
        
        if self.scheme == 'rk4':
            self.rk4(h)
        
        elif self.scheme == 'bork4':
            self.bork4(h)
     
          
    def rk4(self, h=None):
        
        
        if h != None:
            self.h = h
        
        #R = {t, u0, u1, ...}
        R = np.append(self.t, self.u) #pack values
        
        
        
        k1 = self.F(R)
        k2 = self.F(R + k1 * (self.h / 2))
        k3 = self.F(R + k2 * (self.h  /2))
        k4 = self.F(R + k3 * self.h)
        
        R_prev = R
        R = R + ( k1 + k2*2 + k3*2 + k4 ) * (self.h / 6);
 

        #self.t, self.u = R #unpack values
        self.t = R[0]
        self.u = R[1:]
        
        
        self.l += self.h
        
        self.c1 = (self.F(R) - self.F(R_prev)) / self.h
        self.c2 = (self.F(R)*3 - k2*2 - k3*2 + k1) / self.h;
     
        
        if self.cappa == 1:
            c = self.c1
            
        elif self.cappa == 2:
            c = self.c2
        
        self.h = self.h_ / (1 + np.linalg.norm(c)**0.5)
        
        
    
    def reset(self):
        self.u = self.u0
        self.t = 0
        self.l = 0
        
        self.T = []
        self.U = []
        self.C1 = []
        self.C2 = []
        self.H = []
        self.L = []
        
    
    def solve(self):
        pass
    
    
    def first_stage(self,  dl_0=1., DELTA=0.03, k_start=0, k_max=10):
        
        Delta = [0]
        L=[]
        
        for k in xrange(k_start, k_max): 

            dl =  dl_0 / 2**k
            L_prev = L
            #L, H_, U_, T_ = fun.step_solve(h0=dl, step_type='auto', return_list=['L', 'H', 'U', 'T'])

            self.step_solve(h0=dl, step_type='auto')
            
            L = self.L
            
            if k > k_start:

                n = np.min([len(L)/2, len(L_prev)])

                S = 0
                for i in xrange(n):
                    S+= (L[i] - L_prev[i/2])**2 

                Delta.append(np.sqrt(S / n))

                print k, Delta[k]

                if Delta[k] < DELTA:
                    break
                    
        return np.array(self.H), np.array(self.U), np.array(self.T)
    
    
    def get_results(self, return_list=[]):
        
        TMP = []
        for i in return_list:
            
            if i == 'U': TMP.append(np.array(self.U))
            elif i == 'T': TMP.append(np.array(self.T))    
            elif i == 'C1': TMP.append(np.array(self.C1))
            elif i == 'C2': TMP.append(np.array(self.C2))
            elif i == 'H': TMP.append(np.array(self.H)) 
            elif i == 'L' : TMP.append(np.array(self.L))
        
        return TMP
    
    
       
    #automatic solving with given h_ 
    def step_solve(self, h0=0.1, step_type='auto', H=None, return_list=[]):
        
        self.reset()
        self.h_ = h0
        
        self.T = [self.t]
        self.U = [self.u]
        self.C1 = [np.linalg.norm(self.c1)]
        self.C2 = [np.linalg.norm(self.c2)]
        self.H = [self.h]
        self.L = [self.l]
 

        if step_type == 'given':
            self.h = H[0]
            self.H = H
            for i in xrange(0, H.shape[0]):
                
                self.step(H[i])
                self.T.append(self.t)
                self.U.append(self.u)
                self.C1.append(np.linalg.norm(self.c1))
                self.C2.append(np.linalg.norm(self.c2))          
                self.L.append(self.l)
                          

        else:
            
            self.h = h0
             
            while self.l <= self.L_:

                if step_type == 'auto':
                    self.step()
                    self.H.append(self.h)

                elif step_type == 'const':
                    self.step(h0)
                    self.H.append(h0)

                self.T.append(self.t)
                self.U.append(self.u)
                self.C1.append(np.linalg.norm(self.c1))
                self.C2.append(np.linalg.norm(self.c2))          
                self.L.append(self.l)

        return self.get_results(return_list)
    
    def getU(self, i=None):

        if i == None:
            return np.array(self.U)

        else:
            return np.array(self.U)[:,i]

    def getT(self):
        return np.array(self.T)
    
    def getH(self):
        return np.array(self.H)
    
    def getL(self):
        return np.array(self.L)
        
    def getCappa(self, i=1):
        if i==1:
            return np.array(self.C1)
        else:
            return np.array(self.C2)
     
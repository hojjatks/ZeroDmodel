
import numpy as np

class ZeroDFullydynamic:

    def __init__(self, tf, nsteps, u_init,
                 a=0.01, b=0.02, sigma=0.1, drs=0.1, m=1,k=1, c=1,vpl=1,mu_star=0.6,vstar=1e-6):
        self.a = a
        self.b = b
        self.sigma = sigma
        self.drs = drs
        self.m = m
        self.c = c
        self.tf=tf
        self.k=k
        self.nsteps=nsteps
        self.u_init = u_init
        self.vpl=vpl
        self.mu_star=mu_star
        self.vstar=vstar

    def RHS(self, u, t):
        x1, x2, x3 = u
        mu=self.mu_star+self.a*np.log((x2+self.vpl)/self.vstar) + self.b * np.log((self.vstar * x3 )/self.drs)
        f1 = x2

        f2 = -(1/self.m)*(self.c *(x2 +self.vpl) + self.k*x1+self.sigma*mu)
        f3 = 1-(x2+self.vpl)*x3/self.drs
        print(mu)
        f = [f1, f2, f3]
        return f

    def solve(self):
        time = np.linspace(0, self.tf, self.nsteps+1)
        solver = ODESolver(self.RHS)
        solver.set_ics(self.u_init)
        u, t = solver.solve(time)
        return u, t
    
class ODESolver:

    def __init__(self, f):
        self.f = lambda u, t: np.asarray(f(u, t), float)

    def set_ics(self, U0):
        U0 = np.asarray(U0)
        self.neq = U0.size
        self.U0 = U0

    def advance(self):
        u, f, k, t = self.u, self.f, self.k, self.t
        dt = t[k+1] - t[k]
        K1 = dt*f(u[k], t[k])
        K2 = dt*f(u[k] + 0.5*K1, t[k] + 0.5*dt)
        K3 = dt*f(u[k] + 0.5*K2, t[k] + 0.5*dt)
        K4 = dt*f(u[k] + K3, t[k] + dt)
        u_new = u[k] + (1/6.0)*(K1 + 2*K2 + 2*K3 + K4)
        return u_new

    def solve(self, time):
        self.t = np.asarray(time)
        n = self.t.size
        self.u = np.zeros((n,self.neq))
        self.u[0] = self.U0
        for k in range(n-1):
            self.k = k
            self.u[k+1] = self.advance()
        return self.u[:k+2], self.t[:k+2]
    
class ZeroDQuasidynamic:

    def __init__(self, Tmax, u_init,
                 a, b, sigma, drs, G,cs,vpl):
        self.a = a
        self.b = b
        self.sigma = sigma
        self.drs = drs
        self.c = (G/(2*cs))
        self.k=  G
        self.vpl=vpl
        self.Tmax=Tmax
        self.u_init = u_init



    def RHS(self, u):
        x1, x2 = u
        f2 = 1-(x1*x2/self.drs)
        f1 = ((self.c+self.a*self.sigma/x1)**-1)*(-self.k*(x1-self.vpl)-self.sigma*self.b*f2/x2)
        f = np.array([f1, f2])

        return f

    def solve(self):
        time = np.linspace(0, self.tf, self.nsteps+1)
        solver = ODESolver(self.RHS)
        solver.set_ics(self.u_init)
        u, t = solver.solve(time)
        return u, t
    

def advance(dt,u,f):
    K1 = dt*f(u)
    K2 = dt*f(u + 0.5*K1)
    K3 = dt*f(u + 0.5*K2)
    K4 = dt*f(u + K3)
    u_new = u + (1/6.0)*(K1 + 2*K2 + 2*K3 + K4)
    # u_new = u + K1
    return u_new.reshape(1,2)

def solve(f,Tmax,u0,dt_const):
    t=0
    u=u0
    i=1
    while t<Tmax and i<2:
        i+=1
        dt=dt_const/u[-1,0]
        u_new=advance(dt,u[-1,:],f)
        u=np.append(u,u_new,axis=0)









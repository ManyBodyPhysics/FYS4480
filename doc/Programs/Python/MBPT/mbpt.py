from numpy import *

class electronbasis():
    def __init__(self, N, rs, Nparticles):
        ############################################################
        ##
        ##  Initialize basis: 
        ##  N = number of shells
        ##  rs = parameter for volume 
        ##  Nparticles = Number of holes (conflicting naming, sorry)
        ##
        ###########################################################
        
        self.rs = rs
        self.states = []
        self.nstates = 0
        self.nparticles = Nparticles
        self.nshells = N - 1
        self.Nm = N + 1
        
        self.k_step = 2*(self.Nm + 1)
        Nm = N
        n = 0 #current shell
        ene_integer = 0
        while n <= self.nshells:
            is_shell = False
            for x in range(-Nm, Nm+1):
                for y in range(-Nm, Nm+1):
                    for z in range(-Nm,Nm+1):
                        e = x*x + y*y + z*z
                        if e   == ene_integer:
                            is_shell = True
                            self.nstates += 2
                            self.states.append([e, x,y,z,1])
                            self.states.append([e, x,y,z, -1])
                            
            if is_shell:
                n += 1
            ene_integer += 1
        self.L3 = (4*pi*self.nparticles*self.rs**3)/3.0
        self.L2 = self.L3**(2/3.0)
        self.L = pow(self.L3, 1/3.0)
        
        for i in range(self.nstates):
            self.states[i][0] *= 2*(pi**2)/self.L**2 #Multiplying in the missing factors in the single particle energy
        self.states = array(self.states) #converting to array to utilize vectorized calculations    
        
    def hfenergy(self, nParticles):
        #Calculate the HF-energy (reference energy) for nParticles particles
        e0 = 0.0
        if nParticles<=self.nstates:
            for i in range(nParticles):
                e0 += self.h(i,i)
                for j in range(nParticles):
                    if j != i:
                        e0 += .5*self.v(i,j,i,j)
        else:
            #Safety for cases where nParticles exceeds size of basis
            print("Not enough basis states.")
            
        return e0
                
    def h(self, p,q):
        #Return single particle energy
        return self.states[p,0]*(p==q)

    
    def v(self,p,q,r,s):
        #Two body interaction for electron gas
        val = 0
        terms = 0.0
        term1 = 0.0
        term2 = 0.0
        kdpl = self.kdplus(p,q,r,s)
        if kdpl != 0:
            val = 1.0/self.L3
            if self.kdspin(p,r)*self.kdspin(q,s)==1:
                if self.kdwave(p,r) != 1.0:
                    term1 = self.L2/(pi*self.absdiff2(r,p))
            if self.kdspin(p,s)*self.kdspin(q,r)==1:
                if self.kdwave(p,s) != 1.0:
                    term2 = self.L2/(pi*self.absdiff2(s,p))
        return val*(term1-term2)

    
    #The following is a series of kroenecker deltas used in the two-body interactions. 
    #Just ignore these lines unless you suspect an error here
    def kdi(self,a,b):
        #Kroenecker delta integer
        return 1.0*(a==b)
    def kda(self,a,b):
        #Kroenecker delta array
        d = 1.0
        for i in range(len(a)):
            d*=(a[i]==b[i])
        return d
    def kdfullplus(self,p,q,r,s):
        #Kroenecker delta wavenumber p+q,r+s
        return self.kda(self.states[p][1:5]+self.states[q][1:5],self.states[r][1:5]+self.states[s][1:5])
    def kdplus(self,p,q,r,s):
        #Kroenecker delta wavenumber p+q,r+s
        return self.kda(self.states[p][1:4]+self.states[q][1:4],self.states[r][1:4]+self.states[s][1:4])
    def kdspin(self,p,q):
        #Kroenecker delta spin
        return self.kdi(self.states[p][4], self.states[q][4])
    def kdwave(self,p,q):
        #Kroenecker delta wavenumber
        return self.kda(self.states[p][1:4],self.states[q][1:4])
    def absdiff2(self,p,q):
        val = 0.0
        for i in range(1,4):
            val += (self.states[p][i]-self.states[q][i])*(self.states[p][i]-self.states[q][i])
        return val

        
def MBPT2(bs):
    #2. order MBPT Energy 
    Nh = bs.nparticles
    Np = bs.nstates-bs.nparticles #Note the conflicting notation here. bs.nparticles is number of hole states 
    vhhpp = zeros((Nh**2, Np**2))
    vpphh = zeros((Np**2, Nh**2))
    #manual MBPT(2) energy (Should be -0.525588309385 for 66 states, shells = 5, in this code)
    psum2 = 0
    for i in range(Nh):
        for j in range(Nh):
            for a in range(Np):
                for b in range(Np):
                    #val1 = bs.v(i,j,a+Nh,b+Nh)
                    #val2 = bs.v(a+Nh,b+Nh,i,j)
                    vhhpp[i + j*Nh, a+b*Np] = bs.v(i,j,a+Nh,b+Nh)
                    vpphh[a+b*Np,i + j*Nh] = bs.v(a+Nh,b+Nh,i,j)/(bs.states[i,0] + bs.states[j,0] - bs.states[a + Nh, 0] - bs.states[b+Nh,0])
    psum = .25*sum(dot(vhhpp,vpphh).diagonal())
    return psum
    
def MBPT2_fast(bs):
    #2. order MBPT Energy 
    Nh = bs.nparticles
    Np = bs.nstates-bs.nparticles #Note the conflicting notation here. bs.nparticles is number of hole states 
    vhhpp = zeros((Nh**2, Np**2))
    vpphh = zeros((Np**2, Nh**2))
    #manual MBPT(2) energy (Should be -0.525588309385 for 66 states, shells = 5, in this code)
    psum2 = 0
    for i in range(Nh):
        for j in range(i):
            for a in range(Np):
                for b in range(a):
                    val = bs.v(i,j,a+Nh,b+Nh)
                    eps = val/(bs.states[i,0] + bs.states[j,0] - bs.states[a + Nh, 0] - bs.states[b+Nh,0])
                    vhhpp[i + j*Nh, a+b*Np] = val 
                    vhhpp[j + i*Nh, a+b*Np] = -val 
                    vhhpp[i + j*Nh, b+a*Np] = -val
                    vhhpp[j + i*Nh, b+a*Np] = val 
        
                    
                    vpphh[a+b*Np,i + j*Nh] = eps
                    vpphh[a+b*Np,j + i*Nh] = -eps
                    vpphh[b+a*Np,i + j*Nh] = -eps
                    vpphh[b+a*Np,j + i*Nh] = eps
                    
                    
    psum = .25*sum(dot(vhhpp,vpphh).diagonal())
    return psum


#user input here
number_of_shells = 5
number_of_holes = 14 #(particles)


#initialize basis    
bs = electronbasis(number_of_shells,1.0,number_of_holes) #shells, r_s = 1.0, holes

#Print some info to screen
print ("Number of shells:", number_of_shells)
print ("Number of states:", bs.nstates)
print ("Number of holes :", bs.nparticles)
print ("Reference Energy:", bs.hfenergy(number_of_holes), "hartrees ")
print ("                :", 2*bs.hfenergy(number_of_holes), "rydbergs ")

print ("Ref.E. per hole :", bs.hfenergy(number_of_holes)/number_of_holes, "hartrees ")
print ("                :", 2*bs.hfenergy(number_of_holes)/number_of_holes, "rydbergs ")



#calculate MBPT2 energy
print ("MBPT2 energy    :", MBPT2_fast(bs), " hartrees")

## Coupled clusters in CCD approximation
## Implemented for the pairing model of Lecture Notes in Physics 936, Chapter 8.
## Thomas Papenbrock, June 2018

import numpy as np


def init_pairing_v(g,pnum,hnum):
    """
    returns potential matrices of the pairing model in three relevant channels
    
    param g: strength of the pairing interaction, as in Eq. (8.42)
    param pnum: number of particle states
    param hnum: number of hole states
    
    return v_pppp, v_pphh, v_hhhh: np.array(pnum,pnum,pnum,pnum), 
                                   np.array(pnum,pnum,hnum,hnum), 
                                   np.array(hnum,hnum,hnum,hnum), 
                                   The interaction as a 4-indexed tensor in three channels.
    """
    v_pppp=np.zeros((pnum,pnum,pnum,pnum))
    v_pphh=np.zeros((pnum,pnum,hnum,hnum))
    v_hhhh=np.zeros((hnum,hnum,hnum,hnum))
    
    gval=-0.5*g
    for a in range(0,pnum,2):
        for b in range(0,pnum,2):
            v_pppp[a,a+1,b,b+1]=gval
            v_pppp[a+1,a,b,b+1]=-gval
            v_pppp[a,a+1,b+1,b]=-gval
            v_pppp[a+1,a,b+1,b]=gval
            
    for a in range(0,pnum,2):
        for i in range(0,hnum,2):
            v_pphh[a,a+1,i,i+1]=gval
            v_pphh[a+1,a,i,i+1]=-gval
            v_pphh[a,a+1,i+1,i]=-gval
            v_pphh[a+1,a,i+1,i]=gval
    
    for j in range(0,hnum,2):
        for i in range(0,hnum,2):
            v_hhhh[j,j+1,i,i+1]=gval
            v_hhhh[j+1,j,i,i+1]=-gval
            v_hhhh[j,j+1,i+1,i]=-gval
            v_hhhh[j+1,j,i+1,i]=gval
        
    return v_pppp, v_pphh, v_hhhh
    
    
def init_pairing_fock(delta,g,pnum,hnum):
    """
    initializes the Fock matrix of the pairing model
    
    param delta: Single-particle spacing, as in Eq. (8.41)
    param g: pairing strength, as in eq. (8.42)
    param pnum: number of particle states
    param hnum: number of hole states
    
    return f_pp, f_hh: The Fock matrix in two channels as numpy arrays np.array(pnum,pnum), np.array(hnum,hnum). 
    """
# the Fock matrix for the pairing model. No f_ph needed, because we are in Hartree-Fock basis 
    deltaval=0.5*delta
    gval=-0.5*g
    f_pp = np.zeros((pnum,pnum))
    f_hh = np.zeros((hnum,hnum))

    for i in range(0,hnum,2):
        f_hh[i  ,i  ] = deltaval*i+gval
        f_hh[i+1,i+1] = deltaval*i+gval
        
    for a in range(0,pnum,2):
        f_pp[a  ,a  ] = deltaval*(hnum+a)
        f_pp[a+1,a+1] = deltaval*(hnum+a)
    
    return f_pp, f_hh


def init_t2(v_pphh,f_pp,f_hh):
    """
    Initializes t2 amlitudes as in MBPT2, see first equation on page 345
    
    param v_pphh: pairing tensor in pphh channel
    param f_pp:   Fock matrix in pp channel
    param f_hh:   Fock matrix in hh channel
    
    return t2: numpy array in pphh format, 4-indices tensor
    """
    pnum = len(f_pp)
    hnum = len(f_hh)
    t2_new = np.zeros((pnum,pnum,hnum,hnum))
    for i in range(hnum):
        for j in range(hnum):
            for a in range(pnum):
                for b in range(pnum):
                    t2_new[a,b,i,j] = v_pphh[a,b,i,j] / (f_hh[i,i]+f_hh[j,j]-f_pp[a,a]-f_pp[b,b])
    return t2_new


# CCD equations. Note that the "->abij" assignment is redundant, because indices are ordered alphabetically.
# Nevertheless, we retain it for transparency.
def ccd_iter(v_pppp,v_pphh,v_hhhh,f_pp,f_hh,t2):
    """
    Performs one iteration of the CCD equations (8.34), using also intermediates for the nonliniar terms
    
    param v_pppp: pppp-channel pairing tensor, numpy array
    param v_pphh: pphh-channel pairing tensor, numpy array
    param v_hhhh: hhhh-channel pairing tensor, numpy array
    param f_pp: Fock matrix in pp channel
    param f_hh: Fock matrix in hh channel
    param t2: Initial t2 amplitude, tensor in form of pphh channel
    
    return t2_new: new t2 amplitude, tensor in form of pphh channel
    """
    pnum = len(f_pp)
    hnum = len(f_hh)
    Hbar_pphh = (  v_pphh 
                 + np.einsum('bc,acij->abij',f_pp,t2) 
                 - np.einsum('ac,bcij->abij',f_pp,t2) 
                 - np.einsum('abik,kj->abij',t2,f_hh)
                 + np.einsum('abjk,ki->abij',t2,f_hh)
                 + 0.5*np.einsum('abcd,cdij->abij',v_pppp,t2) 
                 + 0.5*np.einsum('abkl,klij->abij',t2,v_hhhh)
                )

    # hh intermediate, see (8.47)
    chi_hh = 0.5* np.einsum('cdkl,cdjl->kj',v_pphh,t2)

    Hbar_pphh = Hbar_pphh - (  np.einsum('abik,kj->abij',t2,chi_hh) 
                             - np.einsum('abik,kj->abji',t2,chi_hh) )

    # pp intermediate, see (8.46)
    chi_pp = -0.5* np.einsum('cdkl,bdkl->cb',v_pphh,t2)

    Hbar_pphh = Hbar_pphh + (  np.einsum('acij,cb->abij',t2,chi_pp) 
                             - np.einsum('acij,cb->baij',t2,chi_pp) )

    # hhhh intermediate, see (8.48)
    chi_hhhh = 0.5 * np.einsum('cdkl,cdij->klij',v_pphh,t2)

    Hbar_pphh = Hbar_pphh + 0.5 * np.einsum('abkl,klij->abij',t2,chi_hhhh)

    # phph intermediate, see (8.49)
    chi_phph= + 0.5 * np.einsum('cdkl,dblj->bkcj',v_pphh,t2)


    Hbar_pphh = Hbar_pphh + (  np.einsum('bkcj,acik->abij',chi_phph,t2)
                             - np.einsum('bkcj,acik->baij',chi_phph,t2)
                             - np.einsum('bkcj,acik->abji',chi_phph,t2)
                             + np.einsum('bkcj,acik->baji',chi_phph,t2) )
                 
    t2_new=np.zeros((pnum,pnum,hnum,hnum))
    for i in range(hnum):
        for j in range(hnum):
            for a in range(pnum):
                for b in range(pnum):
                    t2_new[a,b,i,j] = (  t2[a,b,i,j] 
                                       + Hbar_pphh[a,b,i,j] / (f_hh[i,i]+f_hh[j,j]-f_pp[a,a]-f_pp[b,b]) )

    return t2_new


def ccd_energy(v_pphh,t2):
    """
    Computes CCD energy. Call as 
    energy = ccd_energy(v_pphh,t2)
    
    param v_pphh: pphh-channel pairing tensor, numpy array
    param t2: t2 amplitude, tensor in form of pphh channel
    
    return energy: CCD correlation energy
    """
    erg = 0.25*np.einsum('abij,abij',v_pphh,t2)
    return erg

###############################
######## Main Program

# set parameters as for model
pnum = 20 # number of particle states
hnum = 10 # number of hole states
delta = 1.0

g = 0.5

print("parameters")
print("delta =", delta, ", g =", g)


# Initialize pairing matrix elements and Fock matrix
v_pppp, v_pphh, v_hhhh = init_pairing_v(g,pnum,hnum)
f_pp, f_hh = init_pairing_fock(delta,g,pnum,hnum)

# Initialize T2 amplitudes from MBPT2
t2 = init_t2(v_pphh,f_pp,f_hh)
erg = ccd_energy(v_pphh,t2)

# Exact MBPT2 for comparison, see last equation on page 365 
exact_mbpt2 = -0.25*g**2*(1.0/(2.0+g) + 2.0/(4.0+g) + 1.0/(6.0+g))
print("MBPT2 energy =", erg, ", compared to exact:", exact_mbpt2)
    
    
# iterate CCD equations niter times
niter=200
mix=0.3
erg_old=0.0
eps=1.e-8
for iter in range(niter):
    t2_new = ccd_iter(v_pppp,v_pphh,v_hhhh,f_pp,f_hh,t2)
    erg = ccd_energy(v_pphh,t2_new)
    myeps = abs(erg-erg_old)/abs(erg)
    if myeps < eps: break
    erg_old=erg
    print("iter=", iter, "erg=", erg, "myeps=", myeps)
    t2 = mix*t2_new + (1.0-mix)*t2 
    
print("Energy = ", erg)
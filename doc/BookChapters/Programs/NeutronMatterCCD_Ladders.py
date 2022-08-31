import numpy as np

##############################################################
# CCD Program for neutron matter with the Minnesota potential.
#
# Thomas Papenbrock, July/August 2018
#
# License: Free Software following Version 3 of GNU General Public License, 
# see https://www.gnu.org/licenses/gpl.html
#######


##########################
# Class for neutron matter 
#######

class MomSpaceBasis:
    """
    momentum-space basis class
    The constructor has the form
    
    MomSpaceBasis(Nmax,kmax)
    
    param Nmax: Number of lattice points in positive kx-direction
    param kmax: Highest lattice momentum (in 1/fm)
    
    return: MomSpaceBasis as a single-partcle basis. 
    attributes of MomSpaceBasis are
    
    dk : lattice spacing in 1/fm
    Lbox : linear dimension (in fm) of cubic box
    nvec : lattice vectors (integers)
    kvec : lattice momentum vectors (floats, in 1/fm)
    ngrid : total number of lattice points 
    """
    def __init__(self,Nmax,kmax,ordered=True):
        """
        the constructor
        
        Generates a cubic lattice in momentum space
        param Nmax: Number of lattice points in positive kx-direction
        param kmax: Highest lattice momentum (in 1/fm)
        param ordered: Optional parameter, True by default, will order lattice points by kinetic energy
        
        return MomSpaceBasis
        """
        self.Nmax = Nmax
        self.dim = 0
        self.ngrid = 0
        self._kvec=[]
        self._nvec=[]

        dk = kmax / Nmax
        self.dk = dk
        self.Lbox = 2.0*np.pi/dk
        
        nx=[]
        nvec=[]
        for i in range(-Nmax,Nmax+1):
            self.dim=self.dim+1
            nx.append(i)
        
        #print('nx=',nx)
        
        for i in nx:
            for j in nx:
                for k in nx:
                    nvec.append(np.array([i,j,k], dtype=int))
                    
        #print('nvec=',nvec)
        self.ngrid=len(nvec)
        
        if ordered:
            #print("ordered")
            norm=np.zeros(self.ngrid,dtype=int)
            for i, vec in enumerate(nvec):
                npvec=np.array(vec,dtype=int)
                norm[i]=np.dot(npvec,npvec)
               # print(i, vec, norm[i])
        
            index=np.argsort(norm)
            #print(index)
            self._nvec=[]
            for i, ind in enumerate(index):
                #print(i, ind, nvec[ind])
                self._nvec.append(nvec[ind])
                
        else: 
            self._nvec=nvec  # a list
            
        self._kvec = np.array(self._nvec)*dk # a numpy array

    
    def kvec(self,indx=-1):
        """
        MomSpaceBasis.kvec(i) returns ith momentum vector
        MomSpaceBasis.kvec() returns all momentum vectors
        
        param indx: index of k-vector to be returned, optional
        return 3-vector (if index non-negative), or all vectors if no index specified
        """
        if indx == -1:
            return self._kvec
        else:
            return self._kvec[indx]
        
    def nvec(self,indx=-1):
        """
        MomSpaceBasis.nvec(i) returns ith lattice vector
        MomSpaceBasis.nvec() returns all lattice vectors
        
        param indx: index of lattice vector to be returned, optional
        return 3-vector (if index non-negative), or all lattice vectors if no index specified
        """
        if indx == -1:
            return self._nvec
        else:
            return self._nvec[indx]
        
    def dens(self,num):
        """
        returns density of system if num particles are present
        param num: int, number of particles
        return dens: float
        """
        return num/(self.Lbox)**3
    
    def update(self,dk):
        """
        Uses dk as new lattice spacing and rescales existing lattice
        param dk: in 1/fm lattice spacing in momentum space
        """
        self.Lbox=2.0*np.pi/dk
        self._kvec = np.array(self._nvec)*dk
    
    def __len__(self):
        """
        overloading of the 'len' function
        """
        return self.ngrid
    
    
############
# useful functions

def magic_numbers(basis):
    """
    param basis: MomSpaceBasis object
    return magic: array of magic numbers
    """
    nvecs = basis.nvec()
    vec=np.array(nvecs[0],dtype=int)
    norm = np.dot(vec,vec)
    magic=[]
    for i in range(1,len(nvecs)):
        vec=np.array(nvecs[i],dtype=int)
        norm2 = np.dot(vec,vec)
        if norm2 > norm: 
            magic.append(2*i)
            norm=norm2
    return magic


def get_dk(rho,Num):
    """
    param rho: desired density
    param Num: magic number of particles
    return dk: grid spacing in momentum space (in 1/fm)
    """
    Lbox = (Num/rho)**(1.0/3.0)
    dk = 2.0*np.pi/Lbox
    return dk

def spbasis_from_MomSpaceBasis(lattice_vecs,st_degen):
    """
    converts a lattice to a single particle basis for spin-isospin degeneracy st_degen
    param lattice_vecs: list of lattice vectors for 1st particle
    param st_degen: spin-isospin degeneracy
    return: basis as a list of momenta
    """
    if st_degen != 2: # for now only neutron matter
        print("Unexpected parameter st_degen")
        return lattice_vecs
    
    basis=[]
    for vec in lattice_vecs:
        for st in range(st_degen):
            basis.append(np.array(vec,dtype=int))
            
    return basis



#########################################################
# Functions for comparisons with infinite free Fermi gas
                         
def kF_from_density(rho,st_degen=2):
    """
    Computes Fermi momentum for given density and spin/isospin degeneracy.
    
    param rho: density in inverse fm cubed
    param st_degen: spin-isospin degeneracy; default 2
    return: Fermi momentum in inverse fm
    """
    res = (6.0*(np.pi)**2*rho/st_degen)**(1.0/3.0)
    return res

def EnergyDensity_FermiGas(kF,st_degen=2):
    """
    Computes energy density of free Fermi gas at Fermi momentum and spin/isospin degeneracy
    param kF: Fermi momentum in inverse fm
    param st_degen: spin-isospin degeneracy; default 2
    return: Energy density in MeV/fm**3
    """
    pvec = np.array([kF,0.0,0.0])
    erg = (st_degen*kF**3/(10.0*np.pi**2)) * Tkin(pvec)
    return erg


########################################################################################
# Functions for CCD of neutron matter
# Implementation uses only pp and hh ladders
# 
########################################################################################


from numba import jit  
# compile a few functions to gain speed; should probably done in Fortran or C++, 
# and called from Python

@jit(nopython=True)
def minnesota_nn(p_out,s1_out,s2_out,p_in,s1_in,s2_in,Lbox):
    """
    The Minnesota potential between two neutrons, not yet anti-symmetrized 
    param p_out: relative out momentum
    param p_in : relative in momentum
    param s1_out, s2_out: spin projections of out particles 1 and 2
    param s1_in, s2_in  : spin projections of in particles 1 and 2
    Lbox : size of momentum box
    return: value of potential in MeV; not anti-symmetrized!
    """
    # parameters. VT is not active between two neutrons (no triplet)
    VR = 200.0
    VS = -91.85  # sign typo in Lecture Notes Physics 936, Chap. 8
    kappaR = 1.487
    kappaS = 0.465
    
    qvec=p_out-p_in
    q2=np.dot(qvec,qvec)
    
    s1_i =spin2spinor(s1_in)
    s2_i =spin2spinor(s2_in)
    s1_o =spin2spinor(s1_out)
    s2_o =spin2spinor(s2_out)
    
    spin_part = 0.5 * ( np.dot(s1_i,s1_o)*np.dot(s2_i,s2_o)
                       -np.dot(s1_i,s2_o)*np.dot(s2_i,s1_o) )
    
    
    pot = spin_part * (  VR*np.exp(-0.25*q2/kappaR) / (Lbox*np.sqrt(kappaR))**3 
                       + VS*np.exp(-0.25*q2/kappaS) / (Lbox*np.sqrt(kappaS))**3 )
    
    pot = pot*(np.sqrt(np.pi))**3 

    
    return pot

@jit
def spin_of_index(i):
    """
    Even indices of the lattive have spin up, odds have spin down
    param i: index of sp_basis
    return: spin as +/- 1
    """
    spin = 1-2*np.remainder(i,2)
    return spin

@jit
def spin2spinor(s):
    """
    Makes a two-component spinor of an integer s
    param s: spin = +/- 1
    return: two-component numpy array [1,0] for up and [0,1] for down
    """
    up  =np.array([1.0,0.0])
    down=np.array([0.0,1.0])
    if s == 1:
        return up
    else:
        return down

@jit
def Tkin(pvec):
    """
    Kinetic energy for a momentum vector
    param pvec: 3-component numpy array in inverse fm
    return: kinetic energy of that momentum in MeV
    """
    nucleon_mass = 938.92
    hbarc = 197.33
# More precise numbers for neutron mass and hbar follow.
# For N=14, this yields E_HF = 10.3337 MeV per nucleon in HF. Benchmarked with Ragnar Stroberg.
#     nucleon_mass = 939.56563
#     hbarc = 197.3269718
    p2 = np.dot(pvec,pvec)
    res = 0.5*hbarc**2*p2/nucleon_mass
    return res
    
@jit
def compute_total_Tkin(Nocc,sp_basis,dk):
    """
    Computes total kinetic energy of reference state
    param Nocc, sp_basis, dk: particle number, integer s.p. lattice, delta k 
    return: total kinetic energy
    """
    erg=0.0
    for i in range(Nocc):
        mom_vec = sp_basis[i]
        vec=np.array(mom_vec)*dk
        erg=erg+Tkin(vec)
        
    return erg



@jit
def Fock(pvec,s,sp_basis,Nocc,dk,Lbox):
    """
    Fock matrix of momentum pvec in hh space
    param pvec: 3-component numpy array in inverse fm
    param s: spin as +/- 1 of state
    param_sp_basis, Nocc, dk, Lbox : parameters of s.p. basis and system
    
    return: Fock matrix element = kinetic energy of that momentum in MeV
    """
    res = Tkin(pvec)
    
    dum=0.0
    for i in range(Nocc):
        vec=sp_basis[i]*dk
        si=spin_of_index(i)
        p_in = 0.5*(vec-pvec)
        p_out= p_in
        dum = dum + ( minnesota_nn(p_out,s,si, p_in, s,si,Lbox) 
                     -minnesota_nn(p_out,s,si,-p_in,si, s,Lbox) )  #antisymmetrized Minnesota
        
    res = res+dum
    return res

def compute_E_HF_simple(Nocc,sp_basis,dk):
    """
    Computes HF energy of reference state
    param Nocc, sp_basis, dk: particle number, integer s.p. lattice, delta k 
    return: total HF energy
    """
    erg=compute_total_Tkin(Nocc,sp_basis,dk)

    pot=0.0
    for i in range(Nocc):
        momi=sp_basis[i]*dk
        si = spin_of_index(i)
        for j in range(Nocc):
            momj=sp_basis[j]*dk
            sj = spin_of_index(j)
            p_rel = 0.5*(momi-momj)
            pot = pot + 0.5* (  minnesota_nn(p_rel,si,sj, p_rel,si,sj,Lbox)
                                - minnesota_nn(p_rel,si,sj,-p_rel,sj,si,Lbox) )
            
    erg = erg+pot
    return erg


def get_channels(sp_basis,start1,end1,start2,end2,identical,other_channels=None):
    """
    Returns channels for coupled cluster based on Minnesota potential
    param sp_Basis: A single-particle basis
    param start1: index to start for particle 1
    param end1: index to end for particle 1
    param start2: index to start for particle 2
    param end2: index to end for particle 2
    param identical: True for hh or pp, False for hp
    param other_channels: list of other channels to compare with
    return: channels, p_rel, t2amp. channels is a list of p12, where p12 is a momentum vector; 
            p_rel is a nested list with relative momenta and spins for each channel
    """
    channel=[]
    p_rel=[]
    for i, mom_vecs1 in enumerate(sp_basis[start1:end1]):
        #vec1=np.array(mom_vecs1,dtype=int)
        vec1=mom_vecs1
        spin1=spin_of_index(i)
            
        for j, mom_vecs2 in enumerate(sp_basis[start2:end2]):
            if identical and i==j: continue  #Fortran cycle
            #vec2=np.array(mom_vecs2,dtype=int)
            vec2=mom_vecs2
            spin2=spin_of_index(j)
                
            p12 = vec1+vec2
            prel= vec1-vec2
            spins=np.array([spin1,spin2],dtype=int)
            ps=[prel,spins]

            new=True
            needed=True
            if other_channels is not None: #check whether we need this channel
                needed=False
                for chan_o in other_channels:
                    if (chan_o==p12).all(): 
                        needed=True
                        break
            if needed: #check whether this channel exists already
                for ipos, chan in enumerate(channel):
                    if (chan==p12).all(): 
                        new=False
                        break
                    
            if needed and new: 
                channel.append(p12)
                p_rel.append([ps])
                
            if needed and not new:
                p_rel[ipos].append(ps)
                
    return channel, p_rel 
    

    
def setup_T2_amplitudes(sp_basis,NN,st_degen):
    """
    returns the t2 amplitudes and t2 channels
    param sp_basis: a sp_basis
    param NN: neutron number
    param st_degen: 2 for the moment, spin-isospin degeneracy
    return: hh_channels, pp_channels, p_rel_hh, p_rel_pp, t2amp
            these are the hh and pp channels of T2, lists of the relative momenta, 
            and t2amps as a list of numpy arrays set to zero 
    """
    num_states = len(sp_basis)
    
    hh_channels, p_rel_hh = get_channels(sp_basis,0,NN,0,NN,True)
    print('hh channels=', len(hh_channels))

    pp_channels, p_rel_pp = get_channels(sp_basis,NN,num_states,NN,num_states,True,hh_channels)
    print('pp channels=', len(pp_channels))
    
    if len(pp_channels) != len(hh_channels): print('pp and hh channels do not match')
     

    ordered_pp_channel=[]
    ordered_p_rel_pp=[]
    for i, chanhh in enumerate(hh_channels):
        for j, chanpp in enumerate(pp_channels):
            if (chanpp==chanhh).all():
                ordered_pp_channel.append(chanpp)
                ordered_p_rel_pp.append(p_rel_pp[j])    
                break
                
    pp_channels = ordered_pp_channel
    p_rel_pp = ordered_p_rel_pp
            
    # set t2 amplitudes to zero in each channel
    t2amp = fill_pot(Lbox, dk, pp_channels, hh_channels, p_rel_pp, p_rel_hh, True)
        
    return hh_channels, pp_channels, p_rel_hh, p_rel_pp, t2amp

def fill_pot(Lbox, dk, channels_out, channels_in, p_rel_out, p_rel_in, T2amp=False):
    """
    Fills lists of matrices such as Vhhhh, Vpphh, Vpppp, t2_pphh
    param Lbox: Lbox
    param dk: dk
    param channels_out, channels_in: the channels we have
    param p_rel_out, p_rel_in: the list of [prel, [s1,s2]]
    param T2amp=False: Set to True if t2_pphh needs to be computed
    return: The object of desire as a list of numpy matrices. 
            Contain matrix elements for potentials, zeros if T2amp=True is requested. 
    """
    Vpot=[]
    for i, chan_in in enumerate(channels_in):
        dim_in = len(p_rel_in[i])
        dim_out= len(p_rel_out[i])
        Vpot_chan=np.zeros((dim_out,dim_in))
        if not T2amp: 
            for ii, ps_i in enumerate(p_rel_in[i]):
                [pii, [s1, s2]] = ps_i
                pii = pii*dk*0.5
                for jj, ps_j in enumerate(p_rel_out[i]):
                    if dim_in == dim_out and jj > ii: continue
                    [pjj, [ss1, ss2]] = ps_j
                    pjj = pjj*dk*0.5
                    Vpot_chan[jj,ii] = ( minnesota_nn( pjj,ss1,ss2, pii,s1,s2,Lbox)
                                        -minnesota_nn(-pjj,ss2,ss1, pii,s1,s2,Lbox) )
                    if dim_in == dim_out : Vpot_chan[ii,jj] = Vpot_chan[jj,ii]
        
        Vpot.append(Vpot_chan)
    return Vpot


def init_V(Lbox, dk, hhchannels, ppchannels, p_relhh, p_relpp,zeros=False):
    """
    Sets up Vhhhh, Vpphh, and Vpppp. 
    
    return: Vhhhh, Vpphh, Vpppp as a lists of numpy arrays
    """
    Vhhhh = fill_pot(Lbox, dk, hhchannels, hhchannels, p_relhh, p_relhh,zeros)
    Vpphh = fill_pot(Lbox, dk, ppchannels, hhchannels, p_relpp, p_relhh,zeros)
    Vpppp = fill_pot(Lbox, dk, ppchannels, ppchannels, p_relpp, p_relpp,zeros)
    
    return Vhhhh, Vpphh, Vpppp  
            
@jit
def make_diagram(obj1,obj2,fac):
    """
    Makes diagrams for pp or hh ladders as matrix-matrix multiplications
    """
    hbar_pphh=[]
    dim1=len(obj1)
    for chan in range(dim1):
        mat1 = obj1[chan]
        mat2 = obj2[chan]
        hbar_pphh.append( fac*np.matmul(mat1,mat2) )
        
    return hbar_pphh


def make_diagrams2_3(t2_pphh,fabij):
    hbar_pphh=[]
    for i, t2_mat in enumerate(t2_pphh):
        f_mat = fabij[i]
        hbar_mat = t2_mat*f_mat
        hbar_pphh.append(hbar_mat)
    return hbar_pphh


def compute_hbar(v_pppp,v_pphh,v_hhhh,t2_pphh,fabij):
    diagram1 = v_pphh.copy()
    diagram23 = make_diagrams2_3(t2_pphh, fabij)
    diagram4 = make_diagram(v_pppp,t2_pphh,0.5)
    diagram5 = make_diagram(t2_pphh,v_hhhh,0.5)
    
    hbar_pphh=[]
    for i in range(len(t2_pphh)):
        mat = (  diagram1[i]
               + diagram23[i]
               + diagram4[i]
               + diagram5[i]  )
        hbar_pphh.append(mat)
    
    return hbar_pphh
    
def get_energy_denominator(hh_channels,p_rel_pp,p_rel_hh,sp_basis,Nocc,dk,Lbox):
    res=[]
    fabij=[]
    for i, Ptot in enumerate(hh_channels):
        dimhh=len(p_rel_hh[i])
        dimpp=len(p_rel_pp[i])
        res_mat = np.zeros((dimpp,dimhh))
        f_mat   = np.zeros((dimpp,dimhh))
        for ii, psh_rel in enumerate(p_rel_hh[i]):
            [pij, [si, sj]] = psh_rel  
            p_i = (Ptot+pij)//2
            p_i = p_i + np.array([Nmax,Nmax,Nmax],dtype=int)
            p_j = (Ptot-pij)//2
            p_j = p_j + np.array([Nmax,Nmax,Nmax],dtype=int)
            ssi = (1-si)//2
            ssj = (1-sj)//2
            fii = fock_mtx4[p_i[0],p_i[1],p_i[2],ssi]
            fjj = fock_mtx4[p_j[0],p_j[1],p_j[2],ssj]
            for jj, psp_rel in enumerate(p_rel_pp[i]):
                [pab, [sa, sb]] = psp_rel
                p_a = (Ptot+pab)//2
                p_a = p_a + np.array([Nmax,Nmax,Nmax],dtype=int)
                p_b = (Ptot-pab)//2
                p_b = p_b + np.array([Nmax,Nmax,Nmax],dtype=int)
                ssa = (1-sa)//2
                ssb = (1-sb)//2
                faa = fock_mtx4[p_a[0],p_a[1],p_a[2],ssa]
                fbb = fock_mtx4[p_b[0],p_b[1],p_b[2],ssb]
                
                res_mat[jj,ii] = 1.0 / (fii + fjj - faa - fbb) 
                f_mat[jj,ii] = faa + fbb - fii - fjj
        res.append(res_mat)
        fabij.append(f_mat)
    return res, fabij


def get_t2_from_mbpt(Vpphh,denom):
    """
    param Vpphh: Vpphh
    param denom: energy denominator in pphh format
    return t2: quotient of both, element for element 
    """
    res = []
    for i, vv in enumerate(Vpphh):
        dd = denom[i]
        res_mat = vv*dd  #how simple in python; element by element multiply
        res.append(res_mat)
    return res


def compute_E_CCD(Vpphh,T2pphh):
    erg=0.0
#     erg2=0.0
    for i, t2mat in enumerate(T2pphh):
        vmat  = Vpphh[i]
        erg = erg + 0.25*np.sum(vmat*t2mat)
    return erg


def compute_Fock_4(sp_basis,Nocc,dk,Lbox,Nmax):
    fock_mtx4=np.zeros(shape=(2*Nmax+1, 2*Nmax+1, 2*Nmax+1, 2))
    for i, vec in enumerate(sp_basis):
        pvec=vec*dk
        spin=spin_of_index(i)
        si = (1 - spin)//2
        px=vec[0]+Nmax
        py=vec[1]+Nmax
        pz=vec[2]+Nmax
        fock_mtx4[px,py,pz,si] = Fock(pvec,spin,sp_basis,Nocc,dk,Lbox)
    return fock_mtx4

#####################################
########### Main Program starts here    

from timeit import default_timer as timer
# for timing purposes

progstart=timer()

Nmax=1
kmax=1.0
mbase = MomSpaceBasis(Nmax,kmax)
lattice=mbase.nvec()

## set particle number
NN=14
st_degen=2  # spin up and down
print("chosen N =", NN)
print("magic numbers", magic_numbers(mbase))

## set density
rho=0.08

dk = get_dk(rho,NN)

mbase.update(dk)
Lbox = mbase.Lbox


## get single particle basis

sp_basis = spbasis_from_MomSpaceBasis(lattice,st_degen)
num_states = len(sp_basis)

print('number of s.p. states:', num_states)

# print out a few facts of the reference state
total_Tkin = compute_total_Tkin(NN,sp_basis,dk)
print('total Tkin per particle =', total_Tkin/NN )

k_fermi = kF_from_density(rho)

print("Fermi momentum =", k_fermi)

E_gas = EnergyDensity_FermiGas(k_fermi)

print("Energy per neutron of infinite free Fermi gas", E_gas/rho)

E_HF = compute_E_HF_simple(NN,sp_basis,dk)
E_HF = E_HF/NN
print("HF energy per neutron =", E_HF)

## now we start our business ...
## get all channels and two-body states within those channels; set T2 to zero
hh_channels, pp_channels, p_rel_hh, p_rel_pp, t2_pphh = setup_T2_amplitudes(sp_basis,NN,st_degen)

# get some insight in how big this all is
count=0
for i, channel in enumerate(p_rel_hh):
    dim=len(p_rel_hh[i])
    count=count+dim

print('hh number of total combinations', count)

count=0
for i, channel in enumerate(p_rel_pp):
    dim=len(p_rel_pp[i])
    count=count+dim
    
print('pp number of total combinations', count)


print("get v_hhhh, v_pphh, v_pppp")
start = timer()
v_hhhh, v_pphh, v_pppp = init_V(Lbox, dk, hh_channels, pp_channels, p_rel_hh, p_rel_pp)
end = timer()
print("what a hog!", end-start, 'seconds')


print("compute energy denominator")
start = timer()
fock_mtx4 = compute_Fock_4(sp_basis,NN,dk,Lbox,Nmax)
denom_pphh, f_abij = get_energy_denominator(pp_channels,p_rel_pp,p_rel_hh,sp_basis,NN,dk,Lbox)
end = timer()
print("that's faster", end-start, 'seconds')

print("Initialize T2 from MBPT2")
t2_pphh = get_t2_from_mbpt(v_pphh,denom_pphh)

erg = compute_E_CCD(v_pphh,t2_pphh)
print('MBPT2 correlation energy per neutron =', erg/NN)


print("start CCD iterations ...")

niter=200
mix=0.99
erg_old=0.0
eps=1.e-8
for iter in range(niter):
    
    start = timer()
    hbar_pphh = compute_hbar(v_pppp,v_pphh,v_hhhh,t2_pphh,f_abij)
    end = timer()
    print("time of making Hbar:", end-start, 'seconds')
    
    t2_new = get_t2_from_mbpt(hbar_pphh,denom_pphh)
    
    for i in range(len(t2_new)):
        t2_new[i] = t2_pphh[i] + t2_new[i]
    
    erg = compute_E_CCD(v_pphh,t2_new)
    
    myeps = abs(erg-erg_old)/abs(erg)
    if myeps < eps: break
    erg_old=erg

    print("iter=", iter, "Correlation energy per neutron=", erg/NN, ", epsilon=", myeps)
    
    for i in range(len(t2_pphh)):
        t2_pphh[i] = mix*t2_new[i] + (1.0-mix)*t2_pphh[i]
    
print("Correlation energy per neutron= ", erg/NN)

progend=timer()
print('total time in seconds', progend-progstart)


# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 15:59:04 2020

@author: Zeng Qiyu
"""
import os  
import numpy as np
import matplotlib.pyplot as plt
import time
import glob

#from matplotlib.ticker import FuncFormatter
#plt.rcParams['font.family'] = ['Times New Roman']

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')

e=1.602e-19                           #charge of a electron, C / electron volt, J
epsilon0 = 8.854e-12                  #vaccum dielectric constant, F/m
me = 9.109e-31                        #mass of a electron , kg
h = 6.626e-34                         #Planck constant, J*s 
hbar = h/(2*np.pi)
NA = 6.022e23
kb = 1.38e-23
aB = 4*np.pi*epsilon0*hbar**2 / (me*e**2)                    #Bohr radius, m
Hatree = 13.6                         #Hatree, eV

# ---------- Metal Units ------------- #
J2eV = 6.2442e18
m2A = 1e10
gmol2g = 1/NA
kg2gmol = 1e3/gmol2g
s2ps = 1e12
s2fs = 1e15
# --------- Atomic Units ------------- #
m2bohr = 1/aB
A2bohr = m2bohr/m2A


# -------- pressure ------------------ #
bar2Pa = 1e5
kbar2GPa = 1e3*bar2Pa/1e9


class Trajectory():
    
    def __init__(self,file_dir,label, SKIP, INTERVAL):
        
        self.dir = file_dir
        self.label = label
        
        if self.label =='raw':
            self.data = np.loadtxt(self.dir+'coord.raw')
            self.N_total = self.data.shape[0] 
            self.N_atom = self.data[0].reshape((-1,3)).shape[0]
        elif self.label =='dump':
            self.file_list = glob.glob(self.dir+'dump.*')
            self.N_total = len(self.file_list)
            self.N_atom = np.loadtxt(self.file_list[0],skiprows=9)[:,2:].shape[0]
            
        temp = np.arange(SKIP,self.N_total, INTERVAL)
        self.N_t = temp.shape[0]
    
    
    def _init_rkt(self, a_label):
        
        self.ave_label = a_label
        if self.ave_label == 'single':
            self.rkt_re = np.zeros((self.N_t,))
            self.rkt_im = np.zeros((self.N_t,))
        if self.ave_label == 'xyz':
            self.rkt_re = np.zeros((self.N_t,3))
            self.rkt_im = np.zeros((self.N_t,3))        
        if self.ave_label == 'all':
            self.rkt_re = np.zeros((self.N_t,6))
            self.rkt_im = np.zeros((self.N_t,6))
    
    def _get_r(self,n_t):
        
        if self.label == 'raw':
            self.position = self.data[n_t].reshape((-1,3))
        elif self.label =='dump':
            self.position = np.loadtxt(self.dir+'dump.'+str(int(n_t)),skiprows=9)[:,2:5]
        
    def _get_rkt(self, vec_k, n_t):
        
        theta_x = vec_k[0]*self.position[:,0]
        theta_y = vec_k[1]*self.position[:,1]
        theta_z = vec_k[2]*self.position[:,2]

        self.rkt_re[n_t] = np.sum(np.cos( theta_x+theta_y+theta_z ))
        self.rkt_im[n_t] = - np.sum(np.sin( theta_x+theta_y+theta_z ))
        
    def _get_xyz_rkt(self, k_value, n_t):
        
        theta_x = k_value*self.position[:,0]
        theta_y = k_value*self.position[:,1]
        theta_z = k_value*self.position[:,2]

        self.rkt_re[n_t,0] = np.sum(np.cos( theta_x ))
        self.rkt_re[n_t,1] = np.sum(np.cos( theta_y ))
        self.rkt_re[n_t,2] = np.sum(np.cos( theta_z ))

        self.rkt_im[n_t,0] = -np.sum(np.sin( theta_x ))
        self.rkt_im[n_t,1] = -np.sum(np.sin( theta_y ))
        self.rkt_im[n_t,2] = -np.sum(np.sin( theta_z ))
        
    def _get_all_rkt(self, k_value, n_t):
        
        theta_x = k_value*self.position[:,0]
        theta_y = k_value*self.position[:,1]
        theta_z = k_value*self.position[:,2]

        self.rkt_re[n_t,0] = np.sum(np.cos( theta_x ))
        self.rkt_re[n_t,1] = np.sum(np.cos( theta_y ))
        self.rkt_re[n_t,2] = np.sum(np.cos( theta_z ))
        self.rkt_re[n_t,3] = np.sum(np.cos( theta_x ))
        self.rkt_re[n_t,4] = np.sum(np.cos( theta_y ))
        self.rkt_re[n_t,5] = np.sum(np.cos( theta_z ))
        
        self.rkt_im[n_t,0] = -np.sum(np.sin( theta_x ))
        self.rkt_im[n_t,1] = -np.sum(np.sin( theta_y ))
        self.rkt_im[n_t,2] = -np.sum(np.sin( theta_z ))
        self.rkt_im[n_t,3] = np.sum(np.sin( theta_x ))
        self.rkt_im[n_t,4] = np.sum(np.sin( theta_y ))
        self.rkt_im[n_t,5] = np.sum(np.sin( theta_z ))        
        
    def _get_skw(self, dt):
        
        freq = np.fft.fftfreq(self.N_t, d=dt)
        
        if self.ave_label == 'single':
            nkt = np.transpose(np.vstack((self.rkt_re,self.rkt_im)))
            n_kt = np.zeros(self.N_t,dtype=np.complex)
            n_kt.real = nkt[:,0]
            n_kt.imag = nkt[:,1]

            n_kw = np.fft.fft(n_kt)
            
            skw = np.abs(n_kw)**2 / (self.N_t*self.N_atom)
            
            #self.skw[int(self.skw.argmax())] = 0
            
        if self.ave_label == 'xyz':

            skw = np.zeros(self.N_t)
            
            for i in range(3):
                nkt = np.transpose(np.vstack((self.rkt_re[:,i],self.rkt_im[:,i])))
                n_kt = np.zeros(self.N_t,dtype=np.complex)
                n_kt.real = nkt[:,0]
                n_kt.imag = nkt[:,1]

                n_kw = np.fft.fft(n_kt)
                
                skw += np.abs(n_kw)**2 / (self.N_t*self.N_atom * 3)
                
        if self.ave_label == 'all':

            skw = np.zeros(self.N_t)
            
            for i in range(6):
                nkt = np.transpose(np.vstack((self.rkt_re[:,i],self.rkt_im[:,i])))
                n_kt = np.zeros(self.N_t,dtype=np.complex)
                n_kt.real = nkt[:,0]
                n_kt.imag = nkt[:,1]

                n_kw = np.fft.fft(n_kt)
                
                skw += np.abs(n_kw)**2 / (self.N_t*self.N_atom * 6)
        
        out = np.transpose(np.vstack((freq,skw)))
        sort = out[np.lexsort(out[:,::-1].T)]
        
        self.dsf = sort
        
if __name__ == '__main__':

    dir1 = '../02.dpmd/dump/'
    
    outdir = os.popen('pwd').readlines()[0][:-1]+'/'
    
    SKIP=0
    INTERVAL=1
    
    work = Trajectory(dir1,'dump', SKIP,INTERVAL)

    # =============== #
    # step 1: specify the length of cubic simulation box (Units: Ang)
    # =============== #
    L = 16*4.05
    
    delta_k = 2*np.pi/L
    delta_t = 1   # fs
    DUMP = 10   

    ave_label='all'
    
    work._init_rkt('all')
    
    for j in range(1):
        init_time = time.time()
        
        NUM = 1
        #NUM = j+1
        
        k_value = delta_k * NUM
    
        vec_t = np.arange(SKIP, work.N_total,INTERVAL)
    
        for i in range(vec_t.shape[0]):
            work._get_r(vec_t[i]*DUMP)
            work._get_all_rkt(k_value, i)
    
        np.savetxt(outdir+'nkt_'+str(NUM)+'k_'+ave_label+'.txt',np.transpose(np.vstack((work.rkt_re,work.rkt_im))))
    
        work._get_skw(dt = delta_t*INTERVAL*DUMP)
    
        np.savetxt(outdir+'dsf_'+str(NUM)+'k_'+ave_label+'.txt',work.dsf)
        
        end_time = time.time()
        print('Time spent for single k value : {:.2f} seconds'.format(end_time - init_time))

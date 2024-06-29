"""
Author: Ziad (https://github.com/ZiadHatab)

Example code for perofming a two-port/one-path SOLT calibration. 

This code is complementary to my blog post on this topic: https://ziadhatab.github.io/posts/two-port-one-path-calibration/
"""

import numpy as np
import matplotlib.pyplot as plt
import skrf as rf
import os

def s2t(S, pseudo=False):
    T = S.copy()
    T[0,0] = -(S[0,0]*S[1,1]-S[0,1]*S[1,0])
    T[0,1] = S[0,0]
    T[1,0] = -S[1,1]
    T[1,1] = 1
    return T if pseudo else T/S[1,0]

def t2s(T, pseudo=False):
    S = T.copy()
    S[0,0] = T[0,1]
    S[0,1] = T[0,0]*T[1,1]-T[0,1]*T[1,0]
    S[1,0] = 1
    S[1,1] = -T[1,0]
    return S if pseudo else S/T[1,1]

def compute_error_terms(measured, ideal):
    '''
    Main function to compute the five error terms of forward direction.

    measured and ideal are in following order:
    SOLT: [short, open, load, thru]
    '''

    short_m, open_m, load_m, thru_m = measured
    short_i, open_i, load_i, thru_i = ideal
    
    # the error terms 
    a11 = []
    a12 = []
    a21 = []
    beta  = []
    alpha = []
    for inx,_ in enumerate(short_m.frequency.f):
        # build system matrix
        H = np.array([[x[1].s[inx,0,0], 1, -x[0].s[inx,0,0]*x[1].s[inx,0,0]] for x in [(short_m, short_i), (open_m, open_i), (load_m, load_i)]])
        # measurement vector
        y = np.array([[x.s[inx,0,0]] for x in [short_m, open_m, load_m]])

        x = np.linalg.inv(H).dot(y).squeeze() # solve for the error terms [a11, a12, a21]
        a11.append(x[0])
        a12.append(x[1])
        a21.append(x[2])

        # solve for the error terms [alpha, beta]
        A = np.array([[x[0], x[1]],
                      [x[2], 1]])
        T = s2t(thru_i.s[inx])
        S11 = thru_m.s[inx,0,0]
        S21 = thru_m.s[inx,1,0]
        xx = np.linalg.inv(T)@np.linalg.inv(A)@np.array([S11/S21, 1/S21])
        alpha.append(xx[0])
        beta.append(xx[1])
        
    return np.array(a11), np.array(a12), np.array(a21), np.array(alpha), np.array(beta)

def apply_cal(netw, netw2, a11, a12, a21, alpha, beta):
    '''
    Apply the calibration to the network netw using the error terms.
    The network is measured in both directions, forwrard and reverse.
    If the network only measured in forward direction, use the same network for both netw and netw2, 
    which enforces the network to be symmetric.
    '''
    S = []
    for inx,_ in enumerate(netw.frequency.f):
        A = np.array([[a11[inx], a12[inx]],
                      [a21[inx], 1]])
        
        S11 = netw.s[inx,0,0]
        S21 = netw.s[inx,1,0]
        S22 = netw2.s[inx,0,0]
        S12 = netw2.s[inx,1,0]

        beta11, alpha11 = np.linalg.inv(A)@np.array([S11, 1])
        alpha21, beta21 = np.array([alpha[inx], beta[inx]])*S21

        beta12, alpha12 = np.linalg.inv(A)@np.array([S22, 1])
        alpha22, beta22 = np.array([alpha[inx], beta[inx]])*S12

        SS = np.array([[beta11, beta22],
                       [beta21, beta12]])@np.linalg.inv([[alpha11, alpha22],
                                                         [alpha21, alpha12]])
        
        S.append(SS) # calibrated S-parameters

    return rf.Network(s=S, frequency=netw.frequency)


def apply_cal_partial(netw, a11, a12, a21, alpha, beta, S12=None, S22=None):    
    '''
    Apply the calibration to the network netw using the error terms under partial condition.
    only S21 and S11 are calibrated, whereas S12 and S22 are assumed to be known.
    You can keep S12=None to enforce the reciprocal condition, i.e., S12 = S21.
    You can keep S22=None to enforce the symmetric condition, i.e., S22 = S11.
    '''
    S12_cal = None if S12 is None else S12*np.ones_like(netw.s[:,1,0])  # None is for S12 = S21
    S22_cal = None if S22 is None else S22*np.ones_like(netw.s[:,1,0])  # None is for S22 = S11
    
    S = []
    for inx, f in enumerate(netw.frequency.f):
        A = np.array([[a11[inx], a12[inx]],
                      [a21[inx], 1]])
        
        S11 = netw.s[inx,0,0]
        S21 = netw.s[inx,1,0]

        beta1, alpha1 = np.linalg.inv(A)@np.array([S11, 1])
        alpha2, beta2 = np.array([alpha[inx], beta[inx]])*S21
        
        # go through the four cases of S12 and S22
        if (S12_cal is None) and (S22_cal is None):
            # this is the syemmetric case, S12=S21 and S22=S11
            s21_cal = (beta2/alpha1 - (beta1/alpha1)*(alpha2/alpha1))/(1 - (alpha2/alpha1)**2)
            s11_cal = beta1/alpha1 - s21_cal*alpha2/alpha1
            S.append([[s11_cal, s21_cal],
                      [s21_cal, s11_cal]])
        elif S12_cal is None:
            # S12=S21 and S22 is assumed known
            s21_cal = beta2/alpha1 - S22_cal[inx]*alpha2/alpha1
            s11_cal = beta1/alpha1 - s21_cal*alpha2/alpha1
            S.append([[s11_cal, s21_cal],
                      [s21_cal, S22_cal[inx]]])
        elif S22_cal is None:
            # S22=S11 and S12 is assumed known
            s11_cal = beta1/alpha1 - S12_cal[inx]*alpha2/alpha1
            s21_cal = beta2/alpha1 - s11_cal*alpha2/alpha1
            S.append([[s11_cal, S12_cal[inx]],
                      [s21_cal, s11_cal]])
        else:
            # S12 and S22 are assumped known
            s21_cal = beta2/alpha1 - S22_cal[inx]*alpha2/alpha1
            s11_cal = beta1/alpha1 - S12_cal[inx]*alpha2/alpha1
            S.append([[s11_cal, S12_cal[inx]],
                      [s21_cal, S22_cal[inx]]])
    
    return rf.Network(s=S, frequency=netw.frequency)
    

if __name__ == "__main__":
    # useful functions
    c0 = 299792458   # speed of light in vacuum (m/s)
    mag2db = lambda x: 20*np.log10(abs(x))
    db2mag = lambda x: 10**(x/20)
    gamma2ereff = lambda x,f: -(c0/2/np.pi/f*x)**2
    ereff2gamma = lambda x,f: 2*np.pi*f/c0*np.sqrt(-(x-1j*np.finfo(float).eps))  # eps to ensure positive square-root
    gamma2dbmm  = lambda x: mag2db(np.exp(x.real*1e-3))  # losses dB/mm
    gamma2dbcm  = lambda x: mag2db(np.exp(x.real*1e-2))  # losses dB/cm
    time2distance = lambda x,er: x*c0/np.sqrt(er.real)

    path = os.path.dirname(os.path.realpath(__file__)) + '\\'

    # Calibration standards
    thru  = rf.Network(path + 'data\\thru.s2p')
    load  = rf.Network(path + 'data\\load.s2p')
    short = rf.Network(path + 'data\\short.s2p')
    open  = rf.Network(path + 'data\\open.s2p')

    thru_ideal  = rf.Network(path + 'data\\thru_ideal.s2p')
    load_ideal  = rf.Network(path + 'data\\load_ideal.s2p')
    short_ideal = rf.Network(path + 'data\\short_ideal.s2p')
    open_ideal  = rf.Network(path + 'data\\open_ideal.s2p')
    
    measured = [short, open, load, thru]
    ideal = [short_ideal, open_ideal, load_ideal, thru_ideal]
    
    a11, a12, a21, alpha, beta = compute_error_terms(measured, ideal)
    
    # DUT measured in forward and reverse directions
    attenF = rf.Network(path + 'data\\attenuator (forward).s2p')
    attenR = rf.Network(path + 'data\\attenuator (reverse).s2p')

    cal_dut = apply_cal(attenF, attenR, a11, a12, a21, alpha, beta)
    cal_dut_par1 = apply_cal_partial(attenF, a11, a12, a21, alpha, beta, S22=0, S12=0)  # only forward direction measurement
    cal_dut_par2 = apply_cal_partial(attenF, a11, a12, a21, alpha, beta, S22=0, S12=None)
    cal_dut_par3 = apply_cal_partial(attenF, a11, a12, a21, alpha, beta, S22=None, S12=None)

    plt.figure()
    cal_dut.s11.plot_s_db(label='Full cal')
    cal_dut_par1.s11.plot_s_db(label='Partial cal, S22=S12=0', linestyle='--')
    cal_dut_par2.s11.plot_s_db(label='Partial cal, S22=0, S12=S21', linestyle='-.')
    cal_dut_par3.s11.plot_s_db(label='Partial cal, S22=S11, S12=S21', linestyle=':')
    plt.title('Calibrated S11')
    plt.legend()

    plt.figure()
    cal_dut.s21.plot_s_db(label='Full cal')
    cal_dut_par1.s21.plot_s_db(label='Partial cal, S22=S12=0', linestyle='--')
    cal_dut_par2.s21.plot_s_db(label='Partial cal, S22=0, S12=S21', linestyle='-.')
    cal_dut_par3.s21.plot_s_db(label='Partial cal, S22=S11, S12=S21', linestyle=':')
    plt.title('Calibrated S21')
    plt.legend()

    plt.show()
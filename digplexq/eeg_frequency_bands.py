# -*- coding:utf-8 -*-

import numpy as np


#Maximum values matrix
def max_matrix(M, n):
    Mx=[]
    for j in range(len(M[0])):
        for k in range(len(M[0])):
            S=[]
            for i in range(len(M)):
                S.append(M[i][j,k])
            S = np.asarray(S)
            Mx.append(S.max())
    L = np.asarray(Mx)
    K = np.asmatrix(L)
    return K.reshape((n, n))


#Set an additional threshold to the weights
def th_weights(M, th):
    for i in range(len(M)):
        for j in range(len(M)):
            if(M[i, j] < th):
                M[i, j] = 0
            else:
                pass
    return M


#NaN to zero
def nan_to_zero(mat, n, estimator):
        
    if estimator == "pdc":
        for i in range(n):
            for j in range(n):
                for k in range(0,30):
                    if np.isnan(mat["pdc2_th"][i][j][k]) == True:
                        mat["pdc2_th"][i][j][k] = 0
                        
    elif estimator == "dtf":
        for i in range(n):
            for j in range(n):
                for k in range(0,30):
                    if np.isnan(mat["dtf2_th"][i][j][k]) == True:
                        mat["dtf2_th"][i][j][k] = 0
    return mat


#Divide by frequency           
def mat_freq(mat, f, n, estimator):
    
    M = np.zeros((n,n))
    
    if estimator == "pdc":
        for i in range(n):
            for j in range(n):
                M[i, j] = round(mat["pdc2_th"][i][j][f], 4)
                    
    elif estimator == "dtf":
        for i in range(n):
            for j in range(n):
                M[i, j] = round(mat["dtf2_th"][i][j][f], 4)
    return M


#----- EEG Frequency Bands -----

def ipdc_delta(M, n, estimator):
    mat = nan_to_zero(M, n, estimator)
    D1 = mat_freq(mat, 0, n, estimator)
    D2 = mat_freq(mat, 1, n, estimator)
    D3 = mat_freq(mat, 2, n, estimator)
    delta = [D1, D2, D3]
    D = max_matrix(delta, n)
    return D

def ipdc_theta(M, n, estimator):
    mat = nan_to_zero(M, n, estimator)
    T1 = mat_freq(mat, 3, n, estimator)
    T2 = mat_freq(mat, 4, n, estimator)
    T3 = mat_freq(mat, 5, n, estimator)
    T4 = mat_freq(mat, 6, n, estimator)
    theta = [T1, T2, T3, T4]
    T = max_matrix(theta, n)
    return T

def ipdc_alphal(M, n, estimator):
    #--Alphal (8hz-10hz):
    mat = nan_to_zero(M, n, estimator)
    Al2 = mat_freq(mat, 7, n, estimator)
    Al3 = mat_freq(mat, 8, n, estimator)
    Al4 = mat_freq(mat, 9, n, estimator)
    alphal = [Al2, Al3, Al4]
    Al = max_matrix(alphal, n)
    return Al

def ipdc_alphah(M, n, estimator):
    #--Alphah (11hz-13hz):
    mat = nan_to_zero(M, n, estimator)
    Ah2 = mat_freq(mat, 10, n, estimator)
    Ah3 = mat_freq(mat, 11, n, estimator)
    Ah4 = mat_freq(mat, 12, n, estimator)
    alphah = [Ah2, Ah3, Ah4]
    Ah = max_matrix(alphah, n)
    return Ah
    
def ipdc_alpha(M, n, estimator):
    #--Alpha (8hz-13hz):
    mat = nan_to_zero(M, n, estimator)
    Al2 = mat_freq(mat, 7, n, estimator)
    Al3 = mat_freq(mat, 8, n, estimator)
    Al4 = mat_freq(mat, 9, n, estimator)
    Ah2 = mat_freq(mat, 10, n, estimator)
    Ah3 = mat_freq(mat, 11, n, estimator)
    Ah4 = mat_freq(mat, 12, n, estimator)
    alpha = [Al2, Al3, Al4, Ah2, Ah3, Ah4]
    A = max_matrix(alpha, n)
    return A

def ipdc_betal(M, n, estimator):
    #--Betal (14-20Hz):
    mat = nan_to_zero(M, n, estimator)
    M2 = mat_freq(mat, 13, n, estimator)
    M3 = mat_freq(mat, 14, n, estimator)
    M4 = mat_freq(mat, 15, n, estimator)
    M5 = mat_freq(mat, 16, n, estimator)
    M6 = mat_freq(mat, 17, n, estimator)
    M7 = mat_freq(mat, 18, n, estimator)
    M8 = mat_freq(mat, 19, n, estimator)
    betal = [M2, M3, M4, M5, M6, M7, M8]
    Bl = max_matrix(betal, n)
    return Bl


def ipdc_betah(M, n, estimator):
    #High Beta (21–30 Hz, "Beta 3"):
    mat = nan_to_zero(M, n, estimator)
    B2 = mat_freq(mat, 20, n, estimator)
    B3 = mat_freq(mat, 21, n, estimator)
    B4 = mat_freq(mat, 22, n, estimator)
    B5 = mat_freq(mat, 23, n, estimator)
    B6 = mat_freq(mat, 24, n, estimator)
    B7 = mat_freq(mat, 25, n, estimator)
    B8 = mat_freq(mat, 26, n, estimator)
    B9 = mat_freq(mat, 27, n, estimator)
    B10 = mat_freq(mat, 28, n, estimator)
    B11 = mat_freq(mat, 29, n, estimator)
    betah = [B2, B3, B4, B5, B6, B7, B8, B9, B10, B11]
    Bh = max_matrix(betah, n)
    return Bh

def ipdc_beta(M, n, estimator):
    #--Beta (14-20Hz):
    mat = nan_to_zero(M, n, estimator)
    B1 = mat_freq(mat, 13, n, estimator)
    B2 = mat_freq(mat, 14, n, estimator)
    B3 = mat_freq(mat, 15, n, estimator)
    B4 = mat_freq(mat, 16, n, estimator)
    B5 = mat_freq(mat, 17, n, estimator)
    B6 = mat_freq(mat, 18, n, estimator)
    B7 = mat_freq(mat, 19, n, estimator)
    B8 = mat_freq(mat, 20, n, estimator)
    B9 = mat_freq(mat, 21, n, estimator)
    B10 = mat_freq(mat, 22, n, estimator)
    B11 = mat_freq(mat, 23, n, estimator)
    B12 = mat_freq(mat, 24, n, estimator)
    B13 = mat_freq(mat, 25, n, estimator)
    B14 = mat_freq(mat, 26, n, estimator)
    B15 = mat_freq(mat, 27, n, estimator)
    B16 = mat_freq(mat, 28, n, estimator)
    B17 = mat_freq(mat, 29, n, estimator)
    beta = [B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12,B13,B14,B15,B16,B17]
    B = max_matrix(beta, n)
    return B


def ipdc_wide_band(M, n, estimator):
    mat = nan_to_zero(M, n, estimator)
    D1 = mat_freq(mat, 0, n, estimator)
    D2 = mat_freq(mat, 1, n, estimator)
    D3 = mat_freq(mat, 2, n, estimator)
    T1 = mat_freq(mat, 3, n, estimator)
    T2 = mat_freq(mat, 4, n, estimator)
    T3 = mat_freq(mat, 5, n, estimator)
    T4 = mat_freq(mat, 6, n, estimator)
    Al2 = mat_freq(mat, 7, n, estimator)
    Al3 = mat_freq(mat, 8, n, estimator)
    Al4 = mat_freq(mat, 9, n, estimator)
    Ah2 = mat_freq(mat, 10, n, estimator)
    Ah3 = mat_freq(mat, 11, n, estimator)
    Ah4 = mat_freq(mat, 12, n, estimator)
    M2 = mat_freq(mat, 13, n, estimator)
    M3 = mat_freq(mat, 14, n, estimator)
    M4 = mat_freq(mat, 15, n, estimator)
    M5 = mat_freq(mat, 16, n, estimator)
    M6 = mat_freq(mat, 17, n, estimator)
    M7 = mat_freq(mat, 18, n, estimator)
    M8 = mat_freq(mat, 19, n, estimator)
    B2 = mat_freq(mat, 20, n, estimator)
    B3 = mat_freq(mat, 21, n, estimator)
    B4 = mat_freq(mat, 22, n, estimator)
    B5 = mat_freq(mat, 23, n, estimator)
    B6 = mat_freq(mat, 24, n, estimator)
    B7 = mat_freq(mat, 25, n, estimator)
    B8 = mat_freq(mat, 26, n, estimator)
    B9 = mat_freq(mat, 27, n, estimator)
    B10 = mat_freq(mat, 28, n, estimator)
    B11 = mat_freq(mat, 29, n, estimator)
    total = [D1, D2, D3, T1, T2, T3, T4, Al2, Al3, Al4, Ah2, Ah3, Ah4, M2, M3, M4, M5, M6, M7, M8, 
             B2, B3, B4, B5, B6, B7, B8, B9, B10, B11]
    Total = max_matrix(total, n)
    return Total



#-------------------------------------
#Convert MATLAB matrix to NumPy matrix:

def resh(s, m, n, l):
    #l = 6 or 14
    S=''
    for i in range(len(s)):
        if s[i] == ' ' and s[i+1] == ' ':
            x = s.replace(' ', '')
    
    for j in range(m*n):
        S = S+x[l*j : l*(j+1)]+' '
    return S


def to_np_matrix(str, m, n):
    X=[]
    arr = str.split()
    for i in range(m*n):
        X.append(float(arr[i]))
        
    Y = np.array(X)
    M = Y.reshape((m,n))
    M = np.asmatrix(M)
    return M


def convert_to_np(A11, A12, A13):
    a1 = to_np_matrix(resh(A11, 24, 8, 6), 24, 8)
    a2 = to_np_matrix(resh(A12, 24, 8, 6), 24, 8)
    a3 = to_np_matrix(resh(A13, 24, 8, 6), 24, 8)
    M = np.hstack((a1, a2, a3))
    return M

def convert_to_np_gct(A11, A12, A13):
    a1 = to_np_matrix(resh(A11, 24, 9, 6), 24, 9)
    a2 = to_np_matrix(resh(A12, 24, 9, 6), 24, 9)
    a3 = to_np_matrix(resh(A13, 24, 6, 6), 24, 6)
    M = np.hstack((a1, a2, a3))
    return M


#if NaN = 0.0001 use adj_zeros(convert_to_np(A11, A12, A13)).
def adj_zeros(M):
    for i in range(len(M)):
        for j in range(len(M)):
            if(M[i, j] == 0.0001):
                M[i, j] = 0.0000
            else:
                pass
    return M


'''
#------------------------
#Divide by frequency band:

def ipdc_delta(A11, A12, A13, A21, A22, A23, A31, A32, A33):
    #--Delta (1hz-3hz):
    D1 = convert_to_np(A11, A12, A13)
    D2 = convert_to_np(A21, A22, A23)
    D3 = convert_to_np(A31, A32, A33)
    #D4 = adj_zeros(convert_to_np(A41, A42, A43))
    delta = [D1, D2, D3]
    D = max_matrix(delta, 24)
    return D


def ipdc_theta(A41, A42, A43, A51, A52, A53, A61, A62, A63, A71, A72, A73):
    #--Theta (4hz-7hz):
    T1 = convert_to_np(A41, A42, A43)
    T2 = convert_to_np(A51, A52, A53)
    T3 = convert_to_np(A61, A62, A63)
    T4 = convert_to_np(A71, A72, A73)
    theta = [T1, T2, T3, T4]
    T = max_matrix(theta, 24)
    return T
    

def ipdc_alphal(A81, A82, A83, A91, A92, A93, A101, A102, A103):
    #--Alphal (8hz-10hz):
    #Al1 = adj_zeros(convert_to_np(A71, A72, A73))
    Al2 = convert_to_np(A81, A82, A83)
    Al3 = convert_to_np(A91, A92, A93)
    Al4 = convert_to_np(A101, A102, A103)
    alphal = [Al2, Al3, Al4]
    Al = max_matrix(alphal, 24)
    return Al

def ipdc_alphah(A111, A112, A113, A121, A122, A123, A131, A132, A133):
    #--Alphah (11hz-13hz):
    #Ah1 = adj_zeros(convert_to_np(A101, A102, A103))
    Ah2 = convert_to_np(A111, A112, A113)
    Ah3 = convert_to_np(A121, A122, A123)
    Ah4 = convert_to_np(A131, A132, A133)
    alphah = [Ah2, Ah3, Ah4]
    Ah = max_matrix(alphah, 24)
    return Ah

def ipdc_beta(A141, A142, A143, A151, A152, A153, A161, A162, A163, A171, A172, 
              A173, A181, A182, A183, A191, A192, A193, A201, A202, A203):
    #--Beta (14-20Hz):
    #M1 = adj_zeros(convert_to_np(A131, A132, A133))
    M2 = convert_to_np(A141, A142, A143)
    M3 = convert_to_np(A151, A152, A153)
    M4 = convert_to_np(A161, A162, A163)
    M5 = convert_to_np(A171, A172, A173)
    M6 = convert_to_np(A181, A182, A183)
    M7 = convert_to_np(A191, A192, A193)
    M8 = convert_to_np(A201, A202, A203)
    beta = [M2, M3, M4, M5, M6, M7, M8]
    B = max_matrix(beta, 24)
    return B


def ipdc_betah(A211, A212, A213, A221, A222, A223, A231, A232, A233, A241, A242, A243, A251, A252, A253,
               A261, A262, A263, A271, A272, A273, A281, A282, A283, A291, A292, A293, A301, A302, A303):
    #High Beta (21–30 Hz, "Beta 3"):
    #B1 = adj_zeros(convert_to_np(A201, A202, A203))
    B2 = convert_to_np(A211, A212, A213)
    B3 = convert_to_np(A221, A222, A223)
    B4 = convert_to_np(A231, A232, A233)
    B5 = convert_to_np(A241, A242, A243)
    B6 = convert_to_np(A251, A252, A253)
    B7 = convert_to_np(A261, A262, A263)
    B8 = convert_to_np(A271, A272, A273)
    B9 = convert_to_np(A281, A282, A283)
    B10 = convert_to_np(A291, A292, A293)
    B11 = convert_to_np(A301, A302, A303)
    betah = [B2, B3, B4, B5, B6, B7, B8, B9, B10, B11]
    Bh = max_matrix(betah, 24)
    return Bh

def ipdc_betaw(A141, A142, A143, A151, A152, A153, A161, A162, A163, A171, A172, 
               A173, A181, A182, A183, A191, A192, A193, A201, A202, A203,
               A211, A212, A213, A221, A222, A223, A231, A232, A233, A241, A242, A243, A251, A252, A253,
               A261, A262, A263, A271, A272, A273, A281, A282, A283, A291, A292, A293, A301, A302, A303):
    M2 = convert_to_np(A141, A142, A143)
    M3 = convert_to_np(A151, A152, A153)
    M4 = convert_to_np(A161, A162, A163)
    M5 = convert_to_np(A171, A172, A173)
    M6 = convert_to_np(A181, A182, A183)
    M7 = convert_to_np(A191, A192, A193)
    M8 = convert_to_np(A201, A202, A203)
    B2 = convert_to_np(A211, A212, A213)
    B3 = convert_to_np(A221, A222, A223)
    B4 = convert_to_np(A231, A232, A233)
    B5 = convert_to_np(A241, A242, A243)
    B6 = convert_to_np(A251, A252, A253)
    B7 = convert_to_np(A261, A262, A263)
    B8 = convert_to_np(A271, A272, A273)
    B9 = convert_to_np(A281, A282, A283)
    B10 = convert_to_np(A291, A292, A293)
    B11 = convert_to_np(A301, A302, A303)
    betaw = [M2, M3, M4, M5, M6, M7, M8, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11]
    Bw = max_matrix(betaw, 24)
    return Bw

def ipdc_wide_band(A11, A12, A13, A21, A22, A23, A31, A32, A33, 
                   A41, A42, A43, A51, A52, A53, A61, A62, A63, A71, A72, A73,
                   A81, A82, A83, A91, A92, A93, A101, A102, A103,
                   A111, A112, A113, A121, A122, A123, A131, A132, A133,
                   A141, A142, A143, A151, A152, A153, A161, A162, A163, A171, A172, 
                   A173, A181, A182, A183, A191, A192, A193, A201, A202, A203,
                   A211, A212, A213, A221, A222, A223, A231, A232, A233, A241, A242, A243, A251, A252, A253,
                   A261, A262, A263, A271, A272, A273, A281, A282, A283, A291, A292, A293, A301, A302, A303):
    
    D1 = convert_to_np(A11, A12, A13)
    D2 = convert_to_np(A21, A22, A23)
    D3 = convert_to_np(A31, A32, A33)
    T1 = convert_to_np(A41, A42, A43)
    T2 = convert_to_np(A51, A52, A53)
    T3 = convert_to_np(A61, A62, A63)
    T4 = convert_to_np(A71, A72, A73)
    Al2 = convert_to_np(A81, A82, A83)
    Al3 = convert_to_np(A91, A92, A93)
    Al4 = convert_to_np(A101, A102, A103)
    Ah2 = convert_to_np(A111, A112, A113)
    Ah3 = convert_to_np(A121, A122, A123)
    Ah4 = convert_to_np(A131, A132, A133)
    M2 = convert_to_np(A141, A142, A143)
    M3 = convert_to_np(A151, A152, A153)
    M4 = convert_to_np(A161, A162, A163)
    M5 = convert_to_np(A171, A172, A173)
    M6 = convert_to_np(A181, A182, A183)
    M7 = convert_to_np(A191, A192, A193)
    M8 = convert_to_np(A201, A202, A203)
    B2 = convert_to_np(A211, A212, A213)
    B3 = convert_to_np(A221, A222, A223)
    B4 = convert_to_np(A231, A232, A233)
    B5 = convert_to_np(A241, A242, A243)
    B6 = convert_to_np(A251, A252, A253)
    B7 = convert_to_np(A261, A262, A263)
    B8 = convert_to_np(A271, A272, A273)
    B9 = convert_to_np(A281, A282, A283)
    B10 = convert_to_np(A291, A292, A293)
    B11 = convert_to_np(A301, A302, A303)

    total = [D1, D2, D3, T1, T2, T3, T4, Al2, Al3, Al4, Ah2, Ah3, Ah4, M2, M3, M4, M5, M6, M7, M8, 
             B2, B3, B4, B5, B6, B7, B8, B9, B10, B11]
    Total = max_matrix(total, 24)
    return Total
    '''
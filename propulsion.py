# -*- coding: utf-8 -*-
#Created on Wed Apr 21 23:01:41 2021

#@author: Thiago Novaes

# include

import os
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

#user

dir_HS = "HS serie"   #diretorio local com os HS
dir_B = "B serie"     #diretorio local com os B

interpol = "linear"      #tipo de interpolação desejado


v = 11.9       #velocidade

D_min = 2.5   #diametro minimo 
D_max = 2.5   #diametro maximo
D_num = 1   #numero de subdivisoes

n_min = 200/60  #rotação minima
n_max = 400/60  #rotação maxima
n_num = 11       #numero de subdivisoes

helices = 4
T = 430800      #força requerida
tol = 0.5       #tolerancia da força requerida em %; float("inf")
rho = 1025      #densidade do fluido

serie = 'HS' # B or HS

h = 3.4 # altura do helice

Rn = 5E9           #B serie only
corr_reynolds = False #B serie only 

output = 'final.csv' #path do output
sigma = 1.81
eixo = 0.99
T = T/helices

# function

sigmas = [0.25,0.3,0.4,0.5,0.6,0.75,1,2.5]
sigmas_str = ["0,25","0,30","0,40","0,50","0,60","0,75","1,00","2,50"]

def lim_inter(sigma):
    for i in range(len(sigmas)):
        if sigma < sigmas[i]:
            return i-1, i 
a,b = lim_inter(sigma)

def J(v,n,D):
    return v/(n*D)

def eta(j,Kt,Kq):
    return j*Kt/(2*np.pi*Kq)

def Trust(Kt,n,D):
    return Kt*(rho*n**2*D**4)

def Torque(Kq,n,D):
    return Kq*(rho*n**2*D**5)

def cavitation(D,u,t,v,n,h,T,rho):
    Ap = np.pi*D**2/4*u*np.cos(np.arctan(t))           
    Vs2 = v**2+(np.pi*n*0.7*D)**2
    sigma = (100000 + 10000*h)/(rho/2*Vs2)
    cav = T/Ap/(rho/2*Vs2)/(0.42*sigma**(0.7))
    return cav, sigma

def main(v,D_min,D_max,D_num,n_min,n_max,n_num,T,tol,serie,output,interpol = "linear",h = 0, Rn = 0,corr_reynolds = False,sigma =0):
    out = pd.DataFrame(index=range(0),columns=range(15))
    out.columns = ["Blades","BAR","P/D", "sigma","D","n", "J", "Kt", "Kq", "Trust", "Q", "Eta","THP","DHP","BHP"]
    for D in np.linspace(D_min,D_max,D_num):
        for n in np.linspace(n_min,n_max,n_num):
            j = J(v,n,D)
            if serie == 'HS':
                out = HSsigma(j,v,D,n,T,tol,interpol,out,sigma)
            elif serie == 'B' and 0 <= j <= 1.6:
                out = B(j,v,D,n,T,tol,interpol,out,h,Rn,corr_reynolds)
            else:
                pass
    out.sort_values(by=["BHP"],ascending=True, inplace=True)
    out.to_csv(output,sep=";")
    return out
                
def HS(j,v,D,n,T,tol,interpol,out):
    for filename in os.listdir(dir_HS):
        df = pd.read_csv(os.path.join(dir_HS,filename))
        if df.head(1)["J"].values <= j <= df.tail(1)["J"].values:
            Kt = interp1d(df["J"], df["Kt"],kind =interpol)(j)
            Kq = interp1d(df["J"], df["Kq"],kind =interpol)(j)
            Eta = interp1d(df["J"], df["Eta"],kind =interpol)(j)
            out = validation_HS(j,v,D,n,T,tol,Kt,Kq,Eta,out,filename)
        else:
            pass
    return out

def HSsigma(j,v,D,n,T,tol,interpol,out,sigma):
    for filename1 in os.listdir(dir_HS):
        if filename1[13:17] == sigmas_str[a]:
            for filename2 in os.listdir(dir_HS):
                if filename2[:17] == filename1[:13] + sigmas_str[b]:
                    df1 = pd.read_csv(os.path.join(dir_HS,filename1))
                    df2 = pd.read_csv(os.path.join(dir_HS,filename2))
                    if df1.head(1)["J"].values <= j <= df1.tail(1)["J"].values and df2.head(1)["J"].values <= j <= df2.tail(1)["J"].values:
                        Kt1 = interp1d(df1["J"], df1["Kt"],kind =interpol)(j)
                        Kq1 = interp1d(df1["J"], df1["Kq"],kind =interpol)(j)
                        Eta1 = interp1d(df1["J"], df1["Eta"],kind =interpol)(j)
                        
                        Kt2 = interp1d(df2["J"], df2["Kt"],kind =interpol)(j)
                        Kq2 = interp1d(df2["J"], df2["Kq"],kind =interpol)(j)
                        Eta2 = interp1d(df2["J"], df2["Eta"],kind =interpol)(j)
                        
                        Kt = interp1d(sigmas[a:a+2], [Kt1,Kt2])(sigma)
                        Kq = interp1d(sigmas[a:a+2], [Kq1,Kq2])(sigma)
                        Eta = interp1d(sigmas[a:a+2], [Eta1,Eta2])(sigma)
                        out = validation_HS(j,v,D,n,T,tol,Kt,Kq,Eta,out,filename1,sigma)
                    else:
                        pass
    return out

def validation_HS(j,v,D,n,T,tol,Kt,Kq,Eta,out,filename,sigma):
    trust = Trust(Kt,n,D)
    Q = Torque(Kq,n,D)
    THP = v*trust
    DHP = THP/Eta
    BHP = DHP/eixo
    if T <= trust <=T*(1+tol):
        blades, bar , pd , lixo, ext = filename.split("-")
        new_row = {"Blades": blades,"BAR":bar,"P/D":pd,"sigma":sigma,"D":D,"n":n, "J":j, "Kt":Kt,"Kq":Kq, "Trust":trust,"Q":Q ,"Eta":Eta,"THP":THP,"DHP":DHP,"BHP":BHP}
        out = out.append(new_row, ignore_index=True)
    else:
        pass
    return out

def B(j,v,D,n,T,tol,interpol,out,h,Rn,corr_reynolds):
    df_kt = pd.read_csv(os.path.join(dir_B,'kt.csv'))
    df_kq = pd.read_csv(os.path.join(dir_B,'kq.csv'))
    df_dkt = pd.read_csv(os.path.join(dir_B,'dkt.csv'))
    df_dkq = pd.read_csv(os.path.join(dir_B,'dkq.csv'))
    for t in np.linspace(0.5,1.4,10):
        for u in np.linspace(0.3,1.05,16):
            for z in np.linspace(2,7,6):
                Kt = (df_kt["coef"]*j**df_kt["s"]*t**df_kt["t"]*u**df_kt["u"]*z**df_kt["v"]).sum() 
                Kq = (df_kq["coef"]*j**df_kq["s"]*t**df_kq["t"]*u**df_kq["u"]*z**df_kq["v"]).sum()
                if corr_reynolds == True and Rn > 2E6:
                    R = np.log10(Rn)-0.301
                    Kt += (df_dkt["coef"]*R**df_dkt["r"]*j**df_dkt["s"]*t**df_dkt["t"]*u**df_dkt["u"]*z**df_dkt["v"]).sum()
                    Kq += (df_dkq["coef"]*R**df_dkq["r"]*j**df_dkq["s"]*t**df_dkq["t"]*u**df_dkq["u"]*z**df_dkq["v"]).sum()
                Eta = eta(j,Kt,Kq)
                out = validation_B(j,v,D,n,T,tol,Kt,Kq,Eta,out,t,u,z,h)
    return out

def validation_B(j,v,D,n,T,tol,Kt,Kq,Eta,out,t,u,z,h):
    cav, sigma = cavitation(D,u,t,v,n,h,T,rho)
    trust = Trust(Kt,n,D)
    Q = Torque(Kq,n,D)
    THP = v*trust
    DHP = THP/Eta
    BHP = DHP/eixo
    if T <= trust <=T*(1+tol) and Kt > 0 and Kq > 0 and cav < 1:
        new_row = {"Blades": z,"BAR":u,"P/D":t,"sigma":sigma,"D":D,"n":n, "J":j, "Kt":Kt,"Kq":Kq, "Trust":trust,"Q":Q,"Eta":Eta,"THP":THP,"DHP":DHP,"BHP":BHP}
        out = out.append(new_row, ignore_index=True)
    else:
        pass
    return out

out = main(v,D_min,D_max,D_num,n_min,n_max,n_num,T,tol,serie,output,interpol,h,Rn,corr_reynolds,sigma)

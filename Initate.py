import OutTxt
import numpy as np
from numpy import cos,sin,angle,arccos,arctan,sqrt
#---------------极坐标形式----------------#
def PolarU(U_Amplitude,U_Angle,**Option):          
    U_Complex = U_Amplitude*cos(U_Angle)+U_Amplitude*sin(U_Angle)*1j
    U = abs(U_Complex)
    Angle = angle(U_Complex)
    # if 'string' not in Option:
    #     Option['string'] = '\n电压初始化为：\n'
    #     OutTxt.Complex(U_Complex,**Option)    
    return(U,Angle)
#---------------直角坐标形式--------------#    
def RectanU(U_Amplitude,U_Angle,**Option):
    U_e = U_Amplitude*cos(U_Angle)
    U_f = U_Amplitude*cos(U_Angle)
    # if 'string' not in Option:
    #     Option['string'] = '\n 电压实部初始化为：\n'
    #     OutTxt.Real(U_e,**Option)
    #     Option['string'] = '\n 电压虚部初始化为：\n'
    #     OutTxt.Real(U_f,**Option)
    return(U_e,U_f)
#--------------------------LCC初始化-------------------------#
def LCC_Polar(DC_NodeData,**Option):
    Nc = DC_NodeData.shape[0] # LCC换流器
    Vd = DC_NodeData[:,3]
    Id = DC_NodeData[:,7]
    W = np.cos(DC_NodeData[:,4]) # 取换流器控制角(rad)
    Kt = DC_NodeData[:,11]
    fi = np.ones([Nc,1])*arccos(0.9)
    # if 'string' not in Option:
    #     Option['string'] = '直流电压初始化为：\n'
    #     OutTxt.Real(Vd,**Option)
    #     Option['string'] = '直流电流初始化为：\n'
    #     OutTxt.Real(Id,**Option)
    #     Option['string'] = '变比初始化为：\n'
    #     OutTxt.Real(Kt,**Option)  
    #     Option['string'] = '控制角初始化为：\n'
    #     OutTxt.Real(57.3*arccos(W),**Option) 
    #     Option['string'] = '功率因数初始化为：\n'
    #     OutTxt.Real(cos(fi),**Option)   
    return(Vd,Id,Kt,W,fi) 
#--------------------------VSC初始化-------------------------#
def VSC_Polar(DC_NodeData,**Option):          
    Psi=DC_NodeData[:,3]                                      # 换流站等效
    Qsi=DC_NodeData[:,4]
    Usi=DC_NodeData[:,5]
    Udc=DC_NodeData[:,6]
    XL=DC_NodeData[:,11]
    N =DC_NodeData[:,12] 
    Idc=Psi/(N*Udc)  
    derta=arctan(Psi/(N*Usi*Usi/XL-N*Qsi))
    M=2*sqrt(2)*Psi*XL/(sqrt(3)*Usi*N*Udc*sin(derta))
    # if 'string' not in Option:
    #     Option['string'] = '直流电压初始化为：\n'
    #     OutTxt.Real(Udc,**Option)
    #     Option['string'] = '直流电流初始化为：\n'
    #     OutTxt.Real(Idc,**Option)
    #     Option['string'] = '功角初始化为：\n'
    #     OutTxt.Real(57.3*derta,**Option)
    #     Option['string'] = '调制比初始化为：\n'
    #     OutTxt.Real(M,**Option)     
    return(Udc,Idc,derta,M)


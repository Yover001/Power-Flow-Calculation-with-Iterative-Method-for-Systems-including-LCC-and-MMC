import datetime
import math
import copy
from File1 import NetData
import numpy as np
from math import pi,sqrt
from numpy import arctan,sin,cos,arccos
from File1.NetData import GetY,GetYdc
from Initate import PolarU,LCC_Polar,VSC_Polar
from AC_NR import PolarNR
from DC_NR import DC_PolarNR 
from OutTxt import Real,SingleTxt,StringTxt
Out_Path = 'out.txt'
starttime = datetime.datetime.now()
#----------------------------网络信息读取-------------------------------#
# 获取直流系统信息
RootPath = 'C:\\Users\\lenovo\\Desktop\\交替迭代 - 双平衡节点\\Data\\'
VSC_Node = RootPath+'VSC_NodeData.txt'
LCC_Node = RootPath+'LCC_NodeData.txt'
DC_Line = RootPath+'DC_LineData.txt'
VSC_NodeData = NetData.VSC_GetNodeData(VSC_Node,show=1)
LCC_NodeData = NetData.LCC_GetNodeData(LCC_Node,show=1)
DC_LineData = NetData.DC_GetLineData(DC_Line,show=1)
Ydc=GetYdc(VSC_NodeData,LCC_NodeData,DC_LineData,path=Out_Path,width=6)
LCC=copy.copy(LCC_NodeData)                                             # 记录矩阵初始值
Ps=copy.copy(VSC_NodeData[:,3])                                         #换流器控制值(迭代过程中不变)
Qs=copy.copy(VSC_NodeData[:,4])
Num_DC = Ydc.shape[0]                                                   # 直流节点数目
AC_LCC=LCC_NodeData[:,0].astype(int)                                    #与直流相连的交流节点(LCC)
AC_VSC=VSC_NodeData[:,0].astype(int)                                    #与直流相连的交流节点(VSC)
# 获取交流系统信息
FilePath_Node = RootPath+'NodeData.txt'
FilePath_Line = RootPath+'LineData.txt'
NodeData = NetData.GetNodeData(FilePath_Node,show=1)
LineData = NetData.GetLineData(FilePath_Line,show=1)
Y = GetY(NodeData,LineData,path=Out_Path,width=6)
Data=copy.copy(NodeData)
StringTxt(path=Out_Path,string='交直流混合系统潮流计算',fmt='w')
#-----------------------------交替迭代计算----------------------------------#
L_V = 0
Tol = 1e-5
while True:
    L_V=L_V+1
    print('交替迭代次数：',L_V)
    S0 = '-'
    for i in range(130):
        S0 = S0+'-' 
    SingleTxt(L_V,path=Out_Path,string='\n交替迭代次数：')
    StringTxt(path=Out_Path,string=S0)
#------------------------------读取VSC直流数据---------------------------------# 
    UNode,QNode,PHNode = NetData.VSC_GetNetData(VSC_NodeData,DC_LineData,path=Out_Path)
    # print('平衡换流器为：',PHNode)
    NodeData = NetData.GetNodeData(FilePath_Node)                       # 重置交流节点矩阵
    for i in QNode:
        NodeData[int(i-1),1]=1
        NodeData[int(i-1),2]=NodeData[int(i-1),2]+VSC_NodeData[np.where(VSC_NodeData[:,0]==i)[0][0],3]
        NodeData[int(i-1),3]=NodeData[int(i-1),3]+VSC_NodeData[np.where(VSC_NodeData[:,0]==i)[0][0],4] 
    if np.any(UNode):
        for i in UNode:
            NodeData[int(i-1),1]=2
            NodeData[int(i-1),2]=NodeData[int(i-1),2]+VSC_NodeData[np.where(VSC_NodeData[:,0]==i),3]
            NodeData[int(i-1),8]=VSC_NodeData[np.where(VSC_NodeData[:,0]==i),5] 
#------------------------------读取LCC直流数据---------------------------------# 
    Nc = LCC_NodeData.shape[0]                                           # LCC换流器
    for i in range(Nc):
        NodeData[int(LCC_NodeData[i,0])-1,2]=NodeData[int(LCC_NodeData[i,0])-1,2]+LCC_NodeData[i,5]
        NodeData[int(LCC_NodeData[i,0])-1,3]=NodeData[int(LCC_NodeData[i,0])-1,3]+LCC_NodeData[i,6]
#------------------------------交流潮流计算---------------------------------# 
    U,Angle = PolarU(NodeData[:,8],NodeData[:,9],path=Out_Path,width=9)  # 交流初始化
    Iter = 0
    MaxIter = 10
    PQNode,PVNode,SlackNode,P_Real,Q_Real = NetData.GetNetData(NodeData,LineData,path=Out_Path)
    while True:
        Iter = Iter + 1
        U,Angle,MaxError,P,Q = PolarNR(U,Angle,Y,PQNode,PVNode,SlackNode,P_Real,Q_Real,Tol,path=Out_Path,width=9)
        if Iter>MaxIter or MaxError<Tol:
            break
        # 结束交流循环
    if MaxError<Tol:
        SingleTxt(Iter-1,path=Out_Path,string=S0+'\n交流迭代完成，更新次数为：')
        SingleTxt(MaxError,path=Out_Path,string='交流最大误差为：')
        # Real(U,path=Out_Path,string=S0+'\n电压幅值为：\n')
        # Real(Angle,path=Out_Path,string='相角为：\n')
    else:
        SingleTxt(MaxError,path=Out_Path,string='交流结果不收敛!')
# --------------------------- 交流数据 P,Q,U,angle传输给直流(VSC,LCC)-------------------#
    NodeData[:,2]=P.squeeze(-1)
    NodeData[:,3]=Q.squeeze(-1)
    NodeData[:,8]=U
    NodeData[:,9]=Angle
    for i in AC_VSC:
        VSC_NodeData[np.where(VSC_NodeData[:,0]==i),3]=-NodeData[i-1,2]-Data[i-1,2]+Data[i-1,4] #传输功率
        VSC_NodeData[np.where(VSC_NodeData[:,0]==i),4]=-NodeData[i-1,3]-Data[i-1,3]+Data[i-1,5]
        VSC_NodeData[np.where(VSC_NodeData[:,0]==i),5]=NodeData[i-1,8]
        VSC_NodeData[np.where(VSC_NodeData[:,0]==i),9]=NodeData[i-1,9]
    PAC=copy.copy(VSC_NodeData[:,3])
    QAC=copy.copy(VSC_NodeData[:,4])
    for i in AC_LCC:
        LCC_NodeData[np.where(LCC_NodeData[:,0]==i),5]=-NodeData[i-1,2]-Data[i-1,2]+Data[i-1,4] #传输功率
        LCC_NodeData[np.where(LCC_NodeData[:,0]==i),6]=-NodeData[i-1,3]-Data[i-1,3]+Data[i-1,5]
    PLC=copy.copy(LCC_NodeData[:,5])
    QLC=copy.copy(LCC_NodeData[:,6])
#-----------------------------------（LCC-MMC）直流潮流计算----------------------------#
    DC_Iter = 0
    DC_MaxIter = 5
    # StringTxt(path=Out_Path,string=S0+'\nLCC直流初始化参数：\n')
    Vd,Id,kt,W,fi = LCC_Polar(LCC_NodeData,path=Out_Path,width=9)                   # 直流初始化
    # StringTxt(path=Out_Path,string=S0+'\nVSC直流初始化参数：\n')
    Udc,Idc,derta,M = VSC_Polar(VSC_NodeData,path=Out_Path,width=9)                 # 直流初始化
    while True:
        DC_Iter = DC_Iter+1                             
        Vd,Id,kt,W,fi,P_LCC,Q_LCC,Udc,Idc,derta,M,DC_MaxError,P_VSC,Q_VSC = DC_PolarNR(Vd,Id,kt,W,fi,Udc,Idc,derta,M,Ydc,NodeData,VSC_NodeData,LCC_NodeData,LCC,Ps,Qs,Tol,path=Out_Path)
        if DC_Iter>DC_MaxIter or DC_MaxError<Tol:
            break
        # 结束直流循环
    if DC_MaxError<Tol:
        StringTxt(path=Out_Path,string=S0)
        SingleTxt(DC_Iter-1,path=Out_Path,string='直流迭代完成，更新次数为：')
        SingleTxt(DC_MaxError,path=Out_Path,string='直流最大误差为：')
#------------------------------LCC计算结果-------------------------------------------#
        # Real(Vd,path=Out_Path,string=S0+'\nLCC直流电压：\n')
        # Real(Id,path=Out_Path,string='LCC直流电流：\n')
        # Real(kt,path=Out_Path,string='换流变变比：\n')
        # Real(57.3*arccos(W),path=Out_Path,string='控制角：\n')
        # Real(cos(fi),path=Out_Path,string='功率因数：\n')
#------------------------------VSC计算结果-------------------------------------------#
        # Real(Udc,path=Out_Path,string=S0+'\nVSC的直流电压为：\n')
        # Real(Idc,path=Out_Path,string='VSC的直流电流为：\n')
        # Real(57.3*derta,path=Out_Path,string='功角：\n')
        # Real(M,path=Out_Path,string='调制比：\n')
    else:
        SingleTxt(DC_MaxError,path=Out_Path,string='直流结果不收敛!')   
#-----------------------直流结果(P,Q)写入DC_Nodedata传入交流-----------------------#
    P_VSC=P_VSC.squeeze(-1)
    Q_VSC=Q_VSC.squeeze(-1)
    VSC_NodeData[:,3]=P_VSC
    VSC_NodeData[:,4]=Q_VSC
    P_LCC=P_LCC.reshape(-1)
    Q_LCC=Q_LCC.reshape(-1)
    LCC_NodeData[:,3]=Vd
    LCC_NodeData[:,5]=P_LCC
    LCC_NodeData[:,6]=Q_LCC
# #---------------------------------------收敛判据------------------------------------#
    D_PVSC=np.max(np.abs(PAC-P_VSC))
    D_QVSC=np.max(np.abs(QAC-Q_VSC))                                                # VSC交界处功率差值
    D_PLCC=np.max(np.abs(PLC-P_LCC))
    D_QLCC=np.max(np.abs(QLC-Q_LCC))
    # print('交替迭代偏差：',D_PVSC,D_QVSC,D_PLCC,D_QLCC)                           # LCC交界处功率差值
    AD_MaxError=max(D_PVSC,D_QVSC,D_PLCC,D_QLCC)
    print(AD_MaxError) 
    if L_V>7 or AD_MaxError<Tol:
        break
if AD_MaxError<Tol:
    SingleTxt(L_V,path=Out_Path,string='交替迭代迭代完成，更新次数为：')
    SingleTxt(AD_MaxError,path=Out_Path,string='最大误差为：')
    Real(U,path=Out_Path,string=S0+'\n电压幅值为：\n')
    Real(Angle,path=Out_Path,string='相角为：\n')
    Real(Vd,path=Out_Path,string=S0+'\nLCC直流电压：\n')
    Real(Id,path=Out_Path,string='LCC直流电流：\n')
    Real(kt,path=Out_Path,string='换流变变比：\n')
    Real(57.3*arccos(W),path=Out_Path,string='控制角：\n')
    Real(cos(fi),path=Out_Path,string='功率因数：\n')
    Real(Udc,path=Out_Path,string=S0+'\nVSC的直流电压为：\n')
    Real(Idc,path=Out_Path,string='VSC的直流电流为：\n')
    Real(57.3*derta,path=Out_Path,string='功角：\n')
    Real(M,path=Out_Path,string='调制比：\n')
else:
    SingleTxt(AD_MaxError,path=Out_Path,string='结果不收敛!')
# print('交替迭代偏差：',AD_MaxError)
# TIME
endtime = datetime.datetime.now()
print (endtime - starttime)
np. savetxt(r"C:\Users\lenovo\Desktop\python代码\XML_e\对比数据\Angle.csv",Angle,delimiter=',')
np. savetxt(r"C:\Users\lenovo\Desktop\python代码\XML_e\对比数据\U.csv",U,delimiter=',')
np. savetxt(r"C:\Users\lenovo\Desktop\python代码\XML_e\对比数据\LCC_Vd.csv",Vd,delimiter=',')
np. savetxt(r"C:\Users\lenovo\Desktop\python代码\XML_e\对比数据\LCC_Id.csv",Id,delimiter=',')
np. savetxt(r"C:\Users\lenovo\Desktop\python代码\XML_e\对比数据\VSC_Udc.csv",Udc,delimiter=',')
np. savetxt(r"C:\Users\lenovo\Desktop\python代码\XML_e\对比数据\VSC_Idc.csv",Idc,delimiter=',')

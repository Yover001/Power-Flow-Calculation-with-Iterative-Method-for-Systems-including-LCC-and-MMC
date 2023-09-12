import math
import OutTxt
import numpy as np
import scipy.linalg  # 矩阵对角组合
from OutTxt import Real
from math import pi,sqrt
from numpy import arctan,sin,cos,arccos
# ------------------------牛顿法计算---------------------------#
def DC_PolarNR(Vd,Id,Kt,W,fi,Udc,Idc,derta,M,Ydc,NodeData,VSC_NodeData,LCC_NodeData,LCC,Ps,Qs,Tol,**Option):
#---------------------------LCC不平衡量------------------------#
	Ids = LCC[:,7]
	Vds = LCC[:,3]
	Pds = LCC[:,5]
	Ws = np.cos(LCC[:,4])
	kts = LCC[:,11]	
	kr = 0.995
	DC_Node = LCC_NodeData[:,0]
	X = LCC_NodeData[:,8]
	Vt = np.zeros([len(DC_Node),1])
	for i in range(len(DC_Node)):
		Vt[i] = NodeData[int(DC_Node[i]-1),8]
	Nc = LCC_NodeData.shape[0] 					# LCC换流器
	N_LCC = LCC_NodeData[:,12]                  # 换流器组数
	N_VSC = VSC_NodeData[:,12]                 
	N = np.append(N_LCC,N_VSC)
	control1 = LCC_NodeData[:,9]
	control2 = LCC_NodeData[:,10]
	Delta_D1 = np.zeros([Nc,1]) # 初始化
	Delta_D2 = np.zeros([Nc,1])
	Delta_D3 = np.zeros([Nc,1])
	Delta_D4 = np.zeros([Nc,1])
	Delta_D5 = np.zeros([Nc,1])
	control1 = np.int64(control1)
	control2 = np.int64(control2)
	Uc=np.append(Vd,Udc)
	Pdc = np.zeros([Nc,1])
	Qdc = np.zeros([Nc,1])
	for i in range(Nc): 
		Delta_D1[i] = Vd[i]-2.7*Kt[i]*Vt[i]*W[i]+1.9*X[i]*Id[i]  # I为标量
		Delta_D2[i] = Vd[i]-2.7*kr*Kt[i]*Vt[i]*cos(fi[i])
		if LCC_NodeData[i,2]==1:
			Pdc[i] = N_LCC[i]*Vd[i]*Id[i] # 对交流的影响
			Qdc[i] = N_LCC[i]*Vd[i]*Id[i]*np.tan(fi[i])
			Delta_D3[i] = Id[i]-N_LCC[i]*np.sum(Ydc[i,:]*Uc)
		elif LCC_NodeData[i,2]==2:
			Pdc[i] = -N_LCC[i]*Vd[i]*Id[i] # 对交流的影响  #8月31日 加了负号
			Qdc[i] = N_LCC[i]*Vd[i]*Id[i]*np.tan(fi[i])   #8月31日 去了负号
			Delta_D3[i] = -Id[i]-N_LCC[i]*np.sum(Ydc[i,:]*Uc)
		else:
			print("直流节点类型参数输入错误！")     
		if control1[i]==1:
			Delta_D4[i] = Id[i]-Ids[i] # 定电流
		elif control1[i]==2:
			Delta_D4[i] = Vd[i]-Vds[i] # 定电压
		elif control1[i]==3:
			Delta_D4[i] = Vd[i]*Id[i]-Pds[i] # 定功率
		elif control1[i]==4:
			Delta_D4[i] = W[i]-Ws[i] # 定控制角
		elif control1[i]==5:
			Delta_D4[i] = Kt[i]-kts[i] # 定变比
		else:
			print("直流节点控制类型参数输入错误！")
		if control2[i]==1:
			Delta_D5[i] = Id[i]-Ids[i]
		elif control2[i]==2:
			Delta_D5[i] = Vd[i]-Vds[i]
		elif control2[i]==3:
			Delta_D5[i] = Vd[i]*Id[i]-Pds[i]
		elif control2[i]==4:
			Delta_D5[i] = W[i]-Ws[i]
		elif control2[i]==5:
			Delta_D5[i]= Kt[i]-kts[i]
		else:
			print("直流节点控制类型参数输入错误！")
#---------------------------VSC不平衡量--------------------------------#
	R=VSC_NodeData[:,10]
	Xl=VSC_NodeData[:,11]
	a=arctan(R/Xl)
	Y=1/np.sqrt(R*R+Xl*Xl)
	Usi=VSC_NodeData[:,5]
	VSC_Num = VSC_NodeData.shape[0] # 节点数目
	P = np.zeros([VSC_Num,1])
	Q = np.zeros([VSC_Num,1])
	Deltad1 = np.zeros([VSC_Num,1])
	Deltad2 = np.zeros([VSC_Num,1])
	Deltad3 = np.zeros([VSC_Num,1])
	Deltad4 = np.zeros([VSC_Num,1])
	iter = 0 
	for i in range(VSC_Num):  # 求解功率不平衡量
	    P[i] = (sqrt(6)/4)*N_VSC[i]*M[i]*Usi[i]*Udc[i]*Y[i]*sin(derta[i]-a[i]) + N_VSC[i]*Usi[i]*Usi[i]*Y[i]*sin(a[i])
	    Q[i] = -(sqrt(6)/4)*N_VSC[i]*M[i]*Usi[i]*Udc[i]*Y[i]*cos(derta[i]-a[i]) + N_VSC[i]*Usi[i]*Usi[i]*Y[i]*cos(a[i])  
	    Deltad1[iter] = Ps[i]- P[i] 
	    Deltad2[iter] = Qs[i]- Q[i] 
	    Deltad3[iter] = N_VSC[i]*Udc[i]*Idc[i]-0.5*sqrt(3/2)*N_VSC[i]*M[i]*Usi[i]*Udc[i]*Y[i]*sin(derta[i]+a[i]) + N_VSC[i]*Udc[i]*Udc[i]*(3/8*M[i]*M[i])*Y[i]*sin(a[i]) 
	    Deltad4[iter] = Idc[i]-N_VSC[i]*np.sum((Ydc[i+Nc,:]*Uc))
	    iter = iter+1
#------------------------------------整合不平衡量------------------------------------------#
	DeltaD = np.vstack([Delta_D1,Delta_D2,Delta_D3,Delta_D4,Delta_D5,Deltad1,Deltad2,Deltad3,Deltad4])  # 功率不平衡量
	Ng = 0
	for i in range(VSC_Num):
	    if VSC_NodeData[i,2]==3:
		    DeltaD = np.delete(DeltaD,5*Nc+i-Ng,0)   #去掉MMC平衡换流器Udc  
		    Ng = Ng+1
	# Option['string'] = '功率不平衡量为：\n'
	# Real(DeltaD,**Option)
	DC_MaxError = np.max(np.abs(DeltaD))
	if DC_MaxError<Tol:
	    print('直流偏差：',DC_MaxError)
	    return(Vd,Id,Kt,W,fi,Pdc,Qdc,Udc,Idc,derta,M,DC_MaxError,P,Q)
#---------------------------------------雅克比矩阵-----------------------------------------# 
	F11 = np.eye(Nc)
	F21 = np.eye(Nc)
	F31 = np.zeros([Nc,Nc])
	F32 = np.eye(Nc)
	for i in range(Nc):
	    if LCC_NodeData[i,2]==2:
	        F32[i,i]=F32[i,i]*(-1)
	F12 = np.diag(1.9*X)
	F13 = -np.diag(Vt.reshape(Nc)*W*2.7)
	F23 = -np.diag(kr*Vt.reshape(Nc)*cos(fi.reshape(Nc))*2.7)
	F14 = -np.diag(Kt*Vt.reshape(Nc)*2.7)
	F25 = np.diag(kr*Kt*Vt.reshape(Nc)*sin(fi.reshape(Nc))*2.7)
	F41 = np.zeros([Nc,Nc])
	F42 = np.zeros([Nc,Nc])
	F43 = np.zeros([Nc,Nc])
	F44 = np.zeros([Nc,Nc])
	F45 = np.zeros([Nc,Nc])
	F51 = np.zeros([Nc,Nc])
	F52 = np.zeros([Nc,Nc])
	F53 = np.zeros([Nc,Nc])
	F54 = np.zeros([Nc,Nc])
	F55 = np.zeros([Nc,Nc])
	A   = np.zeros([Nc,Nc])
	# F41 F42 F43 F44
	F4_iter = -1
	for i in range(Nc):
		F4_iter = F4_iter+1 # 记录行数
		F4_iter_y = -1 # 记录列数
		for j in range(4*Nc):
			F4_iter_y = F4_iter_y+1
			if j>=0 and j<Nc: # F41
				if j==i and control1[j]==2:
					F41[F4_iter,F4_iter_y]=1
				elif j==i and control1[j]==3:
					F41[F4_iter,F4_iter_y]=Id[i]
			elif j>=Nc and j<2*Nc: # F42
				if j-Nc==i and control1[j-Nc]==1:
					F42[F4_iter,F4_iter_y-Nc]=1
				elif j-Nc==i and control1[j-Nc]==3:
					F42[F4_iter,F4_iter_y-Nc]=Vd[i]
			elif j>=2*Nc and j<3*Nc: # F43
				if j-2*Nc==i and control1[j-2*Nc]==5:
					F43[F4_iter,F4_iter_y-2*Nc]=1
			elif j>=3*Nc and j<4*Nc: # F44
				if j-3*Nc==i and control1[j-3*Nc]==4:
					F44[F4_iter,F4_iter_y-3*Nc]=1
	# F51 F52 F53 F54
	F5_iter = -1
	for i in range(Nc):
		F5_iter = F5_iter+1 # 记录行数
		F5_iter_y = -1 # 记录列数
		for j in range(4*Nc):
			F5_iter_y = F5_iter_y+1
			if j>=0 and j<Nc: # F51
				if j==i and control2[j]==2:
					F51[F5_iter,F5_iter_y]=1
				elif j==i and control2[j]==3:
					F51[F5_iter,F5_iter_y]=Id[i]
			elif j>=Nc and j<2*Nc: # F52
				if j-Nc==i and control2[j-Nc]==1:
					F52[F5_iter,F5_iter_y-Nc]=1
				elif j-Nc==i and control2[j-Nc]==3:
					F52[F5_iter,F5_iter_y-Nc]=Vd[i]
			elif j>=2*Nc and j<3*Nc: # F53
				if j-2*Nc==i and control2[j-2*Nc]==5:
					F53[F5_iter,F5_iter_y-2*Nc]=1
			elif j>=3*Nc and j<4*Nc: # F54
				if j-3*Nc==i and control2[j-3*Nc]==4:
					F54[F5_iter,F5_iter_y-3*Nc]=1
	F1 = np.concatenate((F11,F12,F13,F14,A),axis=1)
	F2 = np.concatenate((F21,A,F23,A,F25),axis=1)
	F3 = np.concatenate((F31,F32,A,A,A),axis=1)
	F4 = np.concatenate((F41,F42,F43,F44,F45),axis=1)
	F5 = np.concatenate((F51,F52,F53,F54,F55),axis=1)
	F = np.concatenate((F1,F2,F3,F4,F5)) 
#-----------------------------------VSC------------------#
	D11=np.zeros([VSC_Num,VSC_Num])
	D12=np.zeros([VSC_Num,VSC_Num]) # 0阵 
	D13=np.zeros([VSC_Num,VSC_Num])
	D14=np.zeros([VSC_Num,VSC_Num])
	D21=np.zeros([VSC_Num,VSC_Num])
	D22=np.zeros([VSC_Num,VSC_Num]) # 0阵 
	D23=np.zeros([VSC_Num,VSC_Num])
	D24=np.zeros([VSC_Num,VSC_Num])
	D31=np.zeros([VSC_Num,VSC_Num])
	D32=np.zeros([VSC_Num,VSC_Num]) 
	D33=np.zeros([VSC_Num,VSC_Num])
	D34=np.zeros([VSC_Num,VSC_Num])
	D41=np.zeros([VSC_Num,VSC_Num])
	D42=np.eye(VSC_Num)   # 单位阵 
	D43=np.zeros([VSC_Num,VSC_Num]) # 0阵 
	D44=np.zeros([VSC_Num,VSC_Num]) # 0阵
	for i in range(VSC_Num):
		D11[i,i]=-(sqrt(6)/4)*N_VSC[i]*M[i]*Usi[i]*Y[i]*sin(derta[i]-a[i])
		D13[i,i]=-(sqrt(6)/4)*N_VSC[i]*M[i]*Usi[i]*Udc[i]*Y[i]*cos(derta[i]-a[i])
		D14[i,i]=-(sqrt(6)/4)*N_VSC[i]*Usi[i]*Udc[i]*Y[i]*sin(derta[i]-a[i])
		D21[i,i]=(sqrt(6)/4)*N_VSC[i]*M[i]*Usi[i]*Y[i]*cos(derta[i]-a[i])
		D23[i,i]=-(sqrt(6)/4)*N_VSC[i]*M[i]*Usi[i]*Udc[i]*Y[i]*sin(derta[i]-a[i])
		D24[i,i]=(sqrt(6)/4)*N_VSC[i]*Usi[i]*Udc[i]*Y[i]*cos(derta[i]-a[i])
		D31[i,i]=N_VSC[i]*Idc[i]-N_VSC[i]*(sqrt(6)/4)*M[i]*Usi[i]*Y[i]*sin(derta[i]+a[i])+N_VSC[i]*Udc[i]*(3/4)*M[i]*M[i]*Y[i]*sin(a[i]) 
		D32[i,i]=N_VSC[i]*Udc[i]
		D33[i,i]=-(sqrt(6)/4)*N_VSC[i]*M[i]*Usi[i]*Udc[i]*Y[i]*cos(derta[i]+a[i])
		D34[i,i]=-(sqrt(6)/4)*N_VSC[i]*Usi[i]*Udc[i]*Y[i]*sin(derta[i]+a[i])+N_VSC[i]*Udc[i]*Udc[i]*(3/4*M[i])*Y[i]*sin(a[i])
	D = np.vstack([np.hstack([D11,D12,D13,D14]),np.hstack([D21,D22,D23,D24]),np.hstack([D31,D32,D33,D34]),np.hstack([D41,D42,D43,D44])])
#---------------------------整合雅克比矩阵------------------------------#
	J = scipy.linalg.block_diag(F,D)
#--------------------------LCC与MMC直流网络方程偏导---------------------# 
	for i in range(Nc):
		for j in range(Nc):
		    J[2*Nc+i,j]=-N_LCC[j]*Ydc[i,j]                                    #注意：LCC标号靠前
		for j in range(VSC_Num):
		    J[2*Nc+i,5*Nc+j]=-N_VSC[j]*Ydc[i,Nc+j]
	for i in range(VSC_Num):
		for j in range(Nc):
		    J[5*Nc+3*VSC_Num+i,j]=-N_LCC[j]*Ydc[Nc+i,j]
		for j in range(VSC_Num):
		    J[5*Nc+3*VSC_Num+i,5*Nc+j]=-N_VSC[j]*Ydc[Nc+i,Nc+j]
#---------------------------去掉MMC平衡换流器P方程---------------------#                                                                                     # MMC平衡换流器个数
	NP = 0
	for i in range(VSC_Num):
		if VSC_NodeData[i,2]==3:
			J = np.delete(J,5*Nc+i-NP,0)
			J = np.delete(J,5*Nc+i-NP,1)
			NP = NP+1
#------------------------------修正------------------------------------#
	# Option['string'] = 'jacobi矩阵为：\n'
	# Real(J,**Option)
	Delta = np.linalg.solve(J,DeltaD)
	# Option['string'] = '方程组求解结果：\n'
	# Real(Delta,**Option)
	Vd = Vd-Delta[0:Nc].reshape(Nc)                                   # 对直流的影响
	Id = Id-Delta[Nc:2*Nc].reshape(Nc)
	Kt = Kt-Delta[2*Nc:3*Nc].reshape(Nc)
	W= W-Delta[3*Nc:4*Nc].reshape(Nc)
	fi = fi-Delta[4*Nc:5*Nc]
#---------------------------------------------------MMC修正-----------------------------------------#
	k=0
	for i in range(VSC_Num):
	    if VSC_NodeData[i,2]==3:
	        Udc[i] = Udc[i]
	    else:
	        Udc[i] = Udc[i]-Delta[5*Nc+k]
	        k=k+1 
	NV=VSC_Num-NP                  
	Idc = Idc-Delta[5*Nc+NV:5*Nc+NV+VSC_Num].reshape(VSC_Num)
	derta = derta-Delta[5*Nc+NV+VSC_Num:5*Nc+NV+2*VSC_Num].reshape(VSC_Num)
	M = M-Delta[5*Nc+NV+2*VSC_Num:5*Nc+NV+3*VSC_Num].reshape(VSC_Num)
	return(Vd,Id,Kt,W,fi,Pdc,Qdc,Udc,Idc,derta,M,DC_MaxError,P,Q)
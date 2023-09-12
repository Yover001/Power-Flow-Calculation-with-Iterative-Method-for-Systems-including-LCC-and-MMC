# 牛顿拉夫逊求解
import numpy as np
from OutTxt import Real
def PolarNR(U,Angle,Y,PQNode,PVNode,SlackNode,P_Real,Q_Real,Tol,**Option):
    P_iter = 0  # 为形成雅可比矩阵
    Q_iter = 0  # 为形成雅可比矩阵
    PQNode = PQNode-1
    PVNode = PVNode-1
    SlackNode = SlackNode-1
    NumNode = Y.shape[0] # 节点数目
    NumPQ = max(PQNode.shape) # PQ节点数目
    G = Y.real
    B = Y.imag
    P = np.zeros([NumNode,1])
    Q = np.zeros([NumNode,1])
    DeltaP = np.zeros([NumNode-2,1])
    DeltaQ = np.zeros([NumPQ,1])
    for i in range(NumNode):  # 求解功率不平衡量
        P[i] = U[i]*np.sum(U*(G[i,:]*np.cos(Angle[i]-Angle) +  B[i,:]*np.sin(Angle[i]-Angle)))
        Q[i] = U[i]*np.sum(U*(G[i,:]*np.sin(Angle[i]-Angle) -  B[i,:]*np.cos(Angle[i]-Angle)))
        if i not in SlackNode:    # 不是平衡节点
            DeltaP[P_iter] = P_Real[i]-P[i]  # NumPQ+NumPV
            if i in PQNode:    # PQ节点
                DeltaQ[Q_iter] = Q_Real[i]-Q[i] # NumPQ
                Q_iter = Q_iter+1
            P_iter = P_iter+1
    DeltaPQ = np.vstack([DeltaP,DeltaQ])  # 功率不平衡量
    # Option['string'] = '功率不平衡量为：\n'
    # Real(DeltaPQ,**Option)
    MaxError = np.max(np.abs(DeltaPQ))
    if MaxError<Tol:
        print('交流偏差：',MaxError)
        return(U,Angle,MaxError,P,Q)
    ## H and N and J and L
    SlackNode_1level = [int(j) for i in SlackNode for j in i]
    PVNode_1level = [int(j) for i in PVNode for j in i]
    # H
    H = 2 * np.ones([NumNode,NumNode])
    np.fill_diagonal(H,1)
    H = np.delete(H,SlackNode_1level,axis=0)
    H = np.delete(H,SlackNode_1level,axis=1)
    U_delete = np.delete(U,SlackNode_1level)
    Angle_delete = np.delete(Angle,SlackNode_1level)
    Q_delete = np.delete(Q,SlackNode_1level)
    G_delete = np.delete(G,SlackNode_1level,axis=1) #axis=1代表列；axis=0代表行
    G_delete = np.delete(G_delete,SlackNode_1level,axis=0)
    B_delete = np.delete(B,SlackNode_1level,axis=1)
    B_delete = np.delete(B_delete,SlackNode_1level,axis=0)
    np.copyto(H, G_delete, where=((G_delete == 0) & (B_delete == 0) & (H != 1)))
    coordinate = np.where(H != 0)
    coordinate_1 = np.where(H == 1)
    H_ones_coordinate = list(zip(coordinate[0],coordinate[1]))
    H_ones_coordinate_1 = list(zip(coordinate_1[0],coordinate_1[1]))
    for element in H_ones_coordinate:
        if element in H_ones_coordinate_1:
            i = element[0]
            j = element[1]
            Angleij = Angle_delete[i]-Angle_delete[j]
            Anglei_j = Angle_delete[i]-Angle_delete
            # H[i,j] = U_delete[i]*(np.dot(U_delete*G_delete[i],np.sin(Anglei_j))-np.dot(U_delete*B_delete[i],np.cos(Anglei_j)))-U_delete[i]*U_delete[j]*(G_delete[i,j]*np.sin(Angleij)-B_delete[i,j]*np.cos(Angleij))
            H[i,j] = Q_delete[i]+U_delete[i]**2*B_delete[i,j]
        else:
            i = element[0]
            j = element[1]
            Angleij = Angle_delete[i]-Angle_delete[j]
            H[i,j] = -U_delete[i]*U_delete[j]*(G_delete[i,j]*np.sin(Angleij)-B_delete[i,j]*np.cos(Angleij))
    # N
    N = 2 * np.ones([NumNode,NumNode])
    np.fill_diagonal(N,1)
    N = np.delete(N,SlackNode_1level,axis=0)
    N = np.delete(N,SlackNode_1level+PVNode_1level,axis=1)
    U_delete_i = np.delete(U,SlackNode_1level)
    U_delete_j = np.delete(U,SlackNode_1level+PVNode_1level)
    Angle_delete_i = np.delete(Angle,SlackNode_1level)
    Angle_delete_j = np.delete(Angle,SlackNode_1level+PVNode_1level)
    P_delete = np.delete(P,SlackNode_1level)
    G_delete = np.delete(G,SlackNode_1level,axis=0)
    G_delete = np.delete(G_delete,SlackNode_1level+PVNode_1level,axis=1)
    B_delete = np.delete(B,SlackNode_1level,axis=0)
    B_delete = np.delete(B_delete,SlackNode_1level+PVNode_1level,axis=1)
    np.copyto(N, G_delete, where=((G_delete == 0) & (B_delete == 0) & (N != 1)))
    coordinate = np.where(N != 0)
    coordinate_1 = np.where(N == 1)
    N_ones_coordinate = list(zip(coordinate[0],coordinate[1]))
    N_ones_coordinate_1 = list(zip(coordinate_1[0],coordinate_1[1]))
    for element in N_ones_coordinate:
        if element in N_ones_coordinate_1:
            i = element[0]
            j = element[1]
            Angleij = Angle_delete[i]-Angle_delete[j]
            Anglei_j = Angle_delete[i]-Angle_delete_j
            # N[i,j] = -U_delete_i[i]*(np.dot(U_delete_j*G_delete[i],np.cos(Anglei_j))+np.dot(U_delete_j*B_delete[i],np.sin(Anglei_j)))+U_delete_i[i]*U_delete_j[j]*(G_delete[i,j]*np.cos(Angleij)+B_delete[i,j]*np.sin(Angleij))-2*G_delete[i,i]*U_delete_i[i]**2
            N[i,j] = -P_delete[i]-G_delete[i,j]*U_delete_i[i]**2
        else:
            i = element[0]
            j = element[1]
            Angleij = Angle_delete_i[i]-Angle_delete_j[j]
            N[i,j] = -U_delete_i[i]*U_delete_j[j]*(G_delete[i,j]*np.cos(Angleij)+B_delete[i,j]*np.sin(Angleij))
    # J
    J = 2 * np.ones([NumNode,NumNode])
    np.fill_diagonal(J,1)
    J = np.delete(J,SlackNode_1level+PVNode_1level,axis=0)
    J = np.delete(J,SlackNode_1level,axis=1)
    U_delete_i = np.delete(U,SlackNode_1level+PVNode_1level)
    U_delete_j = np.delete(U,SlackNode_1level)
    Angle_delete_i = np.delete(Angle,SlackNode_1level+PVNode_1level)
    Angle_delete_j = np.delete(Angle,SlackNode_1level)
    P_delete = np.delete(P,SlackNode_1level+PVNode_1level)
    G_delete = np.delete(G,SlackNode_1level+PVNode_1level,axis=0)
    G_delete = np.delete(G_delete,SlackNode_1level,axis=1)
    B_delete = np.delete(B,SlackNode_1level+PVNode_1level,axis=0)
    B_delete = np.delete(B_delete,SlackNode_1level,axis=1)
    np.copyto(J, G_delete, where=((G_delete == 0) & (B_delete == 0) & (J != 1)))
    coordinate = np.where(J != 0)
    coordinate_1 = np.where(J == 1)
    J_ones_coordinate = list(zip(coordinate[0],coordinate[1]))
    J_ones_coordinate_1 = list(zip(coordinate_1[0],coordinate_1[1]))
    for element in J_ones_coordinate:
        if element in J_ones_coordinate_1:
            i = element[0]
            j = element[1]
            Angleij = Angle_delete[i]-Angle_delete[j]
            Anglei_j = Angle_delete[i]-Angle_delete_j
            # J[i,j] = -U_delete_i[i]*(np.dot(U_delete_j*G_delete[i],np.cos(Anglei_j))+np.dot(U_delete_j*B_delete[i],np.sin(Anglei_j)))+U_delete_i[i]*U_delete_j[j]*(G_delete[i,j]*np.cos(Angleij)+B_delete[i,j]*np.sin(Angleij))
            J[i,j] = -P_delete[i]+G_delete[i,j]*U_delete_i[i]**2
        else:
            i = element[0]
            j = element[1]
            Angleij = Angle_delete_i[i]-Angle_delete_j[j]
            J[i,j] = U_delete_i[i]*U_delete_j[j]*(G_delete[i,j]*np.cos(Angleij)+B_delete[i,j]*np.sin(Angleij))
    # L
    L = 2 * np.ones([NumNode,NumNode])
    np.fill_diagonal(L,1)
    L = np.delete(L,SlackNode_1level+PVNode_1level,axis=0)
    L = np.delete(L,SlackNode_1level+PVNode_1level,axis=1)
    U_delete = np.delete(U,SlackNode_1level+PVNode_1level)
    Angle_delete = np.delete(Angle,SlackNode_1level+PVNode_1level)
    Q_delete = np.delete(Q,SlackNode_1level+PVNode_1level)
    G_delete = np.delete(G,SlackNode_1level+PVNode_1level,axis=0)
    G_delete = np.delete(G_delete,SlackNode_1level+PVNode_1level,axis=1)
    B_delete = np.delete(B,SlackNode_1level+PVNode_1level,axis=0)
    B_delete = np.delete(B_delete,SlackNode_1level+PVNode_1level,axis=1)
    np.copyto(L, G_delete, where=((G_delete == 0) & (B_delete == 0) & (L != 1)))
    coordinate = np.where(L != 0)
    coordinate_1 = np.where(L == 1)
    L_ones_coordinate = list(zip(coordinate[0],coordinate[1]))
    L_ones_coordinate_1 = list(zip(coordinate_1[0],coordinate_1[1]))
    for element in L_ones_coordinate:
        if element in L_ones_coordinate_1:
            i = element[0]
            j = element[1]
            Angleij = Angle_delete[i]-Angle_delete[j]
            Anglei_j = Angle_delete[i]-Angle_delete
            # L[i,j] = -U_delete[i]*(np.dot(U_delete*G_delete[i],np.sin(Anglei_j))-np.dot(U_delete*B_delete[i],np.cos(Anglei_j)))+U_delete[i]*U_delete[j]*(G_delete[i,j]*np.sin(Angleij)-B_delete[i,j]*np.cos(Angleij))+2*B_delete[i,i]*U_delete[i]**2
            L[i,j] = -Q_delete[i]+B_delete[i,j]*U_delete[i]**2
        else:
            i = element[0]
            j = element[1]
            Angleij = Angle_delete[i]-Angle_delete[j]
            L[i,j] = -U_delete[i]*U_delete[j]*(G_delete[i,j]*np.sin(Angleij)-B_delete[i,j]*np.cos(Angleij))
    # 修正
    Jaccobi = np.vstack([np.hstack([H,N]),np.hstack([J,L])])
    # Option['string'] = 'jacobi矩阵为：\n'
    # Real(Jaccobi,**Option)
    Delta = np.linalg.solve(Jaccobi,DeltaPQ)
    # Option['string'] = '方程组求解结果：\n'
    # Real(Delta,**Option)
    DeltaAngle = Delta[0:NumNode-2]  # 注意下标
    DeltaU_U = Delta[NumNode-2:]
    DA_iter = -1
    U_U_iter = -1
    for i in range(NumNode):
        if i not in SlackNode:
            DA_iter = DA_iter+1
            Angle[i] = Angle[i]-DeltaAngle[DA_iter]
            if i in PQNode:
                U_U_iter = U_U_iter+1
                U[i] = U[i]-U[i]*DeltaU_U[U_U_iter] #不理解
    return(U,Angle,MaxError,P,Q)

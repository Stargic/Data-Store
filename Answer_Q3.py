import numpy as np
import math
from scipy.optimize import minimize
from time import time

#var
index = 624
MT = [0]*index
# MT[0] ->seed

def inter(t):
    return(0xFFFFFFFF & t) #取最后32位->t

def twister():
    global index
    for i in range(624):
        y = inter((MT[i] & 0x80000000) +(MT[(i + 1) % 624] & 0x7fffffff))
        MT[i] = MT[(i + 397) % 624] ^ y >> 1
        if y % 2 != 0:
            MT[i] = MT[i] ^ 0x9908b0df
    index = 0

def exnum():
    global index
    if index >= 624:
        twister()
    y = MT[index]
    y = y ^ y >> 11
    y = y ^ y << 7 & 2636928640
    y = y ^ y << 15 & 4022730752
    y = y ^ y >> 18
    index = index + 1
    return inter(y)

def mainset(seed):
    MT[0] = seed    #seed
    for i in range(1,624):
        MT[i] = inter(1812433253 * (MT[i - 1] ^ MT[i - 1] >> 30) + i)
    return exnum()

def GetRandom(minnum,maxnum):
    return minnum + int((maxnum-minnum)*(mainset(int(time())) / (2**32-1)))

# 产生随机的误差数据
def RandomData():
    F0 = [13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13]
    t = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    count = GetRandom(0,17)
    for i in range(count):
        F0[GetRandom(1,17)] = GetRandom(10,30)

    for i in range(count):
        if GetRandom(1,17)%2 == 0:
            t[GetRandom(1,17)] = -0.1

    return F0,t

# 获得鼓与水平面可能产生的倾角
def GetQ(F0,t,hs):
    
    n = 17
    PI = 3.1415926
    g = 9.8   # 重力加速度

    # 鼓的参数
    r = 0.2  # 鼓的半径 
    M = 3.6  # 鼓的质量
    l = 1.7  # 绳长
    J = 0.05455 # 转动惯量（自定义的量） 0.0885
    
    # 获取每个人站位的单位向量
    e = []
    for i in range(1,n+1):
        e.append(np.mat([math.cos(i*2*PI/n),math.sin(i*2*PI/n),0]))

    # 角度的初始值及变化后的值
    Q = Q_t = np.mat([0,0,0])

    # 初始时刻的转动惯量
    wt0 = wt0_t = 0

    # 初始状态下鼓倾斜的上升距离
    s = s_t = 0

    # 初始状态下鼓倾斜时的上升速度
    v = v_t = 0
    
    # 表示人手的位置
    Hand = []
    for k in range(0,n):
        Hand.append(math.sqrt(l * l + hs * hs) * e[k] + np.mat([0,0,hs]))

    FList = []
    # 初始状态下所有人拉住鼓的基础拉力的列表
    for i in range(0,n):
        FList.append(13)

    # 时间片划分的标注，在循环中，值为100时表示将0~0.1的时间进行等间隔的100次剖分
    TimeSlice = 100
    for i in range(0,n):
        # 若有人提前发力，则剖分的时间片为 0~0.2
        if t[i] != 0:
            TimeSlice = 200
            # 若有人提前发力则将对应人的拉力更改为他发力时的拉力
            FList[i] = F0[i]

    # 若有人提前发力（TimeSlice != 100）则将更新后的基础拉力列表提供给变量 F
    # 若无人提前发力，则将同时发力时发力列表提供给变量 F
    if TimeSlice != 100:
        F = FList
    else:
        F = F0

    # 划分 0 ~ 0.1 时间为100份
    # 若有人提前施加拉力，则将施加拉力的时刻定为0，并将总时间 0~0.2 部分划分200份
    dt = 0.001

    for t0 in range(0, TimeSlice):
        # 若TimeSlice == 100时所有人开始发力，此时则将拉力情况更换为F0
        if TimeSlice == 100:
            F = F0
        # 绳到质心的矢量
        sl = []
        for k in range(0,n):
            sl.append(r * np.cross(Q,e[k]) + r * e[k] + np.mat([0,0,np.linalg.norm(s)]))

        Ft = np.mat([0,0,0])
        for k in range(0,n):
            Ft = Ft + ((Hand[k] - sl[k]) / np.linalg.norm(Hand[k] - sl[k])) * F[k]

        Ft = Ft + M * g * np.mat([0,0,-1])
        v_t = dt * Ft / M + v
        s_t = s + v * dt

        Mt = np.mat([0,0,0])
        for k in range(0,n):
            Mt = Mt + np.cross((r * np.cross(Q,e[k]) + r * e[k]) , ((Hand[k] - sl[k]) / np.linalg.norm(Hand[k] - sl[k])) * F[k])

        wt0_t = wt0 + dt*Mt/J
        Q_t = Q + dt*wt0

        # 下一时刻将计算出来的值作为已知量
        s = s_t
        v = v_t
        Q = Q_t
        wt0 = wt0_t
    
    return np.linalg.norm(Q) * 180 / PI

def Func(hs):
    # 随机获得带有误差的数据
    F0,t = RandomData()
    print("F0 = ",F0)
    print("t = ",t)
    return F0,t,GetQ(F0,t,hs)

def TestQ(args):
    F0,t = args
    Q = lambda hs:GetQ(F0,t,hs)
    return Q

file_name = 'Question_3_data.txt' 
with open(file_name,'w') as file_obj:
    for i in range(50,51):
        print("第",i+1,"组测试数据——————————————————————————————\n")
        file_obj.write("第"+str(i+1)+"组测试数据——————————————————————————————\n")
        F0,t,Q = Func(0.11)
        file_obj.write(str(F0)+"\n")
        print("优化前因误差产生的倾斜角度为："+str(Q))
        file_obj.write("优化前因误差产生的倾斜角度为："+str(Q))

        args = (F0,t)

        print("优化开始！")
        file_obj.write("优化开始！\n")
        v = TestQ
        x0 = np.asarray((10)) 
        cons = ({'type': 'ineq', 'fun': lambda s: 1 - s},{'type': 'ineq', 'fun': lambda s: s})
        res = minimize(v(args), x0, method='SLSQP',constraints=cons)
        print("是否完成优化：",res.success)
        if res.success == True:
            file_obj.write("是否完成优化：True\n")
        else:
            file_obj.write("是否完成优化：False\n")
        print("改变策略的方法：令鼓的初始位置升高到",res.x)
        file_obj.write("改变策略的方法：令鼓的初始位置升高到"+str(res.x)+"\n")
        Q = GetQ(F0,t,float(res.x))
        print("优化后会产生的倾斜角度为：",Q)
        file_obj.write("优化后会产生的倾斜角度为："+str(Q)+"\n")
        print("—————————————————————————————————————————\n")
        file_obj.write("—————————————————————————————————————————\n")
        


import numpy as np
import math

# F0 = 所有人同时发力时的发力列表
# t = 每个人的发力时机
# n = 参加的人数
def GetAnswer_2(F0,t,n):
    # 默认值
    PI = 3.1415926
    g = 9.8   # 重力加速度

    # 鼓的参数
    r = 0.2  # 鼓的半径 
    M = 3.6  # 鼓的质量
    l = 1.7  # 绳长
    J = 0.05455 # 转动惯量（自定义的量） 0.0885

    # BF为基础拉力 0.11 为鼓身高度的半 
    BF = l * M *10 / (0.11*n)
    
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
        Hand.append((r + l) * e[k] + np.mat([0,0,0.11]))

    FList = []
    # 初始状态下所有人拉住鼓的基础拉力的列表
    for i in range(0,n):
        FList.append(BF)

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


# 设置表1中的数据
InputList = [[90,80,80,80,80,80,80,80],[0,0,0,0,0,0,0,0],
             [90,90,80,80,80,80,80,80],[0,0,0,0,0,0,0,0],
             [90,80,80,90,80,80,80,80],[0,0,0,0,0,0,0,0],
             [80,80,80,80,80,80,80,80],[-0,1,0,0,0,0,0,0,0],
             [80,80,80,80,80,80,80,80],[-0.1,-0.1,0,0,0,0,0,0],
             [80,80,80,80,80,80,80,80],[-0.1,0,0,-0.1,0,0,0,0],
             [90,80,80,80,80,80,80,80],[-0.1,0,0,0,0,0,0,0],
             [90,80,80,90,80,80,80,80],[0,-0.1,0,0,-0.1,0,0,0],
             [90,80,80,90,80,80,80,80],[0,0,0,0,-0.1,0,0,-0.1]]

# 循环读入表1中的数据并输出结果
for i in range(0,18,2):
    print(GetAnswer_2(InputList[i],InputList[i + 1],8))
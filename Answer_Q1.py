import math
import numpy as np
PI = 3.1415926
# 函数声明
def GetMaxNum_of_People():
    n = 8    # 至少8人
    List1 = []
    List2 = []
    while True:
        # np.arccos方法在遇到无法求解的情况时值为NaN
        a = np.arccos(0.056 * n)
        a = a * (180/PI) 
        List1.append([a,n])
        
        try:
            # math.acos方法在遇到无法求解的情况时会报错
            b = math.acos(0.056 * n)  
            b = b * (180/PI)  
            List2.append([b,n])
        except:
            pass
        n += 1
        if n > 20:
            break
    print(List1)
    print(List2)

# 运行程序
GetMaxNum_of_People()

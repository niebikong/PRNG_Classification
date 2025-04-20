import math
import scipy.special as ss

# RANDOM EXCURSION TEST
def random_excursion_test(input, n):
    """
    执行随机偏移测试，用于评估输入序列的随机性。

    :param input: 输入的二进制字符串。
    :param n: 输入序列的长度。

    :return: 返回包含测试结果的列表，包括以下信息：
        - n: 输入序列的长度。
        - J: 周期数。
        - chi_sq: 卡方统计量列表。
        - plist: p 值列表。
        - p_total: 平均 p 值。
        - success: 是否通过检验的布尔值。
    """

    # Convert to +1,-1
    x = list()
    for i in range(n):
        x.append(int(input[i],2)*2 -1 )

    # Build the partial sums
    pos = 0
    s = list()
    for e in x:
        pos = pos+e
        s.append(pos)    
    sprime = [0]+s+[0] # Add 0 on each end
    

    # Build the list of cycles
    pos = 1
    cycles = list()
    while (pos < len(sprime)):
        cycle = list()
        cycle.append(0)
        while sprime[pos]!=0:
            cycle.append(sprime[pos])
            pos += 1
        cycle.append(0)
        cycles.append(cycle)
        pos = pos + 1
    
    J = len(cycles)  
    
    vxk = [['a','b','c','d','e','f'] for y in [-4,-3,-2,-1,1,2,3,4] ]

    # Count Occurances  
    for k in range(6):
        for index in range(8):
            mapping = [-4,-3,-2,-1,1,2,3,4]
            x = mapping[index]
            cyclecount = 0
            #count how many cycles in which x occurs k times
            for cycle in cycles:
                oc = 0
                #Count how many times x occurs in the current cycle
                for pos in cycle:
                    if (pos == x):
                        oc += 1
                # If x occurs k times, increment the cycle count
                if (k < 5):
                    if oc == k:
                        cyclecount += 1
                else:
                    if k == 5:
                        if oc >=5:
                            cyclecount += 1
            vxk[index][k] = cyclecount
    
    # Table for reference random probabilities 
    pikx=[[0.5     ,0.25   ,0.125  ,0.0625  ,0.0312 ,0.0312],
          [0.75    ,0.0625 ,0.0469 ,0.0352  ,0.0264 ,0.0791],
          [0.8333  ,0.0278 ,0.0231 ,0.0193  ,0.0161 ,0.0804],
          [0.875   ,0.0156 ,0.0137 ,0.012   ,0.0105 ,0.0733],
          [0.9     ,0.01   ,0.009  ,0.0081  ,0.0073 ,0.0656],
          [0.9167  ,0.0069 ,0.0064 ,0.0058  ,0.0053 ,0.0588],
          [0.9286  ,0.0051 ,0.0047 ,0.0044  ,0.0041 ,0.0531]]
    
    success = True
    plist = list()
    chi_sq = list()
    p_total = 0.0
    for index in range(8):
        #list of states
        mapping = [-4,-3,-2,-1,1,2,3,4] 
        x = mapping[index]
        chisq = 0.0
        for k in range(6):
            top = float(vxk[index][k]) - (float(J) * (pikx[abs(x)-1][k]))
            top = top*top
            bottom = J * pikx[abs(x)-1][k]
            chisq += top/bottom

        p = ss.gammaincc(5.0/2.0,chisq/2.0)
        p_total += p
        plist.append(p)
        chi_sq.append(chisq)
        if p < 0.01:
            success = False

    return [n, J, chi_sq, plist, p_total/8, success]



if __name__ == "__main__":
    # 定义输入序列和参数
    input_sequence = "1100101010101010110001111000010101010110101010101010111101010101010101010101010101010101010101"
    sequence_length = len(input_sequence)

    # 执行测试
    result = random_excursion_test(input_sequence, sequence_length)

    # 输出测试结果
    print("输入序列的长度:", result[0])
    print("周期数:", result[1])
    print("卡方统计量列表:", result[2])
    print("p 值列表:", result[3])
    print("平均 p 值:", result[4])
    print("是否通过检验:", result[5])
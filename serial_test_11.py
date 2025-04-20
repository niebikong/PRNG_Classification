from __future__ import print_function
import math
import scipy.special as ss

def int2patt(n,m):
    """
    将整数转换为模式列表。

    :param n: 整数表示的模式。
    :param m: 模式长度。

    :return: 返回模式的列表表示。
    """
    pattern = list()
    for i in range(m):
        pattern.append((n >> i) & 1)
    return pattern
    
def countpattern(patt,padded_input,n):
    """
    计算给定模式在填充输入中的出现次数。

    :param patt: 要匹配的模式。
    :param padded_input: 填充后的输入序列。
    :param n: 输入序列的长度。

    :return: 返回模式在输入中的出现次数。
    """
    thecount = 0
    for i in range(n):
        match = True
        for j in range(len(patt)):
            if str(patt[j]) != padded_input[i+j]:
                match = False
                break
        if match:
            thecount += 1
    return thecount

def psi_sq_mv1(m, n, padded_input):
    """
    计算 psi-square 统计量的值。

    :param m: 模式长度。
    :param n: 输入序列的长度。
    :param padded_input: 填充后的输入序列。

    :return: 返回 psi-square 统计量的值。
    """
    counts = [0 for i in range(2**m)]
    for i in range(2**m):

        pattern = int2patt(i,m)
        count = countpattern(pattern,padded_input,n)
        counts.append(count)

        # pattern = padding(bin(i)[2:], m)


        # count = padded_input.count(pattern)
        
        # counts.append(count)
        
    psi_sq_m = 0.0
    for count in counts: 
        psi_sq_m += (count**2)
    psi_sq_m = psi_sq_m * (2**m)/n 
    psi_sq_m -= n
    return psi_sq_m            
         
def serial_test(input, n, patternlen=None):
    """
    执行统计测试，用于评估输入序列的随机性。

    :param input: 输入的二进制字符串。
    :param n: 输入序列的长度。
    :param patternlen: 模板长度，可选参数，默认为 None。如果不提供，则根据输入长度选择默认值。

    :return: 返回包含测试结果的列表，包括以下信息：
        - psi_sq_m: 当前模板长度的 psi-square 统计量。
        - psi_sq_mm1: 模板长度减一的 psi-square 统计量。
        - psi_sq_mm2: 模板长度减二的 psi-square 统计量。
        - delta1: 当前模板长度与模板长度减一的差异。
        - delta2: 当前模板长度与模板长度减二的差异。
        - p1: 模板长度与模板长度减一之间差异的 p 值。
        - p2: 模板长度与模板长度减二之间差异的 p 值。
        - average_p: p1 和 p2 的平均值。
        - success: 是否通过检验的布尔值。
    """
    # Pattern length
    if patternlen != None:
        m = patternlen  
    else:  
        m = int(math.floor(math.log(n,2)))-2
    
        if m < 4:
            print("Error. Not enough data for m to be 4")
            return [0]*8
        m = 4

    # Step 1
    padded_input=input[0:n]+input[0:m-1]
    
    # Step 2
    psi_sq_m   = psi_sq_mv1(m, n, padded_input)
    psi_sq_mm1 = psi_sq_mv1(m-1, n, padded_input)
    psi_sq_mm2 = psi_sq_mv1(m-2, n, padded_input)    
    
    delta1 = psi_sq_m - psi_sq_mm1
    delta2 = psi_sq_m - (2*psi_sq_mm1) + psi_sq_mm2

    p1 = ss.gammaincc(2**(m-2),delta1/2.0)
    p2 = ss.gammaincc(2**(m-3),delta2/2.0)
     
    success = (p1 >= 0.01) and (p2 >= 0.01)

    return [psi_sq_m, psi_sq_mm1, psi_sq_mm2, delta1, delta2, p1, p2, (p1+p2)/2, success]

if __name__ == "__main__":
    # 定义输入序列和参数
    input_sequence = "1100101010101010110001111000010101010110101010101010111101010101010101010101010101010101010101"
    sequence_length = len(input_sequence)

    # 执行测试
    result = serial_test(input_sequence, sequence_length)

    # 输出测试结果
    print("当前模板长度的 psi-square 统计量:", result[0])
    print("模板长度减一的 psi-square 统计量:", result[1])
    print("模板长度减二的 psi-square 统计量:", result[2])
    print("当前模板长度与模板长度减一的差异:", result[3])
    print("当前模板长度与模板长度减二的差异:", result[4])
    print("模板长度与模板长度减一之间差异的 p 值:", result[5])
    print("模板长度与模板长度减二之间差异的 p 值:", result[6])
    print("p1 和 p2 的平均值:", result[7])
    print("是否通过检验:", result[8])
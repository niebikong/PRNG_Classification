import math
import scipy.special as ss

def normcdf(n):
    """
    计算标准正态分布函数的累积分布函数值。

    :param n: 输入值。

    :return: 标准正态分布函数的累积分布函数值。
    """
    return 0.5 * math.erfc(-n * math.sqrt(0.5))

def p_value(n,z):
    """
    计算给定统计值 z 的 p 值。

    :param n: 样本大小。
    :param z: 统计值。

    :return: p 值。
    """
    sum_a = 0.0
    startk = int(math.floor((((float(-n)/z)+1.0)/4.0)))
    endk   = int(math.floor((((float(n)/z)-1.0)/4.0)))
    for k in range(startk,endk+1):
        c = (((4.0*k)+1.0)*z)/math.sqrt(n)
        #d = scipy.stats.norm.cdf(c)
        d = normcdf(c)
        c = (((4.0*k)-1.0)*z)/math.sqrt(n)
        #e = scipy.stats.norm.cdf(c)
        e = normcdf(c)
        sum_a = sum_a + d - e

    sum_b = 0.0
    startk = int(math.floor((((float(-n)/z)-3.0)/4.0)))
    endk   = int(math.floor((((float(n)/z)-1.0)/4.0)))
    for k in range(startk,endk+1):
        c = (((4.0*k)+3.0)*z)/math.sqrt(n)
        #d = scipy.stats.norm.cdf(c)
        d = normcdf(c)
        c = (((4.0*k)+1.0)*z)/math.sqrt(n)
        #e = scipy.stats.norm.cdf(c)
        e = normcdf(c)
        sum_b = sum_b + d - e 

    p = 1.0 - sum_a + sum_b
    return p
    
def cumulative_sums_test(input, n):
    """
    执行统计测试，用于评估输入序列的随机性。

    :param input: 输入的二进制字符串。
    :param n: 输入序列的长度。

    :return: 返回包含测试结果的列表，包括以下信息：
        - p_forward: 正向最大偏移的 p 值。
        - p_backward: 反向最大偏移的 p 值。
        - success: 是否通过检验的布尔值。
    """
    # Step 1
    x = list()             # Convert to +1,-1
    for i in range(n):
        #if bit == 0:
        x.append(int(input[i],2)*2-1)
        
    # Steps 2 and 3 Combined
    # Compute the partial sum and records the largest excursion.
    pos = 0
    forward_max = 0
    for e in x:
        pos = pos+e
        if abs(pos) > forward_max:
            forward_max = abs(pos)
    pos = 0
    backward_max = 0
    for e in reversed(x):
        pos = pos+e
        if abs(pos) > backward_max:
            backward_max = abs(pos)
     
    # Step 4
    p_forward  = p_value(n, forward_max)
    p_backward = p_value(n,backward_max)
    
    success = ((p_forward >= 0.01) and (p_backward >= 0.01))

    return [p_forward, p_backward, success]

if __name__ == "__main__":
    # 定义输入序列和参数
    input_sequence = "1100101010101010110001111000010101010110101010101010111101010101010101010101010101010101010101"
    sequence_length = len(input_sequence)

    # 执行测试
    result = cumulative_sums_test(input_sequence, sequence_length)

    # 输出测试结果
    print("正向最大偏移的 p 值:", result[0])
    print("反向最大偏移的 p 值:", result[1])
    print("是否通过检验:", result[2])
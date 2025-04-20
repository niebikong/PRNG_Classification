import scipy.special as ss

def longest_run_ones_in_a_block_test(input, n):
    """
    执行最长重复序列检测（Longest Run Test），用于检测输入序列中连续重复比特的长度。

    :param input: 输入的序列，应为字符串，包含 '0' 和 '1'。
    :param n: 输入序列的长度。

    :return: 返回一个包含检验结果的列表，包括以下信息：
        - 卡方统计量
        - p 值
        - 是否通过检验（布尔值）
    """

    M8 = [0.2148, 0.3672, 0.2305, 0.1875]

    # Length of blocks
    M = 8 
                
    K = 3

    N = 16
            
    # Table of frequencies
    v = [0,0,0,0,0,0,0]

    for i in range(N): # over each block
        #find the longest run
        if (i+1)*M < n:  # 不越界
            block = input[i*M:((i+1)*M)] # Block i
        
        run = 0
        longest = 0
        for j in range(M): # Count the bits.
            if block[j] == '1':
                run += 1
                if run > longest:
                    longest = run
            else:
                run = 0

        if longest <= 1:    v[0] += 1
        elif longest == 2:  v[1] += 1
        elif longest == 3:  v[2] += 1
        else:               v[3] += 1
    
    # Compute Chi-Sq
    chi_sq = 0.0
    for i in range(K+1):
        p_i = M8[i]
        upper = (v[i] - N*p_i)**2
        lower = N*p_i
        chi_sq += upper/lower
    # p-value
    p = ss.gammaincc(K/2.0, chi_sq/2.0)
    
    success = (p>=0.01)

    return [chi_sq, p, success]

if __name__ == "__main__":
    # 定义输入序列和长度
    input_sequence = "1100101010101010110001111000010101010110101010101010111101010101010101010101010101010101010101"
    sequence_length = len(input_sequence)

    # 执行最长重复序列检测测试
    result = longest_run_ones_in_a_block_test(input_sequence, sequence_length)

    # 输出测试结果
    print("卡方统计量:", result[0])
    print("p 值:", result[1])
    print("是否通过检验:", result[2])
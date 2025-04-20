import math
            
def maurers_universal_test(input, n, patternlen=None, initblocks=None):
    """
    执行重叠计数测试（Overlapping Template Matching Test），用于检测输入序列中是否存在与模板匹配的模式。

    :param input: 输入的序列，应为一个二进制字符串，例如 '00101010111010101101'。
    :param n: 输入序列的总长度。
    :param patternlen: 模板的长度，可选参数，默认为 None。如果不提供，则根据输入序列长度自动选择合适的模板长度。
    :param initblocks: 初始块数量，可选参数，默认为 None。如果不提供，则根据输入序列长度自动选择初始块数量。

    :return: 返回一个包含测试结果的列表，包括以下信息：
        - 输入序列长度（n）
        - 块的数量（nblocks）
        - 模板的长度（L）
        - 自由度的数量（K）
        - 初始块数量（Q）
        - 标准差（sigma）
        - p 值
        - 是否通过检验（布尔值）
    """
    # Step 1. Choose the block size
    if patternlen != None:
        L = patternlen  
    else: 
        ns = [904960,2068480,4654080,10342400,
              22753280,49643520,107560960,
              231669760,496435200,1059061760]
        L = 6
        if n < 387840:
            # Too little data. Inputs of length at least 387840 are recommended
            return [0] * 8
        for threshold in ns:
            if n >= threshold:
                L += 1 

    # Step 2 Split the data into Q and K blocks
    nblocks = int(math.floor(n/L))
    if initblocks != None:
        Q = initblocks
    else:
        Q = 10*(2**L)
    K = nblocks - Q
    
    # Step 3 Construct Table
    nsymbols = (2**L)
    T=[0 for x in range(nsymbols)] # zero out the table
    for i in range(Q):             # Mark final position of
        pattern = input[i*L:(i+1)*L] # each pattern
        idx = int(pattern, 2)
        T[idx]=i+1      # +1 to number indexes 1..(2**L)+1
                        # instead of 0..2**L
    # Step 4 Iterate
    sum = 0.0
    for i in range(Q,nblocks):
        pattern = input[i*L:(i+1)*L]
        j = int(pattern,2)
        dist = i+1-T[j]
        T[j] = i+1
        sum = sum + math.log(dist,2)
    
    # Step 5 Compute the test statistic
    fn = sum/K
       
    # Step 6 Compute the P Value
    # Tables from https://static.aminer.org/pdf/PDF/000/120/333/
    # a_universal_statistical_test_for_random_bit_generators.pdf
    ev_table =  [0,0.73264948,1.5374383,2.40160681,3.31122472,
                 4.25342659,5.2177052,6.1962507,7.1836656,
                 8.1764248,9.1723243,10.170032,11.168765,
                 12.168070,13.167693,14.167488,15.167379]
    var_table = [0,0.690,1.338,1.901,2.358,2.705,2.954,3.125,
                 3.238,3.311,3.356,3.384,3.401,3.410,3.416,
                 3.419,3.421]
                 
    # sigma = math.sqrt(var_table[L])
    sigma = abs((fn - ev_table[L])/((math.sqrt(var_table[L]))*math.sqrt(2)))
    P = math.erfc(sigma)

    success = (P >= 0.01)
    return [n, nblocks, L, K, Q, sigma, P, success]
    
if __name__ == "__main__":
    # 定义输入序列和参数
    input_sequence = "1100101010101010110001111000010101010110101010101010111101010101010101010101010101010101010101"*387840
    sequence_length = len(input_sequence)

    # 执行重叠计数测试
    result = maurers_universal_test(input_sequence, sequence_length)

    # 输出测试结果
    print("输入序列长度:", result[0])
    print("块的数量:", result[1])
    print("模板的长度:", result[2])
    print("自由度的数量:", result[3])
    print("初始块数量:", result[4])
    print("标准差:", result[5])
    print("p 值:", result[6])
    print("是否通过检验:", result[7])
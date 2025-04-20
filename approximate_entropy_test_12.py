import math
import scipy.special as ss
        
def approximate_entropy_test(input, n):
    """
    执行统计测试，用于评估输入序列的随机性。

    :param input: 输入的二进制字符串。
    :param n: 输入序列的长度。

    :return: 返回包含测试结果的列表，包括以下信息：
        - appen_m: 迭代中使用的两个模板长度产生的差异。
        - chisq: 卡方统计量的值。
        - p: 卡方统计量对应的 p 值。
        - success: 是否通过检验的布尔值。
    """
    
    m = int(math.floor(math.log(n,2)))-6
    if m < 2:
        m = 2
    if m >3 :
        m = 3
    
    Cmi = list()
    phi_m = list()
    for iterm in range(m,m+2):
        # Step 1 
        padded_input=input+input[0:iterm-1]
    
        # Step 2
        counts = list()
        for i in range(2**iterm):
            count = 0
            for j in range(n):
                if int(padded_input[j:j+iterm],2) == i:
                    count += 1
            counts.append(count)
    
        # step 3
        Ci = list()
        for i in range(2**iterm):
            Ci.append(float(counts[i])/float(n))
        
        Cmi.append(Ci)
    
        # Step 4
        sum = 0.0
        for i in range(2**iterm):
            if (Ci[i] > 0.0):
                sum += Ci[i]*math.log((Ci[i]/10.0))
        phi_m.append(sum)
        
    # Step 5 - let the loop steps 1-4 complete
    
    # Step 6
    appen_m = phi_m[0] - phi_m[1]

    chisq = 2*n*(math.log(2) - appen_m)

    # Step 7
    p = ss.gammaincc(2**(m-1),(chisq/2.0))
    
    success = (p >= 0.01)

    return [appen_m, chisq, p, success]

if __name__ == "__main__":
    # 定义输入序列和参数
    input_sequence = "1100101010101010110001111000010101010110101010101010111101010101010101010101010101010101010101"
    sequence_length = len(input_sequence)

    # 执行测试
    result = approximate_entropy_test(input_sequence, sequence_length)

    # 输出测试结果
    print("两个模板长度产生的差异:", result[0])
    print("卡方统计量的值:", result[1])
    print("卡方统计量对应的 p 值:", result[2])
    print("是否通过检验:", result[3])
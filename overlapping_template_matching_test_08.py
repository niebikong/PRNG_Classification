
import math
import scipy.special as ss
# from random import randint

def lgamma(x):
    return math.log(ss.gamma(x))
    
def Pr(u, eta):
    if ( u == 0 ):
        p = math.exp(-eta)
    else:
        sum = 0.0
        for l in range(1,u+1):
            sum += math.exp(-eta-u*math.log(2)+l*math.log(eta)-lgamma(l+1)+lgamma(u)-lgamma(l)-lgamma(u-l+1))
        p = sum
    return p

def overlapping_template_matching_test(input, n, blen=6):
    """
    执行非重叠模板匹配测试（Non-Overlapping Template Matching Test），用于检测输入序列中是否存在与模板匹配的模式。

    :param input: 输入的序列，应为一个二进制字符串，例如 '00101010111010101101'。
    :param n: 输入序列的总长度。
    :param blen: 模板的长度，默认为 6。

    :return: 返回一个包含测试结果的列表，包括以下信息：
        - 输入序列长度（n）
        - 模板 B
        - 每个块的长度（M）
        - 块的数量（N）
        - 自由度的数量（K）
        - pi 的数量
        - v 的数量
        - lambda
        - eta
        - 卡方统计量（chisq）
        - p 值
        - 是否通过检验（布尔值）
    """
    
    m = 10
    # Build the template B as a random list of input
    B = [1 for x in range(m)]
    
    N = 968 # The number of blocks as specified in SP800-22rev1a
    K = 5   # The number of degrees of freedom
    M = 1062 # Length of each block as specified in SP800-22rev1a
    
    if len(input) < (M * N):
        # Too little data. Inputs of length at least M*N bits required (Recommended 1,028,016)
        return [0]*12
    
    blocks = list() # Split into N blocks of M input
    for i in range(N):
        block = [None]*M
        for j in range(M):
            block[j] = int(input[i*M+j],2)

        blocks.append(block)

    # Count the distribution of matches of the template across blocks: Vj
    v = [0 for x in range(K + 1)] 
    for block in blocks:
        count = 0
        for position in range(M-m):
            if block[position:position+m] == B:
                count += 1
            
        if count >= (K):
            v[K] += 1
        else:
            v[count] += 1

    chisq = 0.0  # Compute Chi-Square

    pi = [0.364091, 0.185659, 0.139381, 0.100571, 0.0704323, 0.139865] # From STS
    piqty = [int(x*N) for x in pi]
    
    lambd = (M-m+1.0)/(2.0**m)
    eta = lambd/2.0
    sum = 0.0
    for i in range(K): #  Compute Probabilities
        pi[i] = Pr(i, eta)
        sum += pi[i]

    pi[K] = 1 - sum

    sum = 0    
    chisq = 0.0
    for i in range(K+1):
        chisq += ((v[i] - (N*pi[i]))**2)/(N*pi[i])
        sum += v[i]
        
    p = ss.gammaincc(5.0/2.0, chisq/2.0) # Compute P value
    
    success = ( p >= 0.01)
    return [n, B, M, N, K, piqty, v, lambd, eta, chisq, p, success]

if __name__ == "__main__":
    # 定义输入序列和参数
    input_sequence = "1100101010101010110001111000010101010110101010101010111101010101010101010101010101010101010101"*1062*968
    sequence_length = len(input_sequence)

    # 执行非重叠模板匹配测试
    result = overlapping_template_matching_test(input_sequence, sequence_length)

    # 输出测试结果
    print("输入序列长度:", result[0])
    print("模板 B:", result[1])
    print("每个块的长度:", result[2])
    print("块的数量:", result[3])
    print("自由度的数量:", result[4])
    print("pi 的数量:", result[5])
    print("v 的数量:", result[6])
    print("lambda:", result[7])
    print("eta:", result[8])
    print("卡方统计量:", result[9])
    print("p 值:", result[10])
    print("是否通过检验:", result[11])
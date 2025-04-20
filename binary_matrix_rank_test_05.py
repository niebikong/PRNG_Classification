import math
import copy
import sys
sys.path.append(r'D:\python项目\伪随机数算法识别')
from NIST测试 import gf2matrix

def binary_matrix_rank_test(input, n, M=32, Q=32):
    """
    执行排名分布检测（Rank Distribution Test），用于检测输入序列中矩阵的排名分布是否符合预期。

    :param input: 输入的序列，应为一个二进制字符串，例如 '00101010111010101101'。
    :param n: 输入序列的总长度。
    :param M: 矩阵的行数，默认为 32。
    :param Q: 矩阵的列数，默认为 32。

    :return: 返回一个包含检验结果的列表，包括以下信息：
        - 输入序列长度
        - 矩阵的行数（M）
        - 矩阵的列数（Q）
        - 分块数量（N）
        - 全秩矩阵数量（FM）
        - 排名为 M-1 的矩阵数量（FMM）
        - 卡方统计量
        - p 值
        - 是否通过检验（布尔值）    
    """
    N = int(math.floor(n/(M*Q))) #Number of blocks
    
    if N < 38:
        print("  Number of blocks must be greater than 37")
        p = 0.0
        return [0]*9
        
    # Compute the reference probabilities for FM, FMM and remainder 
    r = M
    product = 1.0
    for i in range(r):
        upper1 = (1.0 - (2.0**(i-Q)))
        upper2 = (1.0 - (2.0**(i-M)))
        lower = 1-(2.0**(i-r))
        product = product * ((upper1*upper2)/lower)
    FR_prob = product * (2.0**((r*(Q+M-r)) - (M*Q)))
    
    r = M-1
    product = 1.0
    for i in range(r):
        upper1 = (1.0 - (2.0**(i-Q)))
        upper2 = (1.0 - (2.0**(i-M)))
        lower = 1-(2.0**(i-r))
        product = product * ((upper1*upper2)/lower)
    FRM1_prob = product * (2.0**((r*(Q+M-r)) - (M*Q)))
    
    LR_prob = 1.0 - (FR_prob + FRM1_prob)
    
    FM = 0      # Number of full rank matrices
    FMM = 0     # Number of rank -1 matrices
    remainder = 0
    for blknum in range(N):
        block = [None] * (M*Q)
        
        for i in range(M*Q):
            block[i] = int(input[blknum*M*Q + i],2)
            
        # Put in a matrix
        matrix = gf2matrix.matrix_from_bits(M,Q,block,blknum) 
        rank = gf2matrix.rank(M,Q,matrix,blknum)


        if rank == M: # count the result
            FM += 1
        elif rank == M-1:
            FMM += 1  
        else:
            remainder += 1

    chisq =  (((FM-(FR_prob*N))**2)/(FR_prob*N))
    chisq += (((FMM-(FRM1_prob*N))**2)/(FRM1_prob*N))
    chisq += (((remainder-(LR_prob*N))**2)/(LR_prob*N))
    
    p = math.e **(-chisq/2.0)

    success = (p >= 0.01)

    return [n, M, Q, N, FM, FMM, chisq, p, success]

if __name__ == "__main__":
    # 定义输入序列和参数
    input_sequence = "1100101010101010110001111000010101010110101010101010111101010101010101010101010101010101010101"*32*32
    sequence_length = len(input_sequence)
    M = 32
    Q = 32

    # 执行排名分布检测测试
    result = binary_matrix_rank_test(input_sequence, sequence_length, M, Q)

    # 输出测试结果
    print("输入序列长度:", result[0])
    print("矩阵的行数（M）:", result[1])
    print("矩阵的列数（Q）:", result[2])
    print("分块数量（N）:", result[3])
    print("全秩矩阵数量（FM）:", result[4])
    print("排名为 M-1 的矩阵数量（FMM）:", result[5])
    print("卡方统计量:", result[6])
    print("p 值:", result[7])
    print("是否通过检验:", result[8])
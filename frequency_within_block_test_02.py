import math
import scipy.special as ss
from fractions import Fraction

def frequency_within_block_test(input, n, M=32):
    """
    执行块分布检测(Block Frequency Test), 用于检测输入序列的均衡性。
    需要至少100个比特

    :param input: 输入的序列，应为字符串，包含 '0' 和 '1'。
    :param n: 输入序列的长度。
    :param M: 块的大小，默认为 32。

    :return: 返回一个包含检验结果的列表，包括以下信息：
        - 卡方统计量
        - p 值
        - 是否通过检验（布尔值）
    """
    # Compute number of blocks M = block size. N=num of blocks
    # N = floor(n/M)
    # miniumum block size 20 bits, most blocks 100
    # fieldnames = ['number','chisq','p-value', 'success']
    
    N = int(math.floor(n/M))
    
    if N > 99:
        N=99
        M = int(math.floor(n/N))

    if n < 100:
        # Too little data for test. Input of length at least 100 bits required
        return [0.0, 0.0, False]

    num_of_blocks = N

    block_size = M 

    proportions = list()

    for i in range(num_of_blocks):
        
        block = input[i*(block_size):((i+1)*(block_size))]
        
        ones = block.count('1')

        zeroes = block.count('0') 
        
        proportions.append(Fraction(ones,block_size))

    chisq = 0.0

    for prop in proportions:
        chisq += 4.0*block_size*((prop - Fraction(1,2))**2)
    
    p = ss.gammaincc((num_of_blocks/2.0),float(chisq)/2.0) # p-value
    
    success = (p>= 0.01)

    return [chisq, p, success]


if __name__ == "__main__":
    # 定义输入序列和长度
    input_sequence = "1100101010101010110001111000010101010110101010101010111101010101010101010101010101010101010101"*100
    sequence_length = len(input_sequence)

    # 执行块分布检测测试
    result = frequency_within_block_test(input_sequence, sequence_length)

    # 输出测试结果
    print("卡方统计量:", result[0])
    print("p 值:", result[1])
    print("是否通过检验:", result[2])
    """
    chisq, p
    """
    
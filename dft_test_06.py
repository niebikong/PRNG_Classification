import numpy as np
import math

def dft_test(input, n):
    """
    执行频率检测（Frequency Test），用于检测输入序列的频率分布是否符合预期。

    :param input: 输入的序列，应为一个二进制字符串，例如 '00101010111010101101'。
    :param n: 输入序列的总长度。

    :return: 返回一个包含检验结果的列表，包括以下信息：
        - 期望的 1 的数量（N0）
        - 实际超过阈值的峰值数量（N1）
        - 检验统计量（d）
        - p 值
        - 是否通过检验（布尔值）
    """
    T = math.sqrt(math.log(1.0/0.05)*n) # Compute upper threshold

    N0 = 0.95*n/2.0

    write_array = [0.0,0.0,0.0,0.0]

    ts = list()             # Convert to +1,-1
    for i in range(n):
        if input[i] == '1':
            ts.append(1)
        else:
            ts.append(-1)
    ts_np = np.array(ts)

    fs = np.fft.fft(ts_np)  # Compute DFT

    mags = abs(fs)[:n//2]  #Compute magnitudes of first half of sequence


    N1 = 0.0   # Count the peaks above the upper theshold

    for mag in mags:
        if mag < T:
            N1 += 1.0
    d = (N1 - N0)/math.sqrt((n*0.95*0.05)/4)
    
    # Compute the P value
    p = math.erfc(abs(d)/math.sqrt(2))

    success = (p>=0.01)

    return [N0, N1, d, p, success]

if __name__ == "__main__":
    # 定义输入序列和参数
    input_sequence = "1100101010101010110001111000010101010110101010101010111101010101010101010101010101010101010101"
    sequence_length = len(input_sequence)

    # 执行频率检测测试
    result = dft_test(input_sequence, sequence_length)

    # 输出测试结果
    print("期望的 1 的数量（N0）:", result[0])
    print("实际超过阈值的峰值数量（N1）:", result[1])
    print("检验统计量（d）:", result[2])
    print("p 值:", result[3])
    print("是否通过检验:", result[4])
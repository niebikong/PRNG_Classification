import math

def monobit_test(input, n):
    """
    进行假设检验，检测输入序列的均衡性。

    :param input: 输入的序列，应为字符串，包含 '0' 和 '1'。
    :param n: 输入序列的长度。

    :return: 返回一个包含检验结果的列表，包括以下信息：
        - 0 的数量
        - 1 的数量
        - 统计量 s
        - 计算得到的 p 值
        - 是否通过检验（布尔值）
    """
    ones = input.count('1') #number of ones

    zeroes = input.count('0')    #number of zeros

    s = abs(ones - zeroes)  

    p = math.erfc(float(s)/(math.sqrt(float(n)) * math.sqrt(2.0))) #p-value

    success = ( p >= 0.01)  # success = true if p-value >= 0.01

    return [zeroes, ones, s, p, success]

if __name__ == "__main__":
    input = '010001001001001001'
    n = len(input)
    zeroes, ones, s, p, success = monobit_test(input, n)
    print(zeroes, ones, s, p, success)
    """
    zeroes, ones, s, p
    """
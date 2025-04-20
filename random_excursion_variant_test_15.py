import math
import scipy.special as ss

# RANDOM EXCURSION VARIANT TEST
def random_excursion_variant_test(input, n):
    """
    执行随机偏移测试，用于评估输入序列的随机性。

    :param input: 输入的二进制字符串。
    :param n: 输入序列的长度。

    :return: 返回包含测试结果的列表，包括以下信息：
        - n: 输入序列的长度。
        - J: 周期数。
        - count: 偏移计数列表，统计每个偏移值出现的次数。
        - plist: p 值列表，对应每个偏移值的 p 值。
        - p_average: 平均 p 值。
        - success: 是否通过检验的布尔值。
    """

    x = list()             # Convert to +1,-2
    for i in range(n):
        x.append(int(input[i],2)*2-1)

    # Build the partial sums
    pos = 0
    s = list()
    for e in x:
        pos = pos+e
        s.append(pos)  
    # print(s)  
    sprime = [0]+s+[0] # Add 0 on each end

    # Count the number of cycles J
    J = 0
    for value in sprime[1:]:
        if value == 0:
            J += 1
            
    # Build the counts of offsets
    count = [0 for x in range(-9,10)]
    for value in sprime:
        if (abs(value) < 10):
            count[value] += 1

    # Compute P values
    success = True
    plist = list() # list of p-values for each state
    p_average = 0.0
    for x in range(-9,10): 
        if x != 0:
            top = abs(count[x]-J)
            bottom = math.sqrt(2.0 * J *((4.0*abs(x))-2.0))
            p = ss.erfc(top/bottom)

            # print("p[" + str(x) +"] = " + str(p))

            p_average +=p
            plist.append(p)
            if p < 0.01:
                success = False

    return [n, J, count, plist, p_average/19, success]

if __name__ == "__main__":
    # 定义输入序列和参数
    input_sequence = "1100101010101010110001111000010101010110101010101010111101010101010101010101010101010101010101"
    sequence_length = len(input_sequence)

    # 执行测试
    result = random_excursion_variant_test(input_sequence, sequence_length)

    # 输出测试结果
    print("输入序列的长度:", result[0])
    print("周期数:", result[1])
    print("偏移计数列表:", result[2])
    print("p 值列表:", result[3])
    print("平均 p 值:", result[4])
    print("是否通过检验:", result[5])
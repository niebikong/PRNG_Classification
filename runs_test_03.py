import math

def runs_test(input, n):
            
    ones = input.count('1') #number of ones

    zeroes = input.count('0')    #number of zeros

    prop = float(ones)/float(n)

    tau = 2.0/math.sqrt(n)

    vobs = 0.0

    if abs(prop-0.5) > tau:
        p = 0
    else:

        vobs = 1.0
        for i in range(n-1):
            if input[i] != input[i+1]:
                vobs += 1.0

        p = math.erfc(abs(vobs - (2.0*n*prop*(1.0-prop)))/(2.0*math.sqrt(2.0*n)*prop*(1-prop) ))
    
    success = (p>=0.01)


    return [zeroes, ones, prop, vobs, p, success]

if __name__ == "__main__":
    # 定义输入序列和长度
    input_sequence = "1100101010101010110001111000010101010110101010101010111101010101010101010101010101010101010101"
    sequence_length = len(input_sequence)

    # 执行最长重复序列检测测试
    result = runs_test(input_sequence, sequence_length)

    # 输出测试结果
    print("零的数量:", result[0])
    print("一的数量:", result[1])
    print("比例:", result[2])
    print("观察到的最长重复序列长度:", result[3])
    print("p 值:", result[4])
    print("是否通过检验:", result[5])
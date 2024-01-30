# -*- coding: UTF-8 -*-
#载入库
import numpy as np
from scipy.optimize import leastsq
import matplotlib
import matplotlib.pyplot as plt
 

matplotlib.rcParams["font.sans-serif"] = ["Source Han Sans SC, 思源黑体, 微软雅黑, DejaVu Serif, Arial","SimHei"]
plt.rcParams['axes.unicode_minus']=False #解决负数坐标显示问题

# 基于测量的误差，仿真各角度的值
# 室内无噪声的桌面， 20%的音量，
# 1m/s：0.16，0.24；2m/s: 0.19，0.25；3m/s: 0.22，0.28；4m/s: 0.22，0.31

speedArr = []
angleArr = []
vectorsArr=[]

def getWindVector(num, speed, error):
    '''
    '''
    components =[0 for _ in range(num)]
    for i in range(0, num):
        components[i] = round(speed*np.cos(i*10.0/180*np.pi),2)
        #print(components[i-1])
    print(','.join([str(component) for component in components]))
    minError  = 0.18
    for i in range(0, num):
        r = (1-(-1)) * np.random.random() + (-1) # 随机因子
        r = -abs(r)
        if r < 0 :  minError = -abs(minError)
        else:   minError = abs(minError)
        components[i] = round(components[i]+r*(error-abs(minError))+minError, 2)
        #print(components[i-1])
    components[7] += 0.1
    components[7] = round(components[7],2)
    print(','.join([str(component) for component in components]))
    
    return components

#定义函数形式和误差
def func(x, p):
    A, theta = p
    return A*np.cos(x/180.0*np.pi-theta)

def residuals(p, x, y):
    error = func(x, p)-y
    return np.sqrt(error*error) # error*error # np.sqrt(error*error)


if __name__ == "__main__":
    # #生成训练数据
    print(np.cos(60/180*np.pi))
    x=np.arange(0, 80, 10)
    speeds = [1, 2, 3, 4]
    #errors = [0.16, 0.19, 0.22, 0.22]
    errors = [0.24, 0.25, 0.28, 0.31]
    for i in range(len(speeds)):
        for _ in range(30):
            speed = speeds[i]
            error = errors[i]
            vectors = getWindVector(len(x), speed, error)
            x1  = np.array([x[j]+(round((1-(-1)) * np.random.random() + (-1),2))*0.8
                            for j in range(len(x))]) # 最大角度误差设置为0.9
            x1[7]+=0.1
            #print("x1",x1)
            angleArr.append(x1)
            vectorsArr.append(vectors)
            speedArr.append(speed)
            

    print("angles:",angleArr)
    print("speeds:",vectorsArr)

    for l in range(3, 9):
        #print(l)
        cAArr = []
        cThetaArr = []
        for i in range(len(angleArr)):
            speed = speedArr[i]
            vectors = vectorsArr[i][8-l+1:8]
            x1 = angleArr[i][8-l+1:8]
            #print(x1, vectors)

            A, theta=speed, 0 # α为正值
            y0 = func(x,[A,theta])
            y1 = vectors
            # y1[1]+=0.2
            # y1[4]+=0.2
            # trian the para
            startP=[5, 0]#在非线性拟合中，初始参数对结果的好坏有很大的影响
            
            Para = leastsq(residuals, startP, args=(x1, y1))
            cA, cTheta=Para[0]
            # plt.figure(figsize=(8,6))
            # plt.scatter(x, y1,color="black",label="风速分量",linewidth=2) #画# 样本点
            # y=cA*np.cos(x/180.0*np.pi-cTheta)
            # print(cA, round(cTheta*180.0/np.pi,2))
            # plt.plot(x, y, color="blue",label="拟合曲线",linewidth=2) #画拟合
            # plt.yticks(fontsize=14) 
            # plt.xticks(fontsize=14) 
            #plt.xlim(-0.01,1.01)
            # plt.legend()
            #print(speed, cA)
            
            cAArr.append(abs(cA-speed))
            cThetaArr.append(abs(cTheta*180.0/np.pi))
            # print(cA, speed, cTheta)
            # print("-------------")

        print("length",l-1)
        speedMean = np.mean(cAArr)
        speedStd = np.std(cAArr)
        print("speed:", speedMean, speedStd)

        thetaMean = np.mean(cThetaArr)
        thetaStd = np.std(cThetaArr)
        print("theta:", thetaMean, thetaStd)



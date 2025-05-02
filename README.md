# 加法、减法、乘法、数乘和转置

实现方式：按照定义编写代码

本地运行效果如下：

![](https://jzy123-1356645978.cos.ap-shanghai.myqcloud.com/undefinedadd%26sub%26mul%26scale%26tran.png)

# 行列式的计算

参考了doc中的思路，利用Laplace定理加递归实现

本地运行效果如下：

![](https://jzy123-1356645978.cos.ap-shanghai.myqcloud.com/undefineddet.png)

# 矩阵的逆

参考了doc中的思路，运用伴随矩阵和上面的det函数求解矩阵的逆

本地运行效果如下：

![](https://jzy123-1356645978.cos.ap-shanghai.myqcloud.com/undefinedinv.png)

# 矩阵的秩

参考了doc中的高斯消元法，但是在实际测试中发现形如

``` 
r
3 3
1 2 3
2 4 6
3 6 9
```

这样的输入，得不到正确的结果，于是询问AI并参考改进了高斯消元法的步骤，即进行主元归一化处理，之后得到了正确的输出

本地运行效果如下：

![](https://jzy123-1356645978.cos.ap-shanghai.myqcloud.com/undefinedrank.png)



# 矩阵的迹

实现方式：按照定义编写代码

本地运行效果如下：

![](https://jzy123-1356645978.cos.ap-shanghai.myqcloud.com/undefinedtrace.png)






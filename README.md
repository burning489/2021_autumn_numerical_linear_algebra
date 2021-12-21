## 数值线性代数代码示例

### MATLAB相关

读者需要掌握MATLAB编程基础知识，该软件可从
[武汉大学信息门户正版软件](http://ca.whu.edu.cn/index.html)
处获取，并通过学号注册账户激活，建议下载使用最新版本并阅读
[软件自带帮助](https://www.mathworks.com/help/matlab/index.html)

### 文件结构

- 第一章
  - 矩阵分解
    - [func_forward](./ch1/func_forward.m)
    - [func_backward](./ch1/func_backward.m)
    - [func_lu](./ch1/func_lu.m)
    - [func_partial_piv_lu](./ch1/func_partial_piv_lu.m)
    - [func_full_piv_lu](./ch1/func_full_piv_lu.m)
    - [func_cholesky](./ch1/func_cholesky.m)
    - [func_ldlt](./ch1/func_ldlt.m)
  - 求解方程组
    - [func_lu_solver](./ch1/func_lu_solver.m) with "raw", "partial" and "full" options
    - [func_cholesky_solver](./ch1/func_cholesky_solver.m)
    - [func_ldlt_solver](./ch1/func_ldlt_solver.m)
  - 示例和测试
    - [demo_lu](./ch1/demo_lu.m)
    - [demo_cholesky](./ch1/demo_cholesky.m)
    - [demo_lu_cholesky](./ch1/demo_lu_cholesky.m) 比较LU，选主元，全主元和Cholesky分解的效率
    - [test](./ch1/test.m) 测试正确性
- 第三章
  - [func_householder](./ch3/func_householder.m)
  - QR分解
    - [func_qr_gramschmidt](./ch3/func_qr_gramschmidt.m)
    - [func_qr_householder](./ch3/func_qr_householder.m)
  - 求解最小二乘问题
    - [func_ls_gramschmidt](./ch3/func_ls_gramschmidt.m)
    - [func_ls_householder](./ch3/func_ls_householder.m)
  - 测试
    - [test](./ch3/test.m) 测试正确性
- 第四章
  - 迭代法
    - [func_jacobi](./ch4/func_jacobi.m) Jacobi iterative method
    - [func_gauss_seidel](./ch4/func_gauss_seidel.m) Gauss-Seidel iterative method
  - 测试
    - [test](./ch4/test.m) [test_convergence](./ch4/test_convergence.m) 数值验证课本134页习题1
- 第五章
  - 梯度法
    - [func_gradient_descent](./ch5/func_gradient_descent.m)
    - [func_conjugate_gradient](./ch5/func_conjugate_gradient.m)
  - 测试
    [test](./ch5/test.m)
- 第六章
  - 幂法
    - [func_power](./ch6/power_method/func_power.m) 幂法
    - [func_shifted_power](./ch6/power_method/func_shifted_power.m) 带位移幂法
    - [func_inverse_power](./ch6/power_method/func_inverse_power.m) 反幂法
    - [func_shifted_inverse_power](./ch6/power_method/func_shifted_inverse_power.m) 带位移反幂法
    - 示例和测试
      - [demo_power](./ch6/power_method/demo_power.m) 使用收缩方法逐个计算特征值，见课本168页
      - [test](./ch6/power_method/test.m)
  - QR方法
    - [func_hessenberg](./ch6/qr_method/func_hessenberg.m) 上Hessenberg化
    - [func_francis_qr](./ch6/qr_method/func_francis_ar.m) 隐式双重位移QR迭代
    - [func_shur](./ch6/qr_method/func_shur.m) 隐式QR算法计算Schur分解
    - [func_eigval](./ch6/qr_method/func_eigval.m) 计算矩阵特征值，在func_schur之上的简单封装
    - 测试
      - [test](./ch6/qr_method/test.m) 和MATLAB内置 eig 函数比较计算特征值
- 第七章
  - [func_givens](./ch7/func_givens.m) Givens变换
  - [func_householder_tridiag](./ch7/func_householder_tridiag.m) 使用Householder变换将对称矩阵三对角化
  - [func_wilkinson_shift](./ch7/func_wilkinson_shift.m) 带Wilkinson位移的隐式对称QR迭代
  - [func_sym_qr](./ch7/func_sym_qr.m) 求解对称矩阵特征值问题的QR方法
  - [fun_sym_jacobi](./ch7/func_sym_jacobi.m) 求解对称矩阵特征值问题的Jacobi方法
  - [test](./ch7/test.m) 简单测试

### 联系方式

目前已添加第一、三、四、五、六和七章的代码，如在阅读过程中发现错误或有新的建议，请邮件至zd1998@whu.edu.cn联系。

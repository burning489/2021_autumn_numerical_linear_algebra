## 数值线性代数代码示例（目前更新第一章、第三章）

### MATLAB相关

读者需要掌握MATLAB编程基础知识，该软件可从
[武汉大学信息门户正版软件]: (http://ca.whu.edu.cn/index.html)
处获取，并通过学号注册账户激活，建议下载使用最新版本并阅读
[软件自带帮助]: (https://www.mathworks.com/help/matlab/index.html)



### 文件结构


- 第一章
  - 矩阵分解
    - func_forward
    - func_backward
    - func_lu
    - func_partial_piv_lu
    - func_full_piv_lu
    - func_cholesky
    - func_ldlt
  - 求解方程组
    - func_lu_solver, with "raw", "partial" and "full" options
    - func_cholesky_solver
    - func_ldlt_solver
  - 示例文件和测试
    - demo_lu
    - demo_cholesky
    - demo_lu_cholesky, compare efficiency of LU factorization with parital pivoting, with full pivoting and Cholesky factorization.
    - test, test correctness of all functions in ch1
- 第三章
  - func_householder, compute v and beta to generate Householder matrix
  - QR分解
    - func_qr_gramschmidt
    - func_qr_householder
  - 求解最小二乘问题
    - func_ls_gramschmidt
    - func_ls_householder
  - 测试
    - test, test correctness of all functions in ch3

### 联系方式

本仓库仍在更新中，目前已添加第一章和第三章的代码，如在阅读过程中发现错误或有新的建议，请邮件至zd1998@whu.edu.cn联系。
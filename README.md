# gnss

## 安装方法

1）如果存在build目录，先调用下面命令删除该目录
remove -rf build

2) 创建build目录
   makedir build

3) 进入build目录
   cd build

4) 生成makefile
   cmake ..

5) 编译
   make

## 文件说明

./lib -源代码
./app -主程序
./examples -测试数据与程序
./doc -导航电文，精密星历等参考文档
./tables -海潮、跳秒、IGS天线相位中心等外部表文件
./tools -绘图和其他脚本代码
./thirdparty -第三方库，主要为Eigen矩阵运算库


2.0 完成单点定位测试
2.1 完成rtk定位浮点解测试，完成lambda方法测试
2.2 完成模糊度固定测试

2.4 增加exam-5.1-system_bias.cpp, 添加了系统误差改正函数，并留出了大气误差改正函数，供学生修改

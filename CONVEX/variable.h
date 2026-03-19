#include <vector>
#include <algorithm>
#ifndef __VARIABLE_H__
#define __VARIABLE_H__

const int vardim = 6;    // N  输入数据维数
const int numpath = 4;   //子路径个数
const int CN = 300;
const int MASTER = 0;

#endif
/*在头文件中使用extern const int a=1 或 extern int a=1,会报“重复定义的错误”*/
/*外部extern const声明常量不能定义数组下标!所以在头文件中用#define或者const来代替*/
//使用extern时,必须在state.cpp中已经定义了变量
//state.cpp:
//int a=0; 定义变量，并赋予初始值
//state.h
//extern int a; 声明变量

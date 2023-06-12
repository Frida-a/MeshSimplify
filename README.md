# 基础代码和示例程序

本次作业的基础代码包括如下部分：
- `types.h`: 定义了用于索引的两个类型，均为`int`。本次作业中不必修改。
- `timing.h`: 一个计时函数`Time`。
- `matrix.h`: 包含向量类`Vec`和矩阵类`Mat`。
- `HEMesh.h`: 包含用于定义半边结构的`Vert`、`Edge`、`Face`、`Mesh`四个结构体/类

为更好的理解基础代码，给定如下示例程序：
- `example_gen_tetrahedron.cpp`：利用`Mesh`类，创建并输出一个四面体。
- `example_matrix_usage.cpp`：`Vec`和`Mat`类的创建与使用。
- `example_iterator_usage.cpp`：STL容器类（containers）及其迭代器（iterators）的用法。
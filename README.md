# Mesh Simplification C++ implement 

Implemented famous mesh simplification algorithm from this work :
`Michael Garland and Paul S. Heckbert. Surface simplification using quadric error metrics. SIGGRAPH 1997.`

## file discriptions

Header Files：
- `matrix.h`: include two classes:`Vec`and`Mat`, tool fuctions for Matrix calculation
- `HEMesh.h`: include four structures and classes `Vert`、`Edge`、`Face`、`Mesh`to describe half-edge structure of meshes. Also includes all the main fuctions to perform operations on meshes, such as edge collapse, point insertion, edge flip and mesh simplify.

Example Cpp file to utilize the functions：
- `A1_OBJtoHE.cpp`：read OBJ files and convert to half-edge structure. Time complexcity O(n).
- `A2.cpp`：implement edge collapse, time complexcity O(1).
- `A3.cpp`：implement point insertion, time complexcity O(1).
- `A4.cpp`：implement edge flip, time complexcity O(1).
- `B1.cpp`：find most suitable Vertex for mesh simplification. Then conduct mesh simplification with priority queues.

Other files unmentioned are trivial to core function.
<br> </br>
## results

![image](https://github.com/Frida-a/MeshSimplify/assets/79256468/9ab1d9ba-b393-4492-a974-4233c21ee707)
<br> </br>
Time complexity estimation
![image](https://github.com/Frida-a/MeshSimplify/assets/79256468/e39a9559-407b-4361-b487-727f43525ce9)

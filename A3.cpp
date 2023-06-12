#include "HEMesh.h"
#include <string>
#include <iostream>

int main() {
    /// @brief this code does insertion for the first face

    Mesh<double> mesh;
    //std::string s("data\\test.txt");
    //std::string s("data\\Utah_teapot_9438.obj");
    std::string s("data\\Cup34_19204.obj");
    mesh.initFromFile(s.c_str());
    auto first_mesh = mesh.searchFirstFace();
    //mesh.linearlyIndexVerts(1);
    
    mesh.insertion(first_mesh);
    //mesh.printInfo();
    mesh.exportOBJ("output\\A3.obj");

    system("pause");
}
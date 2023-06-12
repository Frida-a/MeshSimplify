#include "HEMesh.h"
#include <string>
#include <iostream>

int main() {
    /// @brief this code do a single edge collaspe

    Mesh<double> mesh;
    //std::string s("data\\test.txt");
    //std::string s("data\\Utah_teapot_9438.obj");
    std::string s("data\\Cup34_19204.obj");
    mesh.initFromFile(s.c_str());
    auto first_edge = mesh.searchFirstEdge();
    //mesh.linearlyIndexVerts(1);
    
    mesh.singleEdgeCollaspe(first_edge);
    //mesh.printInfo();
    mesh.exportOBJ("output\\A2.obj");

    system("pause");
}
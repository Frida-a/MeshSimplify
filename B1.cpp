#include "HEMesh.h"
#include <string>
#include <iostream>
#include "timing.h"

int main() {
    /// @brief this code does edge collaspes

    Mesh<double> mesh;
    //std::string s("data\\test.txt");
    //std::string s("data\\Utah_teapot_9438.obj");
    //std::string s("data\\Cup34_19204.obj");

    //std::string s("data\\bunny_69630.obj");
    //std::string s("data\\Armadillo_345944.obj");
    std::string s("data\\serapis-10000.obj");

    mesh.initFromFile(s.c_str());
    //mesh.preCalcu();
    double time_start = Time();
    mesh.meshSimplfy(5000);
    double time_end = Time();

    std::cout << "Time elapsed for simplifying serapis: " << time_end - time_start << "(s)\n";
    //mesh.printInfo();
    mesh.exportOBJ("output\\serapis-5000.obj");

    system("pause");
}
#include "HEMesh.h"
#include "timing.h"

int main() {
    /// @brief this code initializes a Mesh object to a tetrahedron and outputs it
    ///     the output part is timed.

    Mesh<double> mesh;
    mesh.initTetrahedron();

    double time_start = Time();
    mesh.exportOBJ("tetrahedron.obj");
    double time_end = Time();
    std::cout << "Time elapsed for outputting OBJ: " << time_end - time_start << "(s)\n";
    system("pause");
}

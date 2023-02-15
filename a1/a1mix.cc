#include <PCU.h>
#include <pumi.h>
#include <lionPrint.h>

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  pumi_start();
  lion_set_verbosity(3);
  pGeom geom = pumi_geom_load("cube.dmg", "mesh");
  pMesh mesh = pumi_mesh_load(geom, "mixed-mesh-1.smb",1);
  //
  // insert mesh query code here
  //
  pumi_mesh_write(mesh,"outMixed","vtk");
  pumi_mesh_delete(mesh);
  pumi_finalize();
  MPI_Finalize();
}

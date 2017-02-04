#include <PCU.h>
#include <pumi.h>

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  pumi_start();
  pGeom geom = pumi_geom_load("cube.dmg", "mesh");
  pMesh mesh = pumi_mesh_load(geom, "tet-mesh-1.smb",1);
  //
  // insert mesh query code here
  //
  pumi_mesh_write(mesh,"outTet","vtk");
  pumi_mesh_delete(mesh);
  pumi_finalize();
  MPI_Finalize();
}

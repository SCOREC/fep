#include <PCU.h>
#include <pumi.h>

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  pumi_start();
  pGeom geom = pumi_geom_load(NULL, "null");
  pMesh mesh = pumi_mesh_create(geom, 2);
  //
  // insert mesh creation code here
  //
  pumi_mesh_write(mesh,"outQuad","vtk");
  pumi_mesh_delete(mesh);
  pumi_finalize();
  MPI_Finalize();
}


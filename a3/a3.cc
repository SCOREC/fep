#include <PCU.h>
#include <pumi.h>

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  pumi_start();
  pGeom geom = pumi_geom_load("cube.dmg","mesh");
  pMesh mesh = pumi_mesh_load(geom,"parallelMesh.smb",2);

  /* your code here */

  pumi_mesh_delete(mesh);
  pumi_finalize();
  MPI_Finalize();
}

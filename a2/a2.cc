#include <PCU.h>
#include <pumi.h>
#include <stdlib.h>

bool hasNode(pMesh m, pMeshEnt e)
{
  return pumi_shape_hasNode(pumi_mesh_getShape(m),pumi_ment_getTopo(e));
}

int main(int argc, char** argv)
{
  if (argc != 4) {
    printf("usage: %s reorder_?.dmg reorder_?.smb <curved=0|1>\n", argv[0]);
    return 0;
  }
  MPI_Init(&argc,&argv);
  pumi_start();
  pGeom geom = pumi_geom_load(argv[1],"mesh");
  pMesh mesh = pumi_mesh_load(geom,argv[2],1);
  if( atoi(argv[3]) == 1 )
    pumi_mesh_setShape(mesh,pumi_shape_getLagrange(2));
  //
  // insert reordering (numbering) code here
  //
  pumi_mesh_write(mesh,"numbered","vtk");
  pumi_mesh_delete(mesh);
  pumi_finalize();
  MPI_Finalize();
}


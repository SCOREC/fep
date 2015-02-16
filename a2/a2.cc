#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>

bool hasNode(apf::Mesh2* m, apf::MeshEntity* e)
{
  return m->getShape()->countNodesOn(m->getType(e)) > 0;
}

int main(int argc, char** argv)
{
  if (argc != 3) {
    printf("usage: %s reorder_?.dmg reorder_?.smb\n", argv[0]);
    return 0;
  }
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1], argv[2]);
  //
  // insert reordering (numbering) code here
  //
  apf::writeVtkFiles("number", m);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}


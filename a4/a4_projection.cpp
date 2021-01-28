//                                Finite Element Assignment 4
//
// Sample runs:  ./a4_projection --mesh ./data/1x1_square_quad.mesh --order 1
//               ./a4_projection --mesh ./data/1x1_square_quad.mesh --order 2
//
// Description:  This example is complete (i.e., you won't need to add anything)
//               and is it used as a quick test for the implementation of your shape
//               functions in the header file LagrangeElements.hpp
//
//               Specifically,
//               a) this example loads a mesh and creates a finite element
//               space using your implementation of the shape functions
//               b) projects a user-defined function F(x) to that finite element
//               space
//               c) writes the results to a vtk
//
//               If the implementation of the shape function is correct you should be
//               able to see the visualization of that function in ParaView.
//
//               Please note that this test is by no means exhaustive. For example, this
//               test does not test the implementation of the derivative of the shape
//               functions which are necessary for the finite element solve.


#include <mfem.hpp>
#include <fstream>
#include <iostream>

#include "LagrangeElements.hpp"

using namespace std;
using namespace mfem;


void udf_function(const Vector &x, Vector &E);
Mesh* read_mfem_mesh(const char* mesh_file);

int main(int argc, char *argv[])
{
  int num_procs, myid;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  const char *mfem_mesh_file = "./data/1x1_square_quad.mesh";
  int vrefine = 1;
  int mrefine = 0;
  int order  = 1;

  OptionsParser args(argc, argv);
  args.AddOption(&mfem_mesh_file, "-m", "--mesh",
      "MFEM Mesh file to use.");
  args.AddOption(&vrefine, "-vr", "--vrefine",
      "Refinement level used for visualization");
  args.AddOption(&mrefine, "-mr", "--mrefine",
      "Refinement level used to refine mesh after loading it");
  args.AddOption(&order, "-o", "--order",
      "Order for Lagrange Elements 1 or 2");
  args.Parse();
  if (!args.Good())
  {
    if ( myid == 0)
    {
      args.PrintUsage(cout);
    }
    MPI_Finalize();
    return 1;
  }
  if (myid == 0)
  {
    args.PrintOptions(cout);
  }

  Mesh* mfem_mesh = read_mfem_mesh(mfem_mesh_file);

  for (int i = 0; i < mrefine; i++)
    mfem_mesh->UniformRefinement();

  int dim  = mfem_mesh->Dimension();
  int sdim = mfem_mesh->SpaceDimension();

  FiniteElementCollection *fec = new LG_FECollection(order, dim);
  FiniteElementSpace *fes = new FiniteElementSpace(mfem_mesh, fec, sdim);


  GridFunction gf(fes);
  VectorFunctionCoefficient E(sdim, udf_function);

  gf.ProjectCoefficient(E);

  stringstream ss;
  ss << "a4_projection_order_" << order << ".vtk";
  ofstream ofs;
  ofs.open(ss.str().c_str(), ofstream::out);
  mfem_mesh->PrintVTK(ofs, vrefine);
  gf.SaveVTK(ofs, "field", vrefine);
  ofs.close();

  delete fes;
  delete fec;
  delete mfem_mesh;
  MPI_Finalize();

  return 0;
}

void udf_function(const Vector &x, Vector &E)
{
  E(0) = 100. * x(0) * x(0);
  E(1) =  50. * x(0) * x(1);
}

Mesh* read_mfem_mesh(const char* mesh_file)
{
  // read the mesh and solution files
  named_ifgzstream meshin(mesh_file);
  if (!meshin)
  {
    cerr << "Can not open mesh file " << mesh_file << ". Exit.\n";
    exit(1);
  }

  Mesh* mesh = new Mesh(meshin, 1, 0, false);
  return mesh;
}

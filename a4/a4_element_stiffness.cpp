//                                Finite Element Assignment 4
//
// Sample runs:  ./a4_element_stiffness --mesh ./data/1x1_square_quad.mesh --order 1
//               ./a4_element_stiffness --mesh ./data/1x1_square_quad.mesh --order 2
//
// Description:  For this example you will need to add your code in the sections marked TODO
//
//               In this example you will use MFEM APIs to compute the following:
//
//               1) the element stiffness matrix associated with the diffusive term
//               of the Galerkin weak form (as part of task 1 of the assignment)
//               2) the element stiffness matrix associated with the advective term
//               of the Galerkin weak from (as part of task 2 of the assignment)
//
//
//               As before this example will use your implementation of the shape
//               functions and their derivatives in the header file LagrangeElements.hpp


#include <mfem.hpp>
#include <fstream>
#include <iostream>

#include "LagrangeElements.hpp"

using namespace std;
using namespace mfem;


Mesh* read_mfem_mesh(const char* mesh_file);

int main(int argc, char *argv[])
{
  int num_procs, myid;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  const char *mfem_mesh_file = "./data/1x1_square_quad.mesh";
  int order  = 1;

  OptionsParser args(argc, argv);
  args.AddOption(&mfem_mesh_file, "-m", "--mesh",
      "MFEM Mesh file to use.");
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

  // Read the mesh
  Mesh* mfem_mesh = read_mfem_mesh(mfem_mesh_file);

  int dim  = mfem_mesh->Dimension();
  int sdim = mfem_mesh->SpaceDimension();

  // Construct the Lagrange finite element collection and space
  FiniteElementCollection *fec = new LG_FECollection(order, dim);
  FiniteElementSpace *fes = new FiniteElementSpace(mfem_mesh, fec, sdim);

  // first get the FiniteElement and ElementTransformation objects since
  // we will need this to compute things like shape functions, their derivatives, and
  // Jacobians
  int i = 0; // we will do this for the first element in the mesh
  ElementTransformation* et = mfem_mesh->GetElementTransformation(i);
  const FiniteElement*   fe = fes->GetFE(i);
  int dofs  = order == 1 ? 4 : 9;


  // As part of Task 1: Compute the element stiffness matrix corresponding to the diffusion term.
  // For this part you can assume that kappa = 1.
  //
  // Hints: Some of the APIs needed to complete this part are
  //
  // -- CalcDShape
  // -- InverseJacobian
  // -- Transpose
  //
  // For more info on how to use these APIs  earch mfem documentation here
  // http://mfem.github.io/doxygen/html/index.html

  double kappa = 1.;
  DenseMatrix elmat_diff;
  elmat_diff.SetSize(dofs);

  const IntegrationRule& ir = IntRules.Get(fe->GetGeomType(), 2*order);

  elmat_diff = 0.;
  for (int j = 0; j < ir.GetNPoints(); j++) { // loop over quadrature points

    const IntegrationPoint& ip = ir.IntPoint(j);
    et->SetIntPoint(&ip);

    // TODO: add the code for computation of element matrix corresponding to the diffusive term
    //       and store it into "elmat_diff"

  }

  printf("element stiffness corresponding to the diffusive term\n");
  elmat_diff.Print();
  printf("\n");


  // As part of Task 2: Compute the element stiffness matrix corresponding to the advective term.
  // For this part you can assume that a = [1,1]
  //
  // Hints: Some of the APIs needed to complete this part are
  //
  // -- CalcShape
  // -- CalcDShape
  // -- InverseJacobian
  // -- Transpose
  //
  // For more info on how to use these APIs  earch mfem documentation here
  // http://mfem.github.io/doxygen/html/index.html

  Vector a(2);
  a(0) = 1.;
  a(1) = 1.;
  DenseMatrix elmat_adv;
  elmat_adv.SetSize(dofs);

  elmat_adv = 0.;
  for (int j = 0; j < ir.GetNPoints(); j++) { // loop over quadrature points

    const IntegrationPoint& ip = ir.IntPoint(j);
    et->SetIntPoint(&ip);

    // TODO: add the code for computation of element matrix corresponding to the advective term
    //       and store it into "elmat_adv"

  }

  printf("element stiffness corresponding to the advective term\n");
  elmat_adv.Print();
  printf("\n");


  /* clean ups */
  delete fes;
  delete fec;
  delete mfem_mesh;
  MPI_Finalize();

  return 0;
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

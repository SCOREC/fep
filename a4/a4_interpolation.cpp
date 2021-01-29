//                                Finite Element Assignment 4
//
// Sample runs:  ./a4_interpolation --mesh ./data/1x1_square_quad.mesh --order 1
//               ./a4_interpolation --mesh ./data/1x1_square_quad.mesh --order 2
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
//               c) computes the L2 norm of the error between the exact field
//               and the interpolated field.
//
//               Note that when using an order p finite element space, it is expected
//               that if the function F(x) is a polynomial of order "p" as well, the error
//               should be zero (with in machine precision).


#include <mfem.hpp>
#include <fstream>
#include <iostream>

#include "LagrangeElements.hpp"

using namespace std;
using namespace mfem;


void udf_order1(const Vector &x, Vector &E);
void udf_order2(const Vector &x, Vector &E);

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
      "Refinement level used to refine the mesh after loading it");
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

  for (int i = 0; i < mrefine; i++)
    mfem_mesh->UniformRefinement();


  int dim  = mfem_mesh->Dimension();
  int sdim = mfem_mesh->SpaceDimension();

  // Construct the Lagrange finite element collection and space
  FiniteElementCollection *fec = new LG_FECollection(order, dim);
  FiniteElementSpace *fes = new FiniteElementSpace(mfem_mesh, fec, sdim);


  // Create a GridFunction and project the exact field on to it
  void (*f_ptr) (const Vector&, Vector&) = order == 1 ? udf_order1 : udf_order2;
  GridFunction gf(fes);
  VectorFunctionCoefficient E(sdim, f_ptr);
  gf.ProjectCoefficient(E);


  // Create a linear (order = 1) H1 finite element collection and space
  // for the coordinate field, since the mesh is linear
  FiniteElementCollection *fec_linear = new H1_FECollection(1, dim);
  FiniteElementSpace *fes_linear = new FiniteElementSpace(mfem_mesh, fec_linear, sdim);
  GridFunction coords(fes_linear);
  mfem_mesh->GetNodes(coords);


  // Computing the L2-norm of the error as defined by
  // sqrt(Sum over all elements of Integral_Omega_e |f_interplated - f_exact|^2 dOmega_e)
  double total_error = 0.;
  for (int i = 0; i < mfem_mesh->GetNE(); i++) { // loop over elements

    ElementTransformation* et = mfem_mesh->GetElementTransformation(i);
    const FiniteElement*   fe = fes->GetFE(i);

    const IntegrationRule& ir = IntRules.Get(fe->GetGeomType(), fe->GetOrder());

    double elem_error = 0.;
    for (int j = 0; j < ir.GetNPoints(); j++) { // loop over quadrature points
      const IntegrationPoint& ip = ir.IntPoint(j);
      et->SetIntPoint(&ip);

      // get f_exact first
      Vector x;
      coords.GetVectorValue(*et, ip, x);
      Vector f_exact(2);
      f_ptr(x, f_exact);

      // get f_interpolated
      Vector f_interpolated;
      gf.GetVectorValue(*et, ip, f_interpolated);

      // compute the normal of the difference
      Vector diff = f_exact;
      diff -= f_interpolated;
      elem_error +=
      	(diff*diff) * /*function*/
      	ip.weight   * /*quadrature point weight*/
      	et->Weight()  /*Jacobian Det of the element*/;
    }
    total_error += elem_error;
  }
  if (myid == 0)
    printf("L2 norm of the interpolation error is %e \n", sqrt(total_error));

  /* clean ups */
  delete fes;
  delete fec;
  delete fes_linear;
  delete fec_linear;
  delete mfem_mesh;
  MPI_Finalize();

  return 0;
}

void udf_order1(const Vector &x, Vector &E)
{
  E(0) = 100. * x(1);
  E(1) =  50. * x(0);
}

void udf_order2(const Vector &x, Vector &E)
{
  E(0) = 100. * x(1) * x(0);
  E(1) =  50. * x(0) * x(0);
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

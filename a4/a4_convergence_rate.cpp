//                                Finite Element Assignment 4
//
// Sample runs:  ./a4_convergence_rate --mesh ./data/1x1_square_quad.mesh --order 1 --mrefine 1
//               ./a4_convergence_rate --mesh ./data/1x1_square_quad.mesh --order 1 --mrefine 2
//               ./a4_convergence_rate --mesh ./data/1x1_square_quad.mesh --order 1 --mrefine 3
//               ./a4_convergence_rate --mesh ./data/1x1_square_quad.mesh --order 1 --mrefine 4
//               ./a4_convergence_rate --mesh ./data/1x1_square_quad.mesh --order 1 --mrefine 5
//               ./a4_convergence_rate --mesh ./data/1x1_square_quad.mesh --order 2 --mrefine 1
//               ./a4_convergence_rate --mesh ./data/1x1_square_quad.mesh --order 2 --mrefine 2
//               ./a4_convergence_rate --mesh ./data/1x1_square_quad.mesh --order 2 --mrefine 3
//               ./a4_convergence_rate --mesh ./data/1x1_square_quad.mesh --order 2 --mrefine 4
//               ./a4_convergence_rate --mesh ./data/1x1_square_quad.mesh --order 2 --mrefine 5
//
// Description:  For this example you will need to add your code in the sections marked TODO
//
//               In this example we will set up and solve a diffusion equation with a know
//               exact solution and study the rate of convergence of the fem solutions by using
//               successively refined meshes. The initial mesh has size h = 1./8.
//
//               when you run the code with option --mrefine = 0 => h = 1/8
//               when you run the code with option --mrefine = 1 => h = 1/16
//               when you run the code with option --mrefine = 2 => h = 1/32
//               and so on.
//
//               The equation and the domains are
//
//               -div(kappa.grad(u))= f on a square domain [0,1]x[0,1]
//
//               The boundary conditions are all Dirichlet with value of 0
//
//               u(0,y) = 0
//               u(1,y) = 0
//               u(x,0) = 0
//               u(x,1) = 0
//
//               We will use kappa = 1. for this example, and the following forcing term
//
//               f(x,y) = 2 . sin(4.pi.x) . sin(4.pi.y)
//
//               The above problem has a exact solution
//
//               u(x,y) = sin(4.pi.x) . sin(4.pi.y) / (4.pi)^2
//
//               which we will use to analyse the rate of convergence for finite element
//               spaces of order p=1 and p=2
//
//
//               Notes on boundary conditions:
//               In MFEM boundary conditions can be set by using entity attributes. The MFEM mesh
//               you are provided with has the following attributes (see ./data/1x1_square_details.png):
//               -- edges on the left boundary have attribute 15
//               -- edges on the right boundary have attribute 7
//               -- edges on the top boundary have attribute 17
//               -- edges on the bottom boundary have attribute 13
//
//               As before this example will use your implementation of the shape
//               functions and their derivatives in the header file LagrangeElements.hpp

#include <mfem.hpp>
#include <fstream>
#include <iostream>

#include "LagrangeElements.hpp"

using namespace std;
using namespace mfem;


const double pi = 3.14159265359;

Mesh* read_mfem_mesh(const char* mesh_file);
double f_exact(const Vector& x);
double force_function(const Vector& x);

int main(int argc, char *argv[])
{
  int num_procs, myid;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  const char *mfem_mesh_file = "./data/1x1_square_quad.mesh";
  int order  = 1;
  int mrefine = 0;


  OptionsParser args(argc, argv);
  args.AddOption(&mfem_mesh_file, "-m", "--mesh",
      "MFEM Mesh file to use.");
  args.AddOption(&order, "-o", "--order",
      "Order for Lagrange Elements 1 or 2 (default 1)");
  args.AddOption(&mrefine, "-mr", "--mrefine",
      "Refinement level to get a finer mesh (defaul 0, no refinement)");
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

  // Load the mesh
  Mesh* serial_mesh = read_mfem_mesh(mfem_mesh_file);

  for (int i = 0; i < mrefine; i++) {
    serial_mesh->UniformRefinement();
  }

  int dim  = serial_mesh->Dimension();
  int sdim = serial_mesh->SpaceDimension();

  ParMesh mfem_mesh(MPI_COMM_WORLD, *serial_mesh);

  // Create the Lagrange finite element collection and space for a scalar Temperature field
  FiniteElementCollection *fec = new LG_FECollection(order, dim);
  ParFiniteElementSpace *fes = new ParFiniteElementSpace(&mfem_mesh, fec, 1);

  // Define the solution vector and initialize to 0. This will hold the fem solution
  // Note that this also sets the zero Dirichlet BCs, so there is no need to set them again
  ParGridFunction u_fem(fes);
  u_fem = 0.;


  // TODO use the code in a4_diffusion.cpp as an example and set the diffusion
  // problem here again. This time with all around Dirichlet boundary u = 0
  //
  // Note that you will need to add:
  // - the linear form (corresponding to the forcing term)
  // - the bilinear form (corresponding to the diffusive term)
  // - and specify the essential boundary conditions
  ParLinearForm b(fes);
  ParBilinearForm a(fes);
  Array<int> ess_tdof_list;






  // Form the linear system and solve
  OperatorPtr A;
  Vector B, X;
  a.FormLinearSystem(ess_tdof_list, u_fem, b, A, X, B);

  // solve step
#ifdef MFEM_USE_SUPERLU // if MFEM has super lu use a direct solve
  SuperLUSolver *superlu = new SuperLUSolver(MPI_COMM_WORLD);
  Operator *SLU_A = new SuperLURowLocMatrix(*A.As<HypreParMatrix>());
  superlu->SetPrintStatistics(true);
  superlu->SetSymmetricPattern(false);
  superlu->SetColumnPermutation(superlu::METIS_AT_PLUS_A);
  superlu->SetRowPermutation(superlu::LargeDiag_MC64);
  superlu->SetIterativeRefine(superlu::SLU_DOUBLE);
  superlu->SetOperator(*SLU_A);
  superlu->Mult(B, X);
#else // if not use an iterative solve
  HypreBoomerAMG amg(*A.As<HypreParMatrix>());
  amg.SetSystemsOptions(dim);
  GMRESSolver gmres(MPI_COMM_WORLD);
  gmres.SetKDim(150);
  gmres.SetRelTol(1e-12);
  gmres.SetMaxIter(500);
  gmres.SetPrintLevel(1);
  gmres.SetOperator(*A.As<HypreParMatrix>());
  gmres.SetPreconditioner(amg);
  gmres.Mult(B, X);
#endif

  // recover the solution
  a.RecoverFEMSolution(X, b, u_fem);


  // Compare the fem solution and the exact solution
  FunctionCoefficient exact_sol(f_exact);
  double l2_error = u_fem.ComputeL2Error(exact_sol);
  if (order == 1)
  {
    if (myid == 0)
      printf("\nL2 norm of the error between linear fem solve and exact solution || u_h - u ||_{L^2} = %e\n",
	  l2_error);
  }
  else
  {
    if (myid == 0)
      printf("\nL2 norm of the error between quadratic fem solve and exact solution || u_h - u ||_{L^2} = %e\n",
	  l2_error);
  }

  // Write to VTK for visualization
  stringstream ss;
  ss << "a4_convergence_rate_order_" << order <<
                         "_mrefine_" << mrefine << ".vtk";

  ofstream ofs;
  ofs.open(ss.str().c_str(), ofstream::out);
  mfem_mesh.PrintVTK(ofs, 1);
  u_fem.SaveVTK(ofs, "field", 1);
  ofs.close();

  delete serial_mesh;
  delete fes;
  delete fec;

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


double f_exact(const Vector& x)
{
  return sin(4.*pi*x(0))*sin(4.*pi*x(1))/4./4./pi/pi;
}

double force_function(const Vector& x)
{
  return 2.*sin(4.*pi*x(0))*sin(4.*pi*x(1));
}

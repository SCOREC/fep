//                                Finite Element Assignment 4
//
// Sample runs:  ./a4_advection_diffusion --mesh ./data/1x1_square_quad.mesh --order 1 --no-stabilization
//               ./a4_advection_diffusion --mesh ./data/1x1_square_quad.mesh --order 1 --stabilization
//               ./a4_advection_diffusion --mesh ./data/1x1_square_quad.mesh --order 2 --no-stabilization
//               ./a4_advection_diffusion --mesh ./data/1x1_square_quad.mesh --order 2 --stabilization
//
//               ./a4_advection_diffusion --mesh ./data/1x1_square_quad_2by8.mesh --order 1 --no-stabilization
//               ./a4_advection_diffusion --mesh ./data/1x1_square_quad_2by8.mesh --order 1 --stabilization
//               ./a4_advection_diffusion --mesh ./data/1x1_square_quad_2by8.mesh --order 2 --no-stabilization
//               ./a4_advection_diffusion --mesh ./data/1x1_square_quad_2by8.mesh --order 2 --stabilization
//
// Description:  For this example you will need to add your code in the sections marked TODO
//
//               In this example we will set up and solve a simple advection diffusion equation
//               w.a.grad(u) - div(kappa.grad(u))= f on a square domain [0,1]x[0,1] and the following BCs
//
//               u(0,y) = 0
//               u(1,y) = 0
//               du/dy(x,0) = 0 (this is a zero-flux condition)
//               du/dy(x,1) = 0 (this is a zero-flux condition)
//
//               We will use kappa = 0.01 and the vector a = (1.0, 0.0) for this example. The force term
//               on the right hand side will be set to 0, as well.
//
//               The above boundary conditions are such that we will have a 1D solution, and there for
//               the mesh size h you will need to use the edge length in the x-directions.
//
//
//               We will also add a simple stabilization term to the above equation and study the effect
//               of that in to the finite element solution.
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


Mesh* read_mfem_mesh(const char* mesh_file);
double f_exact(const Vector& x);

int main(int argc, char *argv[])
{
  int num_procs, myid;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  const char *mfem_mesh_file = "./data/1x1_square_quad.mesh";
  int order  = 1;
  bool stabilization = false;


  OptionsParser args(argc, argv);
  args.AddOption(&mfem_mesh_file, "-m", "--mesh",
      "MFEM Mesh file to use.");
  args.AddOption(&order, "-o", "--order",
      "Order for Lagrange Elements 1 or 2 (default 1)");
  args.AddOption(&stabilization, "-stab", "--stabilization", "-no_stab", "--no-stabilization",
      "Whether to add the stabilization term (default does not add it)");
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

  int dim  = serial_mesh->Dimension();
  int sdim = serial_mesh->SpaceDimension();

  ParMesh mfem_mesh(MPI_COMM_WORLD, *serial_mesh);

  // Create the Lagrange finite element collection and space for a scalar Temperature field
  FiniteElementCollection *fec = new LG_FECollection(order, dim);
  ParFiniteElementSpace *fes = new ParFiniteElementSpace(&mfem_mesh, fec, 1);

  // Define the solution vector and initialize to 0. This will hold the fem solution
  ParGridFunction u_fem(fes);
  u_fem = 0.;


  // TODO Setup and solve the advection diffusion problem here
  // for both standard Galerkin from and stabilized versions
  //
  // The advection term can be added by using ConvectionIntegrator (refer to mfem documentation
  // here http://mfem.github.io/doxygen/html/index.html for more details)
  //
  // Note that for the version with stabilization you will need to define a Matrix of diffusivity.
  // That can be done by using
  // MatrixCoefficient* kappa_tau = new MatrixConstantCoefficient(k_tau);
  ParLinearForm b(fes);
  ParBilinearForm a(fes);
  Array<int> ess_tdof_list;

  if (stabilization)
  {
    // TODO the stabilization term should be added in this conditional.
  }

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
  gmres.SetKDim(100);
  gmres.SetRelTol(1e-12);
  gmres.SetMaxIter(500);
  gmres.SetPrintLevel(1);
  gmres.SetOperator(*A.As<HypreParMatrix>());
  gmres.SetPreconditioner(amg);
  gmres.Mult(B, X);
  /* GSSmoother M((SparseMatrix&)(*A)); */
  /* GMRES(*A, M, B, X, 0, 200, 50, 1e-24, 0.0); */
#endif

  // recover the solution
  a.RecoverFEMSolution(X, b, u_fem);


  // Compare the fem solution and the exact solution
  if (order == 1)
  {
    FunctionCoefficient exact_sol(f_exact);
    printf("\nL2 norm of the error between linear fem solve and exact solution || u_h - u ||_{L^2} = %e\n",
    	u_fem.ComputeL2Error(exact_sol));
  }
  else
  {
    FunctionCoefficient exact_sol(f_exact);
    printf("\nL2 norm of the error between quadratic fem solve and exact solution || u_h - u ||_{L^2} = %e\n",
    	u_fem.ComputeL2Error(exact_sol));
  }

  // Write to VTK for visualization
  stringstream ss;
  ss << "a4_advection_diffusion_order_" << order;
  if (stabilization)
    ss << "_with_stabilization";
  else
    ss << "_without_stabilization";
  ss << ".vtk";

  ofstream ofs;
  ofs.open(ss.str().c_str(), ofstream::out);
  int vrefine = 1;
  if (order == 2)
    vrefine = 8;
  mfem_mesh.PrintVTK(ofs, vrefine);
  u_fem.SaveVTK(ofs, "field", vrefine);
  ofs.close();

  delete serial_mesh;
  delete fes;
  delete fec;
  /* MPI_Finalize(); */

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
  // TODO return the exact solution here for the companion of the error
  return 0.;
}

//                                Finite Element Assignment 4
//
// Sample runs:  ./a4_diffusion --mesh ./data/1x1_square_quad.mesh --order 1
//               ./a4_diffusion	--mesh ./data/1x1_square_quad.mesh --order 2
//
// Description:  For this example you will need to add your code in the sections marked TODO
//
//               In this example we will set up and solve a simple diffusion equation
//               -div(kappa.grad(u))= f on a square domain [0,1]x[0,1] and the following BCs
//
//               u(0,y) = 0
//               u(1,y) = 0
//               du/dy(x,0) = 0 (this is a zero-flux condition)
//               du/dy(x,1) = 0 (this is a zero-flux condition)
//
//               We will use kappa = 1. for this example.
//
//               The above boundary conditions are such that we will have a 1D "exact" solutions
//               (i.e., u(x,y) = u(x)).
//
//               -- For order 1 we will use f = 0 (you will have to workout the exact solution
//                  yourself and add it to f_exact_linear)
//
//               -- For order 2 we will use f = 1 (you will have to workout the exact solution
//                  yourself and add it to f_exact_quadratic)
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
double f_exact_linear(const Vector& x);
double f_exact_quadratic(const Vector& x);

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
      "Order for Lagrange Elements 1 or 2 (default 1)");
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

  // Setting the Dirichlet boundary conditions
  // 1- left boundary: u = 0 and since the grid function u_fem
  //    is initialized to 0 nothing needs to be
  // 2- right boundary: u = 1
  Array<int> bdr_mask_dirichlet(mfem_mesh.bdr_attributes.Max());
  bdr_mask_dirichlet = 0;
  for (int i = 0; i < bdr_mask_dirichlet.Size(); i++) {
    if (i+1 == 7) // right edge with boundary attribute 7
      bdr_mask_dirichlet[i] = 1;
  }
  ConstantCoefficient one(1.);
  u_fem.ProjectBdrCoefficient(one, bdr_mask_dirichlet);

  // Setting the Neumann boundary conditions
  // top and bottom: zero flux
  //
  // Neumann boundary conditions go to the right hand side of the weak form in an integral form,
  // there for they have to be added as part the assembly. The will be done further below but for
  // now we should create the bdr_mask_neumann to tell mfem which boundaries have this conditions
  Array<int> bdr_mask_neumann(mfem_mesh.bdr_attributes.Max());
  bdr_mask_neumann = 0;
  for (int i = 0; i < bdr_mask_neumann.Size(); i++) {
    if (i+1 == 13) // bottom edge with boundary attribute 13
      bdr_mask_neumann[i] = 1;
    else if (i+1 == 17) // top edge with boundary attribute 17
      bdr_mask_neumann[i] = 1;
    else
      bdr_mask_neumann[i] = 0;
  }

  // Add the right-hand side corresponding to the force term and the zero-flux conditions
  // in mfem terminology this is known as the LinearForm
  ParLinearForm b(fes);
  // force will be zero for order 1 and will be 1 for order 2 (this is because for this we
  // know the exact solutions)
  ConstantCoefficient force(order == 1 ? 0 : 1);
  b.AddDomainIntegrator(new DomainLFIntegrator(force));
  // add the zero flux term to the right hand side
  ConstantCoefficient kappa(1.);
  ConstantCoefficient flux(0.);
  ProductCoefficient kappa_flux(kappa, flux);
  b.AddBoundaryIntegrator(new BoundaryLFIntegrator(kappa_flux), bdr_mask_neumann);
  b.Assemble();


  // Add the right-hand-side corresponding to the diffusive term
  // in mfem terminology this is known as the Bilinear form
  ParBilinearForm a(fes);
  a.AddDomainIntegrator(new DiffusionIntegrator(kappa));
  a.Assemble();

  // Get the essential dofs
  // in this case that would be boundaries with Dirichlet (i.e., left and right boundaries)
  Array<int> ess_tdof_list;
  Array<int> ess_bdr(mfem_mesh.bdr_attributes.Max());
  ess_bdr = 0;
  for (int i = 0; i < ess_bdr.Size(); i++) {
    if (i+1 == 15) // left edge with boundary attribute 15
      ess_bdr[i] = 1;
    else if (i+1 == 7) // right edge with boundary attribute 7
      ess_bdr[i] = 1;
    else
      ess_bdr[i] = 0;
  }
  fes->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);



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
    FunctionCoefficient exact_sol(f_exact_linear);
    double l2_error = u_fem.ComputeL2Error(exact_sol);
    printf("\nL2 norm of the error between linear fem solve and exact solution || u_h - u ||_{L^2} = %e\n",
    	l2_error);
  }
  else
  {
    FunctionCoefficient exact_sol(f_exact_quadratic);
    double l2_error = u_fem.ComputeL2Error(exact_sol);
    printf("\nL2 norm of the error between quadratic fem solve and exact solution || u_h - u ||_{L^2} = %e\n",
    	l2_error);
  }

  // Write to VTK for visualization
  stringstream ss;
  ss << "a4_diffusion_order_" << order << ".vtk";

  ofstream ofs;
  ofs.open(ss.str().c_str(), ofstream::out);
  mfem_mesh.PrintVTK(ofs, 1);
  u_fem.SaveVTK(ofs, "field", 1);
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


double f_exact_linear(const Vector& x)
{
  // TODO add the exact solution value as a function of x(i) for the linear (p=1) case
  return 0.;
}


double f_exact_quadratic(const Vector& x)
{
  // TODO add the exact solution value as a function of x(i) for the quadratic (p=2) case
  return 0.;
}

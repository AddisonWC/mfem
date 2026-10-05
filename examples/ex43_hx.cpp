//                             MFEM Example 43 HX
//
// Compile with: make ex43_hx
// Sample runs: ex43_hx -gmg -r 3 -random-rhs
//              ex43_hx -gmg -e aw -r 3 -random-rhs
//              ex43_hx -gmg -e hzzzr -r 3 -random-rhs
//
// Description: Solve a div-div plus mass problem for a symmetric matrix
// field using lowest-order 2D Johnson--Mercier, Arnold--Winther, Hu--Zhang, or
// (reduced) Huang--Zhang--Zhou--Zhu elements
// and the auxiliary space preconditioner
//
//                 R + Pi B_1 Pi^t + J B_2 J^t.
//
// Here R is a vertex-patch, Jacobi, or symmetric Gauss-Seidel smoother.
// Pi is the canonical interpolant
// from continuous piecewise-linear symmetric matrices, B_1 is the inverse of
// the matrix H1 operator, J is the Airy map from the HCT, Argyris, or Bell space, and B_2 is the
// inverse of the corresponding biharmonic operator.
// With -gmg, B_1 and B_2 are replaced by geometric multigrid V-cycles with
// Gauss-Seidel smoothing and rediscretized operators on each mesh level.

#include "ex43.hpp"
#include <iostream>
#include <memory>
#include <string>

using namespace mfem;
using namespace std;

enum class HXSmoother { VERTEX_PATCH, JACOBI, GAUSS_SEIDEL };

static unique_ptr<Solver> MakeInverse(SparseMatrix &op)
{
#ifdef MFEM_USE_SUITESPARSE
   return unique_ptr<Solver>(new UMFPackSolver(op));
#else
   CGSolver *inverse = new CGSolver;
   inverse->iterative_mode = false;
   inverse->SetOperator(op);
   inverse->SetRelTol(1e-12);
   inverse->SetAbsTol(0.0);
   inverse->SetMaxIter(20000);
   inverse->SetPrintLevel(0);
   return unique_ptr<Solver>(inverse);
#endif
}

static void GetPotentialGaugeDofs(FiniteElementSpace &fes, Array<int> &dofs)
{
   // Hessians annihilate affine functions. Fix value and both first
   // derivatives at one vertex to select a representative modulo P1.
   fes.GetVertexDofs(0, dofs);
   // Argyris and Bell also have three second derivatives at each vertex; those
   // are not in the affine kernel and must remain unconstrained.
   dofs.SetSize(3);
   for (int i = 0; i < dofs.Size(); i++) { dofs[i] = UnsignIndex(dofs[i]); }
}

// Assemble either the weighted matrix H1 operator or, with a null weight,
// the biharmonic operator. Both use a forward GS pre-sweep and its transpose
// for the post-sweep, preserving symmetry for the outer CG solve.
class HXMultigrid : public GeometricMultigrid
{
public:
   HXMultigrid(FiniteElementSpaceHierarchy &hierarchy, MatrixCoefficient *weight)
      : GeometricMultigrid(hierarchy, Array<int>())
   {
      for (int level = 0; level < hierarchy.GetNumLevels(); level++)
      {
         FiniteElementSpace &fes = hierarchy.GetFESpaceAtLevel(level);
         BilinearForm *form = new BilinearForm(&fes);
         if (weight)
         {
            form->AddDomainIntegrator(new VectorMassIntegrator(*weight));
            form->AddDomainIntegrator(new VectorDiffusionIntegrator(*weight));
         }
         else
         {
            form->AddDomainIntegrator(new HessianIntegrator);
            GetPotentialGaugeDofs(fes, *essentialTrueDofs[level]);
            if (level > 0)
            {
               // Preserve the point gauge under both prolongation and
               // restriction, without imposing clamped boundary conditions.
               prolongations[level - 1] = new RectangularConstrainedOperator(
                  hierarchy.GetProlongationAtLevel(level - 1),
                  *essentialTrueDofs[level - 1], *essentialTrueDofs[level]);
               ownedProlongations[level - 1] = true;
            }
         }
         form->SetDiagonalPolicy(Operator::DIAG_ONE);
         form->Assemble();
         bfs.Append(form);

         OperatorPtr system(Operator::MFEM_SPARSEMAT);
         form->FormSystemMatrix(*essentialTrueDofs[level], system);
         system.SetOperatorOwner(false);
         Solver *level_solver = level == 0 ?
                                MakeInverse(*system.As<SparseMatrix>()).release() :
                                new GSSmoother(*system.As<SparseMatrix>(),
                                               GSSmoother::FORWARD, 1);
         AddLevel(system.Ptr(), level_solver, false, true);
      }
   }
};

class HXPreconditioner : public Solver
{
private:
   FiniteElementSpace &h1_fespace;
   FiniteElementSpace &potential_fespace;
   MatrixConstantCoefficient matrix_h1_coefficient;
   BilinearForm matrix_h1_form;
   BilinearForm biharmonic_form;
   DiscreteLinearOperator pi;
   DiscreteLinearOperator airy;
   unique_ptr<Solver> smoother;
   Array<int> potential_gauge_dofs;
   unique_ptr<Solver> matrix_h1_inverse;
   unique_ptr<Solver> biharmonic_inverse;
   mutable Vector matrix_h1_rhs, matrix_h1_solution;
   mutable Vector biharmonic_rhs, biharmonic_solution;
   mutable Vector auxiliary_correction;

   static DenseMatrix MatrixH1Weight()
   {
      DenseMatrix weight(3);
      weight = 0.0;
      weight(0,0) = 1.0;
      weight(1,1) = 2.0;
      weight(2,2) = 1.0;
      return weight;
   }

public:
   HXPreconditioner(const SparseMatrix &op, FiniteElementSpace &stress_fespace,
                    FiniteElementSpaceHierarchy &h1_hierarchy,
                    FiniteElementSpaceHierarchy &potential_hierarchy,
                    HXSmoother smoother_type, bool use_gmg)
      : Solver(op.Height()),
        h1_fespace(h1_hierarchy.GetFinestFESpace()),
        potential_fespace(potential_hierarchy.GetFinestFESpace()),
        matrix_h1_coefficient(MatrixH1Weight()),
        matrix_h1_form(&h1_fespace),
        biharmonic_form(&potential_fespace),
        pi(&h1_fespace, &stress_fespace),
        airy(&potential_fespace, &stress_fespace)
   {
      MFEM_VERIFY(stress_fespace.GetTrueVSize() == stress_fespace.GetVSize(),
                  "HXPreconditioner currently requires a conforming mesh");

      switch (smoother_type)
      {
         case HXSmoother::VERTEX_PATCH:
            smoother.reset(new VertexPatchSmoother(op, stress_fespace));
            break;
         case HXSmoother::JACOBI:
            smoother.reset(new DSmoother(op, DSmoother::JACOBI));
            break;
         case HXSmoother::GAUSS_SEIDEL:
            // Forward/backward sweeps keep the HX preconditioner symmetric
            // for the outer conjugate-gradient solver.
            smoother.reset(new GSSmoother(op, GSSmoother::SYMMETRIC));
            break;
      }

      pi.AddDomainInterpolator(new IdentityInterpolator);
      pi.Assemble();
      pi.Finalize();

      airy.AddDomainInterpolator(new AiryInterpolator);
      airy.Assemble();
      airy.Finalize();

      GetPotentialGaugeDofs(potential_fespace, potential_gauge_dofs);
      if (use_gmg)
      {
         matrix_h1_inverse.reset(new HXMultigrid(h1_hierarchy,
                                                 &matrix_h1_coefficient));
         biharmonic_inverse.reset(new HXMultigrid(potential_hierarchy, nullptr));
      }
      else
      {
         matrix_h1_form.AddDomainIntegrator(
            new VectorMassIntegrator(matrix_h1_coefficient));
         matrix_h1_form.AddDomainIntegrator(
            new VectorDiffusionIntegrator(matrix_h1_coefficient));
         matrix_h1_form.Assemble();
         matrix_h1_form.Finalize();
         matrix_h1_inverse = MakeInverse(matrix_h1_form.SpMat());

         biharmonic_form.AddDomainIntegrator(new HessianIntegrator);
         biharmonic_form.Assemble();
         biharmonic_form.Finalize();

         for (int i = 0; i < potential_gauge_dofs.Size(); i++)
         {
            biharmonic_form.SpMat().EliminateRowCol(potential_gauge_dofs[i]);
         }
         biharmonic_inverse = MakeInverse(biharmonic_form.SpMat());
      }
   }

   void Mult(const Vector &x, Vector &y) const override
   {
      smoother->Mult(x, y);

      matrix_h1_rhs.SetSize(pi.Width());
      pi.MultTranspose(x, matrix_h1_rhs);
      matrix_h1_solution.SetSize(matrix_h1_rhs.Size());
      matrix_h1_inverse->Mult(matrix_h1_rhs, matrix_h1_solution);
      auxiliary_correction.SetSize(pi.Height());
      pi.Mult(matrix_h1_solution, auxiliary_correction);
      y += auxiliary_correction;

      biharmonic_rhs.SetSize(airy.Width());
      airy.MultTranspose(x, biharmonic_rhs);
      for (int i = 0; i < potential_gauge_dofs.Size(); i++)
      {
         biharmonic_rhs(potential_gauge_dofs[i]) = 0.0;
      }
      biharmonic_solution.SetSize(biharmonic_rhs.Size());
      biharmonic_inverse->Mult(biharmonic_rhs, biharmonic_solution);
      auxiliary_correction.SetSize(airy.Height());
      airy.Mult(biharmonic_solution, auxiliary_correction);
      y += auxiliary_correction;
   }

   void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }
   void SetOperator(const Operator &) override
   { MFEM_ABORT("HXPreconditioner does not support SetOperator"); }
};

int main(int argc, char *argv[])
{
   const char *mesh_file = "../data/ref-triangle.mesh";
   int refinements = 2;
   std::string element_name = "jm";
   std::string smoother_name = "vertex-patch";
   bool visualization = false;
   bool random_rhs = false;
   bool use_gmg = false;
   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Input triangle mesh.");
   args.AddOption(&refinements, "-r", "--refinements",
                  "Number of uniform refinements.");
   args.AddOption(&use_gmg, "-gmg", "--geometric-multigrid",
                  "-no-gmg", "--no-geometric-multigrid",
                  "Use a geometric multigrid V-cycle with Gauss-Seidel smoothing "
                  "for each global auxiliary solve.");
   args.AddOption(&visualization, "-vis", "--visualization",
                  "-no-vis", "--no-visualization",
                  "Enable or disable visualization (accepted for consistency).");
   args.AddOption(&random_rhs, "-random-rhs", "--random-rhs",
                  "-constant-rhs", "--constant-rhs",
                  "Use a reproducible random algebraic right-hand side.");
   args.AddOption(&element_name, "-e", "--element",
                  "Stress element: jm, aw, hz, hzzz, or hzzzr (reduced HZZZ).");
   args.AddOption(&smoother_name, "-s", "--smoother",
                  "HX smoother: vertex-patch (default), jacobi, or gauss-seidel "
                  "(symmetric forward/backward sweeps).");
   args.ParseCheck();
   MFEM_VERIFY(refinements >= 0, "Refinement count must be nonnegative.");

   const char *fec_name;
   const char *potential_name;
   if (element_name == "jm")
   {
      fec_name = "JM_2D_P1";
      potential_name = "HCT_2D_P3";
   }
   else if (element_name == "aw")
   {
      fec_name = "AW_2D_P3";
      potential_name = "Argyris_2D_P5";
   }
   else if (element_name == "hz")
   {
      fec_name = "HZ_2D_P3";
      potential_name = "Argyris_2D_P5";
   }
   else if (element_name == "hzzz")
   {
      fec_name = "HZZZ_2D_P3";
      potential_name = "Bell_2D_P5";
   }
   else if (element_name == "hzzzr")
   {
      fec_name = "HZZZr_2D_P3";
      potential_name = "Bell_2D_P5";
   }
   else
   {
      MFEM_ABORT("Unknown stress element '" << element_name
                 << "'. Choose jm, aw, hz, hzzz, or hzzzr.");
   }

   HXSmoother smoother_type;
   if (smoother_name == "vertex-patch")
   {
      smoother_type = HXSmoother::VERTEX_PATCH;
   }
   else if (smoother_name == "jacobi")
   {
      smoother_type = HXSmoother::JACOBI;
   }
   else if (smoother_name == "gauss-seidel")
   {
      smoother_type = HXSmoother::GAUSS_SEIDEL;
   }
   else
   {
      MFEM_ABORT("Unknown HX smoother '" << smoother_name
                 << "'. Choose vertex-patch, jacobi, or gauss-seidel.");
   }

   Mesh mesh(mesh_file);
   MFEM_VERIFY(mesh.Dimension() == 2, "");
   if (!use_gmg)
   {
      for (int level = 0; level < refinements; level++)
      {
         mesh.UniformRefinement();
      }
   }

   H1_FECollection h1_fec(1, 2);
   unique_ptr<FiniteElementCollection> potential_fec(
      FiniteElementCollection::New(potential_name));
   FiniteElementSpaceHierarchy h1_hierarchy(
      &mesh, new FiniteElementSpace(&mesh, &h1_fec, 3, Ordering::byVDIM),
      false, true);
   // Both auxiliary hierarchies and the stress space share the same meshes.
   // Refined meshes are owned by h1_hierarchy, which outlives potential_hierarchy.
   FiniteElementSpaceHierarchy potential_hierarchy(
      &mesh, new FiniteElementSpace(&mesh, potential_fec.get()), false, true);
   for (int level = 0; use_gmg && level < refinements; level++)
   {
      h1_hierarchy.AddUniformlyRefinedLevel(3, Ordering::byVDIM,
                                            Operator::MFEM_SPARSEMAT);
      Mesh *fine_mesh = h1_hierarchy.GetFinestFESpace().GetMesh();
      FiniteElementSpace *fine_potential =
         new FiniteElementSpace(fine_mesh, potential_fec.get());
      OperatorPtr transfer(Operator::MFEM_SPARSEMAT);
      // HCT and Bell are nonnested: use DOF interpolation, as in biharmonic_gmg.
      fine_potential->GetTrueTransferOperator(
         potential_hierarchy.GetFinestFESpace(), transfer);
      potential_hierarchy.AddLevel(fine_mesh, fine_potential, transfer.Ptr(),
                                   false, true, true);
      transfer.SetOperatorOwner(false);
   }

   unique_ptr<FiniteElementCollection> fec(FiniteElementCollection::New(fec_name));
   FiniteElementSpace fespace(h1_hierarchy.GetFinestFESpace().GetMesh(),
                              fec.get());
   cout << "\n" << fec->Name() << " space: " << fespace.GetNE()
        << " elements, " << fespace.GetTrueVSize() << " unknowns\n";

   DenseMatrix identity(2);
   identity = 0.0;
   identity(0,0) = identity(1,1) = 1.0;
   MatrixConstantCoefficient rhs(identity);
   LinearForm b(&fespace);
   b.AddDomainIntegrator(new MatrixFEDomainLFIntegrator(rhs));
   b.Assemble();
   GridFunction solution(&fespace);
   solution = 0.0;

   BilinearForm a(&fespace);
   a.AddDomainIntegrator(new MatrixDivDivIntegrator);
   a.AddDomainIntegrator(new MatrixFEMassIntegrator);
   a.Assemble();
   Array<int> ess_tdof_list;
   SparseMatrix A;
   Vector B, X;
   a.FormLinearSystem(ess_tdof_list, solution, b, A, X, B);
   if (random_rhs) { B.Randomize(1); }

   HXPreconditioner hx(A, fespace, h1_hierarchy, potential_hierarchy,
                       smoother_type, use_gmg);
   CGSolver solver;
   solver.SetOperator(A);
   solver.SetPreconditioner(hx);
   solver.SetRelTol(1e-10);
   solver.SetAbsTol(0.0);
   solver.SetMaxIter(500);
   solver.SetPrintLevel(1);
   solver.Mult(B, X);
   a.RecoverFEMSolution(X, b, solution);
   cout << "PCG iterations: " << solver.GetNumIterations() << '\n';
   return 0;
}

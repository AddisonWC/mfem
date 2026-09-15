// JM geometric multigrid with macro patches or split patches plus HCT Jacobi.
// All level operators and transfers stay in the original moment coordinates.
// See reports/mg_compare.md for the experiment and its measured results.
#include "common.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>

using namespace mfem;
using namespace mfem::jm_hx;
using Clock = std::chrono::steady_clock;

class ScaledInverse : public Solver
{
   const Solver &inverse;
   real_t weight;
public:
   ScaledInverse(const Solver &R, real_t w)
      : Solver(R.Height()), inverse(R), weight(w) { }
   void Mult(const Vector &x, Vector &y) const override
   { inverse.Mult(x, y); y *= weight; }
   void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }
   void SetOperator(const Operator &) override { MFEM_ABORT("fixed smoother"); }
};

// Require the same actual moment-coordinate residual in every comparison.
class ResidualMonitor : public IterativeSolverController
{
   const Operator &A;
   const Vector &b;
   real_t threshold;
   Vector residual;
public:
   ResidualMonitor(const Operator &op, const Vector &rhs, real_t tol)
      : A(op), b(rhs), threshold(tol*rhs.Norml2()), residual(rhs.Size()) { }
   void MonitorResidual(int, real_t, const Vector &r, bool) override
   { converged = r.Norml2() <= threshold; }
   bool RequiresUpdatedSolution() const override { return true; }
   void MonitorSolution(int, real_t, const Vector &x, bool) override
   {
      if (converged)
      {
         A.Mult(x, residual); residual -= b;
         converged = residual.Norml2() <= threshold;
      }
   }
};

struct Level
{
   FiniteElementSpace vertices;
   BilinearForm form;
   std::unique_ptr<SparseMatrix> B, vertex_matrix;
   Level(FiniteElementSpace &moments, JohnsonMercierFECollection &fec)
      : vertices(moments.GetMesh(), &fec), form(&moments)
   {
      form.AddDomainIntegrator(new MatrixFEMassIntegrator);
      form.AddDomainIntegrator(new MatrixDivDivIntegrator);
      form.Assemble(); form.Finalize();
      B = BasisMatrix(moments, vertices);
      // Congruence also works for Galerkin operators, if added later.
      vertex_matrix.reset(RAP(*B, form.SpMat(), *B));
   }
};

int main(int argc, char *argv[])
{
   const char *mesh_file = JM_DEFAULT_MESH;
   const char *choice_arg = "all";
   int refinements = 3, coarse_refinements = 0, steps = 1, max_it = 2000;
   real_t airy_damping = 0.05;
   real_t damping = 0.33, tolerance = 1e-8;
   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Straight connected triangle mesh.");
   args.AddOption(&refinements, "-r", "--refinements", "Refinements above the coarse grid.");
   args.AddOption(&coarse_refinements, "-cr", "--coarse-refinements", "Refine before building the hierarchy.");
   args.AddOption(&choice_arg, "-smoother", "--smoother", "all, macro, macro-vertex, split-patches, or split (patches + Airy Jacobi).");
   args.AddOption(&steps, "-s", "--steps", "Equal pre- and post-smoothing counts.");
   args.AddOption(&damping, "-damping", "--damping", "Patch weight (0 < weight < 0.5).");
   args.AddOption(&airy_damping, "-airy-damping", "--airy-damping",
                  "Airy Jacobi weight; split requires 4*damping + 12*airy-damping < 2.");
   args.AddOption(&tolerance, "-tol", "--tolerance", "True relative residual tolerance.");
   args.AddOption(&max_it, "-max-it", "--max-iterations", "PCG iteration limit.");
   args.ParseCheck(std::cerr);
   const std::string choice(choice_arg);
   MFEM_VERIFY(choice == "all" || choice == "macro" || choice == "macro-vertex" ||
               choice == "split" || choice == "split-patches", "invalid smoother");
   MFEM_VERIFY(refinements >= 0 && coarse_refinements >= 0 && steps > 0 &&
               damping > 0 && damping < 0.5 && tolerance > 0 && max_it > 0,
               "invalid experiment parameters");
   MFEM_VERIFY(airy_damping > 0 &&
               ((choice != "all" && choice != "split") ||
                4*damping + 12*airy_damping < 2), "unsafe combined smoother weights");
   JohnsonMercierFECollection moments_fec, vertices_fec(JMBasis::SplitVertex);
   auto mesh = new Mesh(mesh_file);
   { auto check = AlfeldMesh(*mesh); }
   for (int r = 0; r < coarse_refinements; r++) { mesh->UniformRefinement(); }
   auto coarse = new FiniteElementSpace(mesh, &moments_fec);
   FiniteElementSpaceHierarchy hierarchy(mesh, coarse, true, true);
   for (int r = 0; r < refinements; r++)
   { hierarchy.AddUniformlyRefinedLevel(1, Ordering::byVDIM, Operator::MFEM_SPARSEMAT); }
   std::vector<std::unique_ptr<Level>> levels;
   Array<Operator *> operators, transfers;
   for (int l = 0; l <= refinements; l++)
   {
      levels.emplace_back(new Level(hierarchy.GetFESpaceAtLevel(l), vertices_fec));
      operators.Append(&levels.back()->form.SpMat());
      if (l < refinements)
      {
         transfers.Append(hierarchy.GetProlongationAtLevel(l));
      }
   }
   ExactSolver coarse_inverse(levels[0]->form.SpMat());
   bool success = true;
   std::cout << "refinement,elements,dofs,smoother,damping,airy_damping,steps,iterations,residual,converged,setup_s,solve_s\n"
             << std::setprecision(9);
   // Compare every prefix of ONE hierarchy; all configurations borrow the exact
   // same transfer objects. No transfer is constructed in split coordinates.
   for (int depth = 0; depth <= refinements; depth++)
   {
      auto &A = levels[depth]->form.SpMat();
      Vector rhs(A.Height()), x(A.Height()), residual(A.Height());
      rhs.Randomize(1);
      for (const std::string kind : {"macro", "macro-vertex", "split-patches", "split"})
      {
         if (choice != "all" && choice != kind) { continue; }
         const auto start = Clock::now();
         Array<Operator *> ops, prolongations;
         Array<Solver *> smoothers;
         Array<bool> own_ops, own_smoothers, own_prolongations;
         std::vector<std::unique_ptr<Solver>> local, mapped, scaled, auxiliary,
             weighted_auxiliary, combined;
         for (int l = 0; l <= depth; l++)
         {
            ops.Append(operators[l]); own_ops.Append(false);
            own_smoothers.Append(false);
            if (!l) { smoothers.Append(&coarse_inverse); }
            else
            {
               auto &level = *levels[l];
               const bool vertex = kind != "macro";
               local.emplace_back(new PatchSolver(
                                     vertex ? *level.vertex_matrix : level.form.SpMat(),
                                     vertex ? level.vertices : hierarchy.GetFESpaceAtLevel(l),
                                     (kind == "split" || kind == "split-patches")));
               Solver *R = local.back().get();
               if (vertex)
               {
                  mapped.emplace_back(new MappedSolver(*level.B, *R));
                  R = mapped.back().get();
               }
               if (kind == "split")
               {
                  auxiliary.emplace_back(
                     new AiryJacobi(hierarchy.GetFESpaceAtLevel(l)));
                  weighted_auxiliary.emplace_back(
                     new ScaledInverse(*auxiliary.back(), airy_damping));
                  combined.emplace_back(
                     new HX(*R, *weighted_auxiliary.back(), damping));
                  smoothers.Append(combined.back().get());
               }
               else
               {
                  scaled.emplace_back(new ScaledInverse(*R, damping));
                  smoothers.Append(scaled.back().get());
               }
               prolongations.Append(transfers[l-1]); own_prolongations.Append(false);
            }
         }
         Multigrid mg(ops, smoothers, prolongations, own_ops, own_smoothers, own_prolongations);
         mg.SetCycleType(Multigrid::CycleType::VCYCLE, steps, steps);
         const double setup = std::chrono::duration<double>(Clock::now()-start).count();
         CGSolver cg;
         ResidualMonitor monitor(A, rhs, tolerance);
         cg.SetOperator(A); cg.SetPreconditioner(mg); cg.SetController(monitor);
         cg.SetRelTol(0); cg.SetAbsTol(0); cg.SetMaxIter(max_it); cg.SetPrintLevel(-1);
         x = 0.0;
         const auto solve_start = Clock::now();
         cg.Mult(rhs, x);
         const double solve = std::chrono::duration<double>(Clock::now()-solve_start).count();
         A.Mult(x, residual); residual -= rhs;
         const real_t rel = residual.Norml2()/rhs.Norml2();
         const bool converged = cg.GetConverged() && rel <= tolerance;
         success = success && converged;
         std::cout << coarse_refinements+depth << ',' << hierarchy.GetFESpaceAtLevel(depth).GetNE()
                   << ',' << A.Height() << ',' << kind << ',' << damping
                   << ',' << (kind == "split" ? airy_damping : 0.0) << ',' << steps
                   << ',' << cg.GetNumIterations() << ',' << rel << ',' << converged
                   << ',' << setup << ',' << solve << std::endl;
      }
   }
   return success ? 0 : 3;
}

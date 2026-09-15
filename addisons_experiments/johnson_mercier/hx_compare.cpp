// JM H(div) HX experiment: basis, smoother, and auxiliary H1 mesh comparisons.
//
// Solve (sigma,tau) + (div sigma,div tau) = (f,tau) with no essential
// boundary conditions. Every comparison solves in the SAME moment coordinates;
// vertex-basis smoothers are mapped back as B R B^t. The HCT Airy correction is
// always present. Split patches always denote the intrinsic split-vertex
// subspaces, including when the selected basis is the moment basis.
//
// Examples (from addisons_experiments/build):
//   ./hx_compare -r 1 -levels 3
//   ./hx_compare -r 3 -basis vertex -smoother jacobi -h1 split
//
// Output is CSV. Direct auxiliary solves isolate the decomposition experiment;
// these are not a claim of a scalable auxiliary solver. See reports/hx_compare.md.

#include "common.hpp"
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>

using namespace mfem;
using namespace mfem::jm_hx;
using Clock = std::chrono::steady_clock;

static double Seconds(Clock::time_point start)
{ return std::chrono::duration<double>(Clock::now()-start).count(); }

// MFEM's usual CG tolerance uses the preconditioned norm. Use the same
// unpreconditioned moment-coordinate residual criterion for every configuration.
class CommonResidual : public IterativeSolverController
{
   const Operator &A;
   const Vector &b;
   real_t threshold;
   Vector residual;
public:
   CommonResidual(const Operator &op, const Vector &rhs, real_t relative_tolerance)
      : A(op), b(rhs), threshold(relative_tolerance*rhs.Norml2()), residual(rhs.Size()) { }
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

int main(int argc, char *argv[])
{
   const char *mesh_file = JM_DEFAULT_MESH;
   const char *basis_arg = "all", *smoother_arg = "all", *h1_arg = "all";
   int refinements = 2, levels = 1, repeats = 3, max_iterations = 1000;
   real_t tolerance = 1e-8, damping = 1.0;
   bool random_rhs = true;
   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Straight, connected triangle mesh.");
   args.AddOption(&refinements, "-r", "--refinements", "Initial uniform refinements.");
   args.AddOption(&levels, "-levels", "--levels", "Number of successive refinement levels.");
   args.AddOption(&basis_arg, "-basis", "--basis", "all, moment, or vertex.");
   args.AddOption(&smoother_arg, "-smoother", "--smoother", "all, jacobi, macro, or split.");
   args.AddOption(&h1_arg, "-h1", "--h1", "all, unsplit, or split.");
   args.AddOption(&repeats, "-repeat", "--repeats", "Repeated solves for median timing.");
   args.AddOption(&max_iterations, "-max-it", "--max-iterations", "CG iteration limit.");
   args.AddOption(&tolerance, "-tol", "--tolerance", "Common relative residual tolerance.");
   args.AddOption(&damping, "-damping", "--damping", "Positive smoother weight in HX.");
   args.AddOption(&random_rhs, "-random-rhs", "--random-rhs", "-smooth-rhs",
                  "--smooth-rhs", "Use seeded random moment RHS instead of smooth tensor load.");
   args.ParseCheck(std::cerr);
   const std::string basis_choice(basis_arg), smoother_choice(smoother_arg), h1_choice(h1_arg);
   MFEM_VERIFY(basis_choice == "all" || basis_choice == "moment" || basis_choice == "vertex",
               "invalid basis");
   MFEM_VERIFY(smoother_choice == "all" || smoother_choice == "jacobi" ||
               smoother_choice == "macro" || smoother_choice == "split", "invalid smoother");
   MFEM_VERIFY(h1_choice == "all" || h1_choice == "unsplit" || h1_choice == "split", "invalid H1 mesh");
   MFEM_VERIFY(refinements >= 0 && levels > 0 && repeats > 0 && max_iterations > 0 &&
               tolerance > 0.0 && damping > 0.0, "invalid experiment parameters");
   Mesh mesh(mesh_file);
   // Also validates the mesh restrictions when only unsplit H1 was requested.
   { auto check = AlfeldMesh(mesh); }
   for (int r = 0; r < refinements; r++) { mesh.UniformRefinement(); }
#ifdef MFEM_USE_SUITESPARSE
   std::cerr << "Auxiliary inverses: UMFPACK; timings in seconds; RHS in moment coordinates.\n";
#else
   std::cerr << "Auxiliary inverses: dense LU (limit 2000); enable SuiteSparse for larger runs.\n";
#endif
   std::cout << "refinement,elements,jm_dofs,h1_dofs,basis,smoother,h1,damping,"
             "iterations,residual,converged,moment_assembly_s,vertex_assembly_s,"
             "basis_map_s,aux_setup_s,smoother_setup_s,solve_median_s,apply_median_s\n";
   std::cout << std::setprecision(9);
   bool success = true;
   for (int level = 0; level < levels; level++)
   {
      if (level) { mesh.UniformRefinement(); }
      JohnsonMercierFECollection moment_fec, vertex_fec(JMBasis::SplitVertex);
      FiniteElementSpace moment(&mesh, &moment_fec), vertex(&mesh, &vertex_fec);
      BilinearForm am(&moment), av(&vertex);
      auto start = Clock::now();
      am.AddDomainIntegrator(new MatrixFEMassIntegrator);
      am.AddDomainIntegrator(new MatrixDivDivIntegrator);
      am.Assemble(); am.Finalize();
      const double moment_assembly = Seconds(start);
      start = Clock::now();
      av.AddDomainIntegrator(new MatrixFEMassIntegrator);
      av.AddDomainIntegrator(new MatrixDivDivIntegrator);
      av.Assemble(); av.Finalize();
      const double vertex_assembly = Seconds(start);
      start = Clock::now();
      auto B = BasisMatrix(moment, vertex);
      const double basis_setup = Seconds(start);
      MatrixFunctionCoefficient force(2, [](const Vector &p, DenseMatrix &s)
      {
         s.SetSize(2);
         s(0,0) = std::sin(7*p(0)+3*p(1));
         s(0,1) = s(1,0) = std::cos(5*p(0)-2*p(1));
         s(1,1) = std::sin(4*p(1)-p(0));
      });
      LinearForm load(&moment);
      load.AddDomainIntegrator(new MatrixFEDomainLFIntegrator(force));
      load.Assemble();
      Vector rhs(load);
      if (random_rhs) { rhs.Randomize(1); }
      const real_t rhs_norm = rhs.Norml2();

      for (bool split_h1 : {false, true})
      {
         const std::string h1_name = split_h1 ? "split" : "unsplit";
         if (h1_choice != "all" && h1_choice != h1_name) { continue; }
         start = Clock::now();
         AuxiliaryCorrection auxiliary(moment, vertex, *B, split_h1);
         const double aux_setup = Seconds(start);
         for (bool vertex_basis : {false, true})
         {
            const std::string basis_name = vertex_basis ? "vertex" : "moment";
            if (basis_choice != "all" && basis_choice != basis_name) { continue; }
            for (const std::string kind : {"jacobi", "macro", "split"})
            {
               if (smoother_choice != "all" && smoother_choice != kind) { continue; }
               start = Clock::now();
               std::unique_ptr<Solver> local, mapped;
               const bool use_vertex = vertex_basis || kind == "split";
               auto &fes = use_vertex ? vertex : moment;
               auto &A = use_vertex ? av.SpMat() : am.SpMat();
               if (kind == "jacobi") { local.reset(new DSmoother(A)); }
               else { local.reset(new PatchSolver(A, fes, kind == "split")); }
               if (use_vertex) { mapped.reset(new MappedSolver(*B, *local)); }
               HX hx(mapped ? *mapped : *local, auxiliary, damping);
               const double smoother_setup = Seconds(start);
               Vector solution(moment.GetVSize()), residual(moment.GetVSize());
               std::vector<double> times, apply_times;
               int iterations = 0;
               real_t relative_residual = 0.0;
               bool converged = true;
               for (int repeat = 0; repeat < repeats; repeat++)
               {
                  solution = 0.0;
                  CGSolver cg;
                  CommonResidual monitor(am.SpMat(), rhs, tolerance);
                  cg.SetOperator(am.SpMat()); cg.SetPreconditioner(hx);
                  cg.SetRelTol(0.0); cg.SetAbsTol(0.0); cg.SetController(monitor);
                  cg.SetMaxIter(max_iterations); cg.SetPrintLevel(-1);
                  start = Clock::now(); cg.Mult(rhs, solution);
                  times.push_back(Seconds(start));
                  am.SpMat().Mult(solution, residual); residual -= rhs;
                  relative_residual = residual.Norml2()/rhs_norm;
                  converged = converged && cg.GetConverged() && relative_residual <= tolerance;
                  iterations = cg.GetNumIterations();
                  start = Clock::now(); hx.Mult(rhs, residual);
                  apply_times.push_back(Seconds(start));
               }
               std::sort(times.begin(), times.end());
               std::sort(apply_times.begin(), apply_times.end());
               success = success && converged;
               std::cout << refinements+level << ',' << mesh.GetNE() << ','
                         << moment.GetVSize() << ',' << auxiliary.H1Size() << ','
                         << basis_name << ',' << kind << ',' << h1_name << ',' << damping << ','
                         << iterations << ',' << relative_residual << ',' << converged << ','
                         << moment_assembly << ',' << vertex_assembly << ',' << basis_setup << ','
                         << aux_setup << ',' << smoother_setup << ',' << times[times.size()/2] << ','
                         << apply_times[apply_times.size()/2] << std::endl;
            }
         }
      }
   }
   return success ? 0 : 3;
}

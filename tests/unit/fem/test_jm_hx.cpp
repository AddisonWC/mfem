// Tests for the split-mesh JM HX and multigrid experiments.
#include "unit_tests.hpp"
#include "../../../examples/ex43_hx_compare.hpp"

using namespace mfem;
using namespace mfem::jm_hx;

TEST_CASE("JM HX coordinate maps and split H1 inclusion", "[JohnsonMercier][HX]")
{
   Mesh mesh = Mesh::MakeCartesian2D(2, 2, Element::TRIANGLE, true);
   for (int v = 0; v < mesh.GetNV(); v++)
   {
      real_t *p = mesh.GetVertex(v);
      const real_t x = p[0], y = p[1];
      p[0] = 0.3 + 1.7*x + 0.4*y;
      p[1] = -0.2 + 0.2*x + 0.9*y;
   }
   JohnsonMercierFECollection fm, fv(JMBasis::SplitVertex);
   FiniteElementSpace moments(&mesh, &fm), vertices(&mesh, &fv);
   auto B = BasisMatrix(moments, vertices);
   BilinearForm M(&moments), V(&vertices);
   M.AddDomainIntegrator(new MatrixFEMassIntegrator);
   V.AddDomainIntegrator(new MatrixFEMassIntegrator);
   M.Assemble(); M.Finalize(); V.Assemble(); V.Finalize();
   Vector x(vertices.GetVSize()), bx(moments.GetVSize()), ax(bx.Size()),
          back(x.Size()), vx(x.Size());
   for (int i = 0; i < x.Size(); i++) { x(i) = std::sin(i+1.0); }
   B->Mult(x, bx); M.SpMat().Mult(bx, ax); B->MultTranspose(ax, back);
   V.SpMat().Mult(x, vx); back -= vx;
   REQUIRE(back.Norml2() < 1e-11*vx.Norml2());

   auto patches = Patches(vertices, true);
   Array<int> count(x.Size()); count = 0;
   REQUIRE(patches.size() == size_t(mesh.GetNV()+mesh.GetNE()));
   for (const auto &p : patches) { for (int d : p) { count[d]++; } }
   for (int c : count) { REQUIRE(c == 1); }
   for (int e = 0; e < mesh.GetNE(); e++)
   { REQUIRE(patches[mesh.GetNV()+e].Size() == 3); }

   // Exact block solves are unchanged by this local basis change for macro
   // patches: the moment and vertex DOF lists span the same macro patch space.
   PatchSolver pm(M.SpMat(), moments, false), pv(V.SpMat(), vertices, false);
   MappedSolver mapped(*B, pv);
   Vector ym(bx.Size()), yv(bx.Size());
   pm.Mult(bx, ym); mapped.Mult(bx, yv); yv -= ym;
   REQUIRE(yv.Norml2() < 1e-10*ym.Norml2());

   auto split = AlfeldMesh(mesh);
   REQUIRE(split->GetNE() == 3*mesh.GetNE());
   REQUIRE(split->GetNBE() == mesh.GetNBE());
   H1_FECollection h1_fec(1, 2);
   FiniteElementSpace h1(split.get(), &h1_fec, 3, Ordering::byVDIM);
   auto P = SplitInclusion(h1, vertices);
   DenseMatrix weight(3); weight = 0.0;
   weight(0,0) = weight(2,2) = 1.0; weight(1,1) = 2.0;
   MatrixConstantCoefficient coefficient(weight);
   BilinearForm H(&h1);
   H.AddDomainIntegrator(new VectorMassIntegrator(coefficient));
   H.Assemble(); H.Finalize();
   Vector u(h1.GetVSize()), hu(u.Size()), ptu(u.Size()), pu(x.Size());
   for (int i = 0; i < u.Size(); i++) { u(i) = std::cos(0.7*(i+1)); }
   P->Mult(u, pu); V.SpMat().Mult(pu, vx); P->MultTranspose(vx, ptu);
   H.SpMat().Mult(u, hu); ptu -= hu;
   REQUIRE(ptu.Norml2() < 1e-11*hu.Norml2());

   SECTION("split patch mass energy bounds")
   {
      Vector restricted(x.Size()), action(x.Size());
      for (int sample = 0; sample < 8; sample++)
      {
         for (int i = 0; i < x.Size(); i++)
         { x(i) = std::sin((sample+1)*(i+1.0)); }
         V.SpMat().Mult(x, action);
         const real_t energy = x*action;
         real_t blocks = 0.0;
         for (const auto &p : patches)
         {
            restricted = 0.0;
            for (int d : p) { restricted(d) = x(d); }
            V.SpMat().Mult(restricted, action);
            blocks += restricted*action;
         }
         REQUIRE(energy >= (0.5-1e-12)*blocks);
         REQUIRE(energy <= (2.0+1e-12)*blocks);
      }
   }

   // Check the actual function, including independent values at the centers.
   GridFunction field(&h1); field = u;
   Array<int> dofs;
   const IntegrationRule &rule = IntRules.Get(Geometry::TRIANGLE, 2);
   for (int e = 0; e < mesh.GetNE(); e++)
   {
      vertices.GetElementDofs(e, dofs);
      Vector local; pu.GetSubVector(dofs, local);
      auto &T = *vertices.GetElementTransformation(e);
      InverseElementTransformation inverse(&T);
      for (int c = 0; c < 3; c++)
      {
         auto &S = *h1.GetElementTransformation(3*e+c);
         for (int q = 0; q < rule.GetNPoints(); q++)
         {
            const auto &ip = rule.IntPoint(q);
            Vector point(2), value;
            S.Transform(ip, point);
            IntegrationPoint macro_ip;
            REQUIRE(inverse.Transform(point, macro_ip) == InverseElementTransformation::Inside);
            T.SetIntPoint(&macro_ip);
            DenseTensor shape(2,2,15);
            vertices.GetFE(e)->CalcMShape(T, shape);
            field.GetVectorValue(3*e+c, ip, value);
            for (int comp = 0; comp < 3; comp++)
            {
               const int i = comp == 2 ? 1 : 0, j = comp == 0 ? 0 : 1;
               real_t actual = 0.0;
               for (int k = 0; k < 15; k++) { actual += local(k)*shape(i,j,k); }
               REQUIRE(actual == MFEM_Approx(value(comp)).margin(1e-10));
            }
         }
      }
   }
}

TEST_CASE("JM multigrid patches preserve moment prolongation", "[JohnsonMercier][HX][Multigrid]")
{
   JohnsonMercierFECollection fm, fv(JMBasis::SplitVertex);
   auto mesh = new Mesh(Mesh::MakeCartesian2D(1, 1, Element::TRIANGLE, true));
   auto coarse = new FiniteElementSpace(mesh, &fm);
   FiniteElementSpaceHierarchy hierarchy(mesh, coarse, true, true);
   hierarchy.AddUniformlyRefinedLevel(1, Ordering::byVDIM, Operator::MFEM_SPARSEMAT);
   auto &fine = hierarchy.GetFinestFESpace();
   FiniteElementSpace vertices(fine.GetMesh(), &fv);
   BilinearForm ac(coarse), af(&fine);
   for (auto form : {&ac, &af})
   {
      form->AddDomainIntegrator(new MatrixFEMassIntegrator);
      form->AddDomainIntegrator(new MatrixDivDivIntegrator);
      form->Assemble(); form->Finalize();
   }
   auto B = BasisMatrix(fine, vertices);
   std::unique_ptr<SparseMatrix> av(RAP(*B, af.SpMat(), *B));
   ExactSolver inverse(ac.SpMat());
   PatchSolver macro(af.SpMat(), fine, false), vertex(*av, vertices, false),
               split(*av, vertices, true);
   MappedSolver mapped_macro(*B, vertex), mapped_split(*B, split);
   // Scale symmetrically using a tiny wrapper, as in the example.
   class Scaled : public Solver
   {
      const Solver &R;
   public:
      Scaled(const Solver &r) : Solver(r.Height()), R(r) { }
      void Mult(const Vector &x, Vector &y) const override
      { R.Mult(x, y); y *= 0.33; }
      void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }
      void SetOperator(const Operator &) override { MFEM_ABORT("fixed"); }
   } sm(macro), sv(mapped_macro), ss(mapped_split);
   Array<Operator *> ops, transfers;
   ops.Append(&ac.SpMat()); ops.Append(&af.SpMat());
   transfers.Append(hierarchy.GetProlongationAtLevel(0));
   Array<bool> own_ops(2), own_transfers(1); own_ops = false; own_transfers = false;
   Vector c(coarse->GetVSize()), before(fine.GetVSize()), after(fine.GetVSize());
   for (int i = 0; i < c.Size(); i++) { c(i) = std::cos(i+0.4); }
   transfers[0]->Mult(c, before);
   Vector rhs(fine.GetVSize()), ym(rhs.Size()), yv(rhs.Size()), ys(rhs.Size());
   for (int i = 0; i < rhs.Size(); i++) { rhs(i) = std::sin(i+0.7); }
   Solver *variants[] = {&sm, &sv, &ss};
   Vector *answers[] = {&ym, &yv, &ys};
   for (int k = 0; k < 3; k++)
   {
      Array<Solver *> smoothers; smoothers.Append(&inverse); smoothers.Append(variants[k]);
      Multigrid mg(ops, smoothers, transfers, own_ops, own_ops, own_transfers);
      mg.SetCycleType(Multigrid::CycleType::VCYCLE, 1, 1);
      mg.Mult(rhs, *answers[k]);
      REQUIRE(rhs * *answers[k] > 0.0);
      // Check symmetry of the full V-cycle, including mapped smoothing.
      Vector z(rhs.Size()); mg.Mult(before, z);
      REQUIRE(std::abs(before * *answers[k] - rhs*z) <
              1e-10*before.Norml2()*answers[k]->Norml2());
   }
   yv -= ym;
   REQUIRE(yv.Norml2() < 1e-10*ym.Norml2());
   ys -= ym;
   REQUIRE(ys.Norml2() > 1e-5*ym.Norml2());
   transfers[0]->Mult(c, after); after -= before;
   REQUIRE(after.Norml2() == 0.0);
}

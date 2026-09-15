// Shared implementation for the JM HX comparison example and its tests.
#ifndef MFEM_EX43_HX_COMPARE_HPP
#define MFEM_EX43_HX_COMPARE_HPP

#include "mfem.hpp"
#include <cmath>
#include <memory>
#include <vector>

namespace mfem
{
namespace jm_hx
{

// Split vertices retain the macro vertex numbers; center(T) = NV + T.
inline std::unique_ptr<Mesh> AlfeldMesh(Mesh &macro)
{
   MFEM_VERIFY(macro.Dimension() == 2 && macro.SpaceDimension() == 2 &&
               !macro.Nonconforming() && !macro.GetNodes() && macro.GetNE() > 0,
               "this experiment requires a straight, conforming 2D mesh");
   // One HCT affine gauge is sufficient only on a connected domain.
   const Table &neighbors = macro.ElementToElementTable();
   Array<int> seen(macro.GetNE()); seen = 0;
   std::vector<int> pending(1, 0); seen[0] = 1;
   for (size_t k = 0; k < pending.size(); k++)
   {
      const int e = pending[k];
      for (int j = 0; j < neighbors.RowSize(e); j++)
      {
         const int other = neighbors.GetRow(e)[j];
         if (!seen[other]) { seen[other] = 1; pending.push_back(other); }
      }
   }
   MFEM_VERIFY(pending.size() == size_t(macro.GetNE()),
               "this experiment requires a connected mesh");
   auto split = std::unique_ptr<Mesh>(new Mesh(2, macro.GetNV()+macro.GetNE(),
                                              3*macro.GetNE(), macro.GetNBE()));
   for (int v = 0; v < macro.GetNV(); v++) { split->AddVertex(macro.GetVertex(v)); }
   Array<int> vertices;
   for (int e = 0; e < macro.GetNE(); e++)
   {
      MFEM_VERIFY(macro.GetElementBaseGeometry(e) == Geometry::TRIANGLE,
                  "JM requires triangles");
      macro.GetElementVertices(e, vertices);
      real_t center[2] = {0.0, 0.0};
      for (int v : vertices)
      {
         for (int d = 0; d < 2; d++) { center[d] += macro.GetVertex(v)[d]/3; }
      }
      const int c = split->AddVertex(center);
      for (int j = 0; j < 3; j++)
      {
         split->AddTriangle(vertices[j], vertices[(j+1)%3], c,
                            macro.GetAttribute(e));
      }
   }
   for (int b = 0; b < macro.GetNBE(); b++)
   {
      macro.GetBdrElementVertices(b, vertices);
      split->AddBdrSegment(vertices, macro.GetBdrAttribute(b));
   }
   split->FinalizeTriMesh(1, 0, true);
   return split;
}

/// Map split-vertex coefficients to moment coefficients, on the same mesh.
inline std::unique_ptr<SparseMatrix> BasisMatrix(FiniteElementSpace &moments,
                                                FiniteElementSpace &vertices)
{
   MFEM_VERIFY(moments.GetMesh() == vertices.GetMesh() && moments.GetVDim() == 1 &&
               vertices.GetVDim() == 1, "basis spaces must share a mesh and have vdim=1");
   auto B = std::unique_ptr<SparseMatrix>(new SparseMatrix(moments.GetVSize(),
                                                         vertices.GetVSize()));
   Array<int> rows, cols;
   DenseMatrix C, inverse(15);
   for (int e = 0; e < moments.GetNE(); e++)
   {
      moments.GetElementDofs(e, rows);
      vertices.GetElementDofs(e, cols);
      const auto &fe = dynamic_cast<const JohnsonMercierTriangleFiniteElement &>(
                         *vertices.GetFE(e));
      fe.GetMomentToSplitVertexMatrix(*vertices.GetElementTransformation(e), C);
      DenseMatrixInverse(C).GetInverseMatrix(inverse);
      // Set, rather than add, the shared edge rows. Suppress inversion roundoff
      // in entries that are structurally zero (edge rows cannot use cell data).
      for (int i = 0; i < 15; i++)
      {
         for (int j = 0; j < 15; j++)
         {
            if (i < 12 && (j >= 12 || i/4 != j/4)) { continue; }
            const real_t sign = (rows[i] < 0 ? -1.0 : 1.0)*
                                (cols[j] < 0 ? -1.0 : 1.0);
            B->Set(UnsignIndex(rows[i]), UnsignIndex(cols[j]), sign*inverse(i,j));
         }
      }
   }
   B->Finalize();
   return B;
}

/// Exact inclusion of continuous split P1 tensors into the vertex JM basis.
inline std::unique_ptr<SparseMatrix> SplitInclusion(FiniteElementSpace &h1,
                                                   FiniteElementSpace &jm)
{
   Mesh &macro = *jm.GetMesh();
   MFEM_VERIFY(h1.GetVDim() == 3 && h1.GetOrdering() == Ordering::byVDIM &&
               h1.GetFE(0)->GetOrder() == 1 &&
               h1.GetMesh()->GetNV() == macro.GetNV()+macro.GetNE(),
               "expected an Alfeld P1 symmetric-tensor space");
   MFEM_VERIFY(dynamic_cast<const JohnsonMercierTriangleFiniteElement &>(
                  *jm.GetFE(0)).GetBasisType() == JMBasis::SplitVertex,
               "split inclusion requires the split-vertex JM basis");
   auto P = std::unique_ptr<SparseMatrix>(new SparseMatrix(jm.GetVSize(),
                                                         h1.GetVSize()));
   Array<int> ev, ed, vd;
   for (int e = 0; e < macro.GetNEdges(); e++)
   {
      macro.GetEdgeVertices(e, ev);
      jm.GetEdgeDofs(e, ed);
      const real_t *a = macro.GetVertex(ev[0]), *b = macro.GetVertex(ev[1]);
      const real_t length = std::hypot(b[0]-a[0], b[1]-a[1]);
      const real_t tx = (b[0]-a[0])/length, ty = (b[1]-a[1])/length;
      const real_t nx = ty, ny = -tx;
      const real_t weights[2][3] = {{nx*nx, 2*nx*ny, ny*ny},
         {tx*nx, tx*ny+ty*nx, ty*ny}};
      for (int endpoint = 0; endpoint < 2; endpoint++)
      {
         h1.GetVertexDofs(ev[endpoint], vd);
         for (int row = 0; row < 2; row++)
         {
            for (int c = 0; c < 3; c++)
            {
               P->Set(ed[2*endpoint+row], h1.DofToVDof(vd[0], c), weights[row][c]);
            }
         }
      }
   }
   for (int e = 0; e < macro.GetNE(); e++)
   {
      jm.GetElementInteriorDofs(e, ed);
      h1.GetVertexDofs(macro.GetNV()+e, vd);
      for (int c = 0; c < 3; c++) { P->Set(ed[c], h1.DofToVDof(vd[0], c), 1.0); }
   }
   P->Finalize();
   return P;
}

inline std::vector<Array<int>> Patches(FiniteElementSpace &fes, bool split)
{
   MFEM_VERIFY(!split || dynamic_cast<const JohnsonMercierTriangleFiniteElement &>(
                  *fes.GetFE(0)).GetBasisType() == JMBasis::SplitVertex,
               "split patches require the split-vertex JM basis");
   Mesh &mesh = *fes.GetMesh();
   std::vector<Array<int>> patches(mesh.GetNV() + (split ? mesh.GetNE() : 0));
   Array<int> ev, dofs;
   for (int e = 0; e < mesh.GetNEdges(); e++)
   {
      mesh.GetEdgeVertices(e, ev);
      fes.GetEdgeDofs(e, dofs);
      for (int endpoint = 0; endpoint < 2; endpoint++)
      {
         for (int j = 0; j < (split ? 2 : 4); j++)
         {
            patches[ev[endpoint]].Append(UnsignIndex(dofs[split ? 2*endpoint+j : j]));
         }
      }
   }
   for (int e = 0; e < mesh.GetNE(); e++)
   {
      fes.GetElementInteriorDofs(e, dofs);
      if (split)
      {
         patches[mesh.GetNV()+e] = dofs;
      }
      else
      {
         mesh.GetElementVertices(e, ev);
         for (int v : ev)
         {
            for (int d : dofs) { patches[v].Append(UnsignIndex(d)); }
         }
      }
   }
   return patches;
}

class PatchSolver : public Solver
{
   std::vector<Array<int>> patches;
   std::vector<std::unique_ptr<DenseMatrixInverse>> factors;
   mutable Vector rhs, solution;
public:
   PatchSolver(const SparseMatrix &A, FiniteElementSpace &fes, bool split)
      : Solver(A.Height()), patches(Patches(fes, split))
   {
      for (const auto &p : patches)
      {
         DenseMatrix local(p.Size());
         A.GetSubMatrix(p, p, local);
         // The factorization owns its factors; the source matrix is temporary.
         factors.emplace_back(new DenseMatrixInverse(local, true));
      }
   }
   void Mult(const Vector &x, Vector &y) const override
   {
      y = 0.0;
      for (size_t p = 0; p < patches.size(); p++)
      {
         x.GetSubVector(patches[p], rhs);
         solution.SetSize(rhs.Size());
         factors[p]->Mult(rhs, solution);
         y.AddElementVector(patches[p], solution);
      }
   }
   void SetOperator(const Operator &) override { MFEM_ABORT("fixed patch solver"); }
   void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }
};

/// Express any symmetric inverse in the coordinates of the range of B.
class MappedSolver : public Solver
{
   const SparseMatrix &B;
   const Solver &inverse;
   mutable Vector rhs, solution;
public:
   MappedSolver(const SparseMatrix &map, const Solver &op)
      : Solver(map.Height()), B(map), inverse(op), rhs(map.Width()),
        solution(map.Width()) { }
   void Mult(const Vector &x, Vector &y) const override
   {
      B.MultTranspose(x, rhs);
      inverse.Mult(rhs, solution);
      B.Mult(solution, y);
   }
   void SetOperator(const Operator &) override { MFEM_ABORT("fixed basis map"); }
   void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }
};

class ExactSolver : public Solver
{
#ifdef MFEM_USE_SUITESPARSE
   UMFPackSolver factor;
#else
   DenseMatrix matrix;
   DenseMatrixInverse factor;
#endif
public:
   explicit ExactSolver(SparseMatrix &A) : Solver(A.Height())
#ifdef MFEM_USE_SUITESPARSE
      , factor(A)
#endif
   {
#ifndef MFEM_USE_SUITESPARSE
      MFEM_VERIFY(A.Height() <= 2000,
                  "enable SuiteSparse for auxiliary systems larger than 2000");
      A.ToDenseMatrix(matrix);
      factor.Factor(matrix);
#endif
   }
   void Mult(const Vector &x, Vector &y) const override { factor.Mult(x, y); }
   void SetOperator(const Operator &) override { MFEM_ABORT("fixed direct solve"); }
   void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }
};

/// Auxiliary spaces stay in moment JM coordinates in all experiments.
class AuxiliaryCorrection : public Solver
{
   std::unique_ptr<Mesh> split_mesh;
   H1_FECollection h1_fec;
   FiniteElementSpace h1;
   HCT_FECollection hct_fec;
   FiniteElementSpace hct;
   DenseMatrix weight;
   MatrixConstantCoefficient coefficient;
   BilinearForm h1_form, hct_form;
   std::unique_ptr<SparseMatrix> pi;
   DiscreteLinearOperator airy;
   Array<int> gauge;
   std::unique_ptr<ExactSolver> h1_inverse, hct_inverse;
   mutable Vector rhs1, sol1, rhs2, sol2, correction;

   static DenseMatrix Weight()
   {
      DenseMatrix w(3); w = 0.0;
      w(0,0) = w(2,2) = 1.0; w(1,1) = 2.0;
      return w;
   }
public:
   AuxiliaryCorrection(FiniteElementSpace &moments, FiniteElementSpace &vertices,
                       const SparseMatrix &B, bool split)
      : Solver(moments.GetVSize()),
        split_mesh(split ? AlfeldMesh(*moments.GetMesh()) : nullptr),
        h1_fec(1, 2),
        h1(split ? split_mesh.get() : moments.GetMesh(), &h1_fec, 3, Ordering::byVDIM),
        hct_fec(), hct(moments.GetMesh(), &hct_fec), weight(Weight()),
        coefficient(weight), h1_form(&h1), hct_form(&hct), airy(&hct, &moments),
        rhs1(h1.GetVSize()), sol1(h1.GetVSize()), rhs2(hct.GetVSize()),
        sol2(hct.GetVSize()), correction(height)
   {
      if (split)
      {
         auto inclusion = SplitInclusion(h1, vertices);
         pi.reset(mfem::Mult(B, *inclusion));
      }
      else
      {
         DiscreteLinearOperator interpolation(&h1, &moments);
         interpolation.AddDomainInterpolator(new IdentityInterpolator);
         interpolation.Assemble(); interpolation.Finalize();
         pi.reset(new SparseMatrix(interpolation.SpMat()));
      }
      airy.AddDomainInterpolator(new AiryInterpolator);
      airy.Assemble(); airy.Finalize();
      h1_form.AddDomainIntegrator(new VectorMassIntegrator(coefficient));
      h1_form.AddDomainIntegrator(new VectorDiffusionIntegrator(coefficient));
      h1_form.Assemble(); h1_form.Finalize();
      h1_inverse.reset(new ExactSolver(h1_form.SpMat()));
      hct_form.AddDomainIntegrator(new HessianIntegrator);
      hct_form.Assemble(); hct_form.Finalize();
      hct.GetVertexDofs(0, gauge);
      for (int d : gauge) { hct_form.SpMat().EliminateRowCol(UnsignIndex(d)); }
      hct_inverse.reset(new ExactSolver(hct_form.SpMat()));
   }
   int H1Size() const { return h1.GetVSize(); }
   void Mult(const Vector &x, Vector &y) const override
   {
      pi->MultTranspose(x, rhs1); h1_inverse->Mult(rhs1, sol1); pi->Mult(sol1, y);
      airy.MultTranspose(x, rhs2);
      for (int d : gauge) { rhs2(UnsignIndex(d)) = 0.0; }
      hct_inverse->Mult(rhs2, sol2); airy.Mult(sol2, correction); y += correction;
   }
   void SetOperator(const Operator &) override { MFEM_ABORT("fixed auxiliary spaces"); }
   void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }
};

class HX : public Solver
{
   const Solver &smoother, &auxiliary;
   real_t damping;
   mutable Vector correction;
public:
   HX(const Solver &R, const Solver &aux, real_t scale)
      : Solver(R.Height()), smoother(R), auxiliary(aux), damping(scale),
        correction(height) { }
   void Mult(const Vector &x, Vector &y) const override
   {
      smoother.Mult(x, y); y *= damping;
      auxiliary.Mult(x, correction); y += correction;
   }
   void SetOperator(const Operator &) override { MFEM_ABORT("fixed HX operator"); }
   void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }
};

} // namespace jm_hx
} // namespace mfem
#endif

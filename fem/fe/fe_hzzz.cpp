// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.


#include "fe_hzzz.hpp"
#include "../eltrans.hpp"
#include <cmath>

namespace mfem
{
namespace
{
constexpr real_t vertices[3][2] = {{0,0}, {1,0}, {0,1}};
constexpr int edges[3][2] = {{0,1}, {1,2}, {2,0}};
constexpr int retained[21] = {0,1,2,3,4,5,6,7,8,9,10,11,13,14,15,17,18,19,21,22,23};

// Embed the reduced basis in HZZZ: R(:,i) adds to HZZZ basis function i the
// interior bubbles 18-20 that make its divergence L2-orthogonal to theirs.
// div_shape(ip,div) gives the 21x2 divergences at ip; they are linear, so
// an order-2 rule is exact and the constant map weight cancels.
template <typename DivShape>
void BuildReduction(DivShape div_shape, DenseMatrix &R)
{
   const IntegrationRule &ir = IntRules.Get(Geometry::TRIANGLE,2);
   DenseMatrix div(21,2), inner(21,3);
   inner = 0.0;
   for (int q = 0; q < ir.GetNPoints(); q++)
   {
      const IntegrationPoint &ip = ir.IntPoint(q);
      div_shape(ip,div);
      for (int i = 0; i < 21; i++)
      {
         for (int j = 0; j < 3; j++)
         {
            inner(i,j) += ip.weight*(div(i,0)*div(18+j,0)+div(i,1)*div(18+j,1));
         }
      }
   }
   DenseMatrix gram(3);
   for (int j = 0; j < 3; j++)
   {
      for (int k = 0; k < 3; k++) { gram(j,k) = inner(18+j,k); }
   }
   DenseMatrixInverse gram_inverse(gram);
   Vector rhs(3), alpha(3);
   R.SetSize(21,18);
   R = 0.0;
   for (int i = 0; i < 18; i++)
   {
      R(i,i) = 1.0;
      for (int j = 0; j < 3; j++) { rhs(j) = -inner(i,j); }
      gram_inverse.Mult(rhs,alpha);
      for (int j = 0; j < 3; j++) { R(18+j,i) = alpha(j); }
   }
}

void ReduceMShape(const DenseMatrix &R, const DenseTensor &full,
                  DenseTensor &shape)
{
   shape = 0.0;
   for (int k = 0; k < R.Width(); k++)
   {
      for (int l = 0; l < R.Height(); l++)
      {
         if (R(l,k) == 0.0) { continue; }
         for (int i = 0; i < 2; i++)
         {
            for (int j = 0; j < 2; j++) { shape(i,j,k) += R(l,k)*full(i,j,l); }
         }
      }
   }
}
}

HuangZhangZhouZhuTriangleFiniteElement::HuangZhangZhouZhuTriangleFiniteElement()
   : FiniteElement(2, Geometry::TRIANGLE, 21, 3, FunctionSpace::Pk)
{
   range_type = MATRIX;
   map_type = DOUBLE_CONTRAVARIANT_PIOLA;
   deriv_type = DIV;
   deriv_range_type = VECTOR;
   deriv_map_type = UNKNOWN_MAP_TYPE;
   vdim = 2;
   for (int i = 0; i < dof; i++)
   { Nodes.IntPoint(i) = aw.GetNodes().IntPoint(retained[i]); }
   DenseMatrix J(2);
   J = 0.0; J(0,0) = J(1,1) = 1.0;
   GetEmbedding(J,reference_embedding);
}

void HuangZhangZhouZhuTriangleFiniteElement::GetEmbedding(
   const DenseMatrix &J, DenseMatrix &E) const
{
   E.SetSize(24,21);
   E = 0.0;
   for (int i = 0; i < dof; i++) { E(retained[i],i) = 1.0; }
   for (int e = 0; e < 3; e++)
   {
      const int a = edges[e][0], b = edges[e][1];
      const real_t dx = vertices[b][0]-vertices[a][0];
      const real_t dy = vertices[b][1]-vertices[a][1];
      const real_t tx = J(0,0)*dx+J(0,1)*dy;
      const real_t ty = J(1,0)*dx+J(1,1)*dy;
      const real_t scale = 1.0/(6*std::hypot(tx,ty));
      // For quadratic q(s), integral_0^1 (2s-1)q(s) ds = (q(1)-q(0))/6.
      // This removes the cubic nt trace from the physical AW space.
      const real_t nt[3] = {tx*ty, ty*ty-tx*tx, -tx*ty};
      for (int c = 0; c < 3; c++)
      {
         E(12+4*e,3*a+c) = -scale*nt[c];
         E(12+4*e,3*b+c) = scale*nt[c];
      }
   }
}

void HuangZhangZhouZhuTriangleFiniteElement::ReduceShape(
   const DenseMatrix &E, const DenseTensor &full, DenseTensor &shape) const
{
   shape = 0.0;
   for (int k = 0; k < dof; k++)
   {
      for (int l = 0; l < 24; l++)
      {
         for (int i = 0; i < 2; i++)
         {
            for (int j = 0; j < 2; j++) { shape(i,j,k) += E(l,k)*full(i,j,l); }
         }
      }
   }
}

void HuangZhangZhouZhuTriangleFiniteElement::SelectDofs(
   const DenseMatrix &full, DenseMatrix &I) const
{
   I.SetSize(dof,full.Width());
   for (int i = 0; i < dof; i++)
   {
      for (int j = 0; j < full.Width(); j++) { I(i,j) = full(retained[i],j); }
   }
}

void HuangZhangZhouZhuTriangleFiniteElement::CalcMShape(
   const IntegrationPoint &ip, DenseTensor &shape) const
{
   DenseTensor full(2,2,24);
   aw.CalcMShape(ip,full);
   ReduceShape(reference_embedding,full,shape);
}

void HuangZhangZhouZhuTriangleFiniteElement::CalcMShape(
   ElementTransformation &T, DenseTensor &shape) const
{
   DenseTensor full(2,2,24);
   DenseMatrix E;
   aw.CalcMShape(T,full);
   GetEmbedding(T.Jacobian(),E);
   ReduceShape(E,full,shape);
}

void HuangZhangZhouZhuTriangleFiniteElement::CalcDivShape(
   const IntegrationPoint &ip, DenseMatrix &shape) const
{
   DenseMatrix full(24,2);
   aw.CalcDivShape(ip,full);
   MultAtB(reference_embedding,full,shape);
}

void HuangZhangZhouZhuTriangleFiniteElement::CalcPhysDivShape(
   ElementTransformation &T, DenseMatrix &shape) const
{
   DenseMatrix full(24,2), E;
   aw.CalcPhysDivShape(T,full);
   GetEmbedding(T.Jacobian(),E);
   MultAtB(E,full,shape);
}

void HuangZhangZhouZhuTriangleFiniteElement::Project(
   const FiniteElement &fe, ElementTransformation &T, DenseMatrix &I) const
{
   DenseMatrix full;
   aw.Project(fe,T,full);
   SelectDofs(full,I);
}

void HuangZhangZhouZhuTriangleFiniteElement::GetTransferMatrix(
   const FiniteElement &fe, ElementTransformation &T, DenseMatrix &I) const
{
   MFEM_VERIFY(dynamic_cast<const HuangZhangZhouZhuTriangleFiniteElement *>(&fe),
               "HZZZ transfer requires an HZZZ source element");
   DenseMatrix full, reduced(24,dof);
   aw.GetLocalInterpolation(T,full);
   Mult(full,reference_embedding,reduced);
   SelectDofs(reduced,I);
}

void HuangZhangZhouZhuTriangleFiniteElement::GetPhysicalTransferMatrix(
   const DenseMatrix &, ElementTransformation &child,
   ElementTransformation &fine, DenseMatrix &I) const
{
   const auto &center = Geometries.GetCenter(Geometry::TRIANGLE);
   child.SetIntPoint(&center);
   fine.SetIntPoint(&center);
   MFEM_VERIFY(child.GetSpaceDim() == 2 && fine.GetSpaceDim() == 2 &&
               child.Hessian().FNorm2() < 1e-20 && fine.Hessian().FNorm2() < 1e-20,
               "HZZZ transfer requires affine 2D transformations");
   DenseMatrix K(2), J(2), E, reference, physical, reduced(24,dof);
   CalcInverse(child.Jacobian(),K);
   Mult(fine.Jacobian(),K,J);
   GetEmbedding(J,E);
   // Retain the full AW transfer until after the physical change of basis:
   // the eliminated nt moments depend on the coarse physical geometry.
   aw.GetLocalInterpolation(child,reference);
   aw.GetPhysicalTransferMatrix(reference,child,fine,physical);
   Mult(physical,E,reduced);
   SelectDofs(reduced,I);
}

void HuangZhangZhouZhuTriangleFiniteElement::GetFaceDofs(
   int face, int **dofs, int *ndofs) const
{
   static int indices[3][9] = {{0,1,2,3,4,5,9,10,11},
      {3,4,5,6,7,8,12,13,14},
      {6,7,8,0,1,2,15,16,17}
   };
   MFEM_ASSERT(face >= 0 && face < 3, "invalid face index");
   *dofs = indices[face];
   *ndofs = 9;
}

ReducedHuangZhangZhouZhuTriangleFiniteElement::
ReducedHuangZhangZhouZhuTriangleFiniteElement()
   : FiniteElement(2, Geometry::TRIANGLE, 18, 3, FunctionSpace::Pk)
{
   range_type = MATRIX;
   map_type = DOUBLE_CONTRAVARIANT_PIOLA;
   deriv_type = DIV;
   deriv_range_type = VECTOR;
   deriv_map_type = UNKNOWN_MAP_TYPE;
   vdim = 2;
   for (int i = 0; i < dof; i++) { Nodes.IntPoint(i) = hzzz.GetNodes().IntPoint(i); }
   BuildReduction([&](const IntegrationPoint &ip, DenseMatrix &div)
   { hzzz.CalcDivShape(ip,div); }, reference_reduction);
}

void ReducedHuangZhangZhouZhuTriangleFiniteElement::GetReduction(
   ElementTransformation &T, DenseMatrix &R) const
{
   // Callers keep a pointer to their integration point in T; restore it.
   const IntegrationPoint &caller_ip = T.GetIntPoint();
   BuildReduction([&](const IntegrationPoint &ip, DenseMatrix &div)
   {
      T.SetIntPoint(&ip);
      hzzz.CalcPhysDivShape(T,div);
   }, R);
   T.SetIntPoint(&caller_ip);
}

void ReducedHuangZhangZhouZhuTriangleFiniteElement::CalcMShape(
   const IntegrationPoint &ip, DenseTensor &shape) const
{
   DenseTensor full(2,2,21);
   hzzz.CalcMShape(ip,full);
   ReduceMShape(reference_reduction,full,shape);
}

void ReducedHuangZhangZhouZhuTriangleFiniteElement::CalcMShape(
   ElementTransformation &T, DenseTensor &shape) const
{
   DenseTensor full(2,2,21);
   DenseMatrix R;
   hzzz.CalcMShape(T,full);
   GetReduction(T,R);
   ReduceMShape(R,full,shape);
}

void ReducedHuangZhangZhouZhuTriangleFiniteElement::CalcDivShape(
   const IntegrationPoint &ip, DenseMatrix &shape) const
{
   DenseMatrix full(21,2);
   hzzz.CalcDivShape(ip,full);
   MultAtB(reference_reduction,full,shape);
}

void ReducedHuangZhangZhouZhuTriangleFiniteElement::CalcPhysDivShape(
   ElementTransformation &T, DenseMatrix &shape) const
{
   DenseMatrix full(21,2), R;
   hzzz.CalcPhysDivShape(T,full);
   GetReduction(T,R);
   MultAtB(R,full,shape);
}

void ReducedHuangZhangZhouZhuTriangleFiniteElement::Project(
   const FiniteElement &fe, ElementTransformation &T, DenseMatrix &I) const
{
   // The reduced basis is dual to the HZZZ vertex and edge DOFs, since the
   // added interior bubbles have vanishing boundary DOFs.
   DenseMatrix full;
   hzzz.Project(fe,T,full);
   I.SetSize(dof,full.Width());
   for (int i = 0; i < dof; i++)
   {
      for (int j = 0; j < full.Width(); j++) { I(i,j) = full(i,j); }
   }
}

} // namespace mfem

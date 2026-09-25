// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

// Brezzi-Douglas-Marini finite element classes

#include "fe_bdm.hpp"

namespace mfem
{

const real_t BDM_TriangleElement::nk[14] =
{ 0.,-1.,  1.,1.,  -1.,0.,  1.,0.,  -1.,1.,  0.,-1.,  0.,1. };

BDM_TriangleElement::BDM_TriangleElement(const int p, const int ob_type)
   : VectorFiniteElement(2, Geometry::TRIANGLE, (p + 1)*(p + 2), p,
                         H_DIV, FunctionSpace::Pk),
     dof2nk(dof)
{
   MFEM_VERIFY(p >= 1, "BDM_TriangleElement requires order >= 1.");
   const real_t *bop = poly1d.OpenPoints(p, VerifyOpen(VerifyNodal(ob_type)));
   const int q = p - 1;
   const real_t *eop = (q > 0) ? poly1d.OpenPoints(q - 1) : nullptr;
   const real_t *iop = (q > 1) ? poly1d.OpenPoints(q - 2) : nullptr;

#ifndef MFEM_THREAD_SAFE
   shape_x.SetSize(p + 1); shape_y.SetSize(p + 1); shape_l.SetSize(p + 1);
   dshape_x.SetSize(p + 1); dshape_y.SetSize(p + 1);
   dshape_l.SetSize(p + 1);
   u.SetSize(dof, dim); divu.SetSize(dof);
#else
   Vector shape_x(p + 1), shape_y(p + 1), shape_l(p + 1);
#endif

   int o = 0;
   for (int i = 0; i <= p; i++)
   { Nodes.IntPoint(o).Set2(bop[i], 0.); dof2nk[o++] = 0; }
   for (int i = 0; i <= p; i++)
   { Nodes.IntPoint(o).Set2(bop[p-i], bop[i]); dof2nk[o++] = 1; }
   for (int i = 0; i <= p; i++)
   { Nodes.IntPoint(o).Set2(0., bop[p-i]); dof2nk[o++] = 2; }

   // Cell-local tangential samples, using the ND_{p-1} point layout.
   for (int i = 0; i < q; i++)
   { Nodes.IntPoint(o).Set2(eop[i], 0.); dof2nk[o++] = 3; }
   for (int i = 0; i < q; i++)
   { Nodes.IntPoint(o).Set2(eop[q-1-i], eop[i]); dof2nk[o++] = 4; }
   for (int i = 0; i < q; i++)
   { Nodes.IntPoint(o).Set2(0., eop[q-1-i]); dof2nk[o++] = 5; }
   for (int j = 0; j < q - 1; j++)
      for (int i = 0; i + j < q - 1; i++)
      {
         const real_t w = iop[i] + iop[j] + iop[q-2-i-j];
         Nodes.IntPoint(o).Set2(iop[i]/w, iop[j]/w); dof2nk[o++] = 3;
         Nodes.IntPoint(o).Set2(iop[i]/w, iop[j]/w); dof2nk[o++] = 6;
      }
   MFEM_ASSERT(o == dof, "incorrect number of BDM triangle dofs");

   DenseMatrix T(dof);
   for (int m = 0; m < dof; m++)
   {
      const IntegrationPoint &ip = Nodes.IntPoint(m);
      poly1d.CalcBasis(p, ip.x, shape_x);
      poly1d.CalcBasis(p, ip.y, shape_y);
      poly1d.CalcBasis(p, 1. - ip.x - ip.y, shape_l);
      const real_t *nm = nk + 2*dof2nk[m];
      o = 0;
      for (int j = 0; j <= p; j++)
         for (int i = 0; i + j <= p; i++)
         {
            const real_t s = shape_x(i)*shape_y(j)*shape_l(p-i-j);
            T(o++,m) = s*nm[0]; T(o++,m) = s*nm[1];
         }
   }
   Ti.Factor(T);
}

void BDM_TriangleElement::CalcVShape(const IntegrationPoint &ip,
                                   DenseMatrix &shape) const
{
   const int p = order;
#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p + 1), shape_y(p + 1), shape_l(p + 1);
   DenseMatrix u(dof, dim);
#endif
   poly1d.CalcBasis(p, ip.x, shape_x);
   poly1d.CalcBasis(p, ip.y, shape_y);
   poly1d.CalcBasis(p, 1. - ip.x - ip.y, shape_l);
   int o = 0;
   for (int j = 0; j <= p; j++)
      for (int i = 0; i + j <= p; i++)
      {
         const real_t s = shape_x(i)*shape_y(j)*shape_l(p-i-j);
         u(o,0) = s; u(o,1) = 0.; o++;
         u(o,0) = 0.; u(o,1) = s; o++;
      }
   Ti.Mult(u, shape);
}

void BDM_TriangleElement::CalcDivShape(const IntegrationPoint &ip,
                                     Vector &divshape) const
{
   const int p = order;
#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p + 1), shape_y(p + 1), shape_l(p + 1);
   Vector dshape_x(p + 1), dshape_y(p + 1), dshape_l(p + 1);
   Vector divu(dof);
#endif
   poly1d.CalcBasis(p, ip.x, shape_x, dshape_x);
   poly1d.CalcBasis(p, ip.y, shape_y, dshape_y);
   poly1d.CalcBasis(p, 1. - ip.x - ip.y, shape_l, dshape_l);
   int o = 0;
   for (int j = 0; j <= p; j++)
      for (int i = 0; i + j <= p; i++)
      {
         const int l = p-i-j;
         divu(o++) = (dshape_x(i)*shape_l(l) -
                      shape_x(i)*dshape_l(l))*shape_y(j);
         divu(o++) = (dshape_y(j)*shape_l(l) -
                      shape_y(j)*dshape_l(l))*shape_x(i);
      }
   Ti.Mult(divu, divshape);
}


const real_t BDM_TetrahedronElement::nk[30] =
{ 1.,1.,1.,  -1.,0.,0.,  0.,-1.,0.,  0.,0.,-1.,
  1.,0.,0.,  0.,1.,0.,  0.,0.,1.,  -1.,1.,0.,  -1.,0.,1.,  0.,-1.,1. };

BDM_TetrahedronElement::BDM_TetrahedronElement(const int p,
                                             const int ob_type)
   : VectorFiniteElement(3, Geometry::TETRAHEDRON,
                         (p + 1)*(p + 2)*(p + 3)/2, p,
                         H_DIV, FunctionSpace::Pk),
     dof2nk(dof)
{
   MFEM_VERIFY(p >= 1, "BDM_TetrahedronElement requires order >= 1.");
   const real_t *bop = poly1d.OpenPoints(p, VerifyOpen(VerifyNodal(ob_type)));
   const int q = p - 1;
   const real_t *eop = (q > 0) ? poly1d.OpenPoints(q - 1) : nullptr;
   const real_t *fop = (q > 1) ? poly1d.OpenPoints(q - 2) : nullptr;
   const real_t *iop = (q > 2) ? poly1d.OpenPoints(q - 3) : nullptr;

#ifndef MFEM_THREAD_SAFE
   shape_x.SetSize(p + 1); shape_y.SetSize(p + 1);
   shape_z.SetSize(p + 1); shape_l.SetSize(p + 1);
   dshape_x.SetSize(p + 1); dshape_y.SetSize(p + 1);
   dshape_z.SetSize(p + 1); dshape_l.SetSize(p + 1);
   u.SetSize(dof, dim); divu.SetSize(dof);
#else
   Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_l(p + 1);
#endif

   int o = 0;
   for (int j = 0; j <= p; j++) for (int i = 0; i + j <= p; i++)
   {
      const real_t w = bop[i] + bop[j] + bop[p-i-j];
      Nodes.IntPoint(o).Set3(bop[p-i-j]/w, bop[i]/w, bop[j]/w);
      dof2nk[o++] = 0;
   }
   for (int j = 0; j <= p; j++) for (int i = 0; i + j <= p; i++)
   {
      const real_t w = bop[i] + bop[j] + bop[p-i-j];
      Nodes.IntPoint(o).Set3(0., bop[j]/w, bop[i]/w); dof2nk[o++] = 1;
   }
   for (int j = 0; j <= p; j++) for (int i = 0; i + j <= p; i++)
   {
      const real_t w = bop[i] + bop[j] + bop[p-i-j];
      Nodes.IntPoint(o).Set3(bop[i]/w, 0., bop[j]/w); dof2nk[o++] = 2;
   }
   for (int j = 0; j <= p; j++) for (int i = 0; i + j <= p; i++)
   {
      const real_t w = bop[i] + bop[j] + bop[p-i-j];
      Nodes.IntPoint(o).Set3(bop[j]/w, bop[i]/w, 0.); dof2nk[o++] = 3;
   }

   // Cell-local edge and face tangential samples, as in ND_{p-1}.
   const int edge_dir[6] = {4,5,6,7,8,9};
   for (int e = 0; e < 6; e++) for (int i = 0; i < q; i++)
   {
      const real_t a = eop[i], b = eop[q-1-i];
      switch (e)
      {
         case 0: Nodes.IntPoint(o).Set3(a,0.,0.); break;
         case 1: Nodes.IntPoint(o).Set3(0.,a,0.); break;
         case 2: Nodes.IntPoint(o).Set3(0.,0.,a); break;
         case 3: Nodes.IntPoint(o).Set3(b,a,0.); break;
         case 4: Nodes.IntPoint(o).Set3(b,0.,a); break;
         default: Nodes.IntPoint(o).Set3(0.,b,a); break;
      }
      dof2nk[o++] = edge_dir[e];
   }

   const int face_dir[4][2] = {{7,8},{6,5},{4,6},{5,4}};
   for (int f = 0; f < 4; f++)
      for (int j = 0; j < q - 1; j++)
         for (int i = 0; i + j < q - 1; i++)
         {
            const real_t w = fop[i] + fop[j] + fop[q-2-i-j];
            IntegrationPoint pt;
            switch (f)
            {
               case 0: pt.Set3(fop[q-2-i-j]/w,fop[i]/w,fop[j]/w); break;
               case 1: pt.Set3(0.,fop[j]/w,fop[i]/w); break;
               case 2: pt.Set3(fop[i]/w,0.,fop[j]/w); break;
               default: pt.Set3(fop[j]/w,fop[i]/w,0.); break;
            }
            Nodes.IntPoint(o) = pt; dof2nk[o++] = face_dir[f][0];
            Nodes.IntPoint(o) = pt; dof2nk[o++] = face_dir[f][1];
         }
   for (int k = 0; k < q - 2; k++)
      for (int j = 0; j + k < q - 2; j++)
         for (int i = 0; i + j + k < q - 2; i++)
         {
            const real_t w = iop[i] + iop[j] + iop[k] + iop[q-3-i-j-k];
            IntegrationPoint pt;
            pt.Set3(iop[i]/w, iop[j]/w, iop[k]/w);
            for (int d = 0; d < 3; d++)
            { Nodes.IntPoint(o) = pt; dof2nk[o++] = 4+d; }
         }
   MFEM_ASSERT(o == dof, "incorrect number of BDM tetrahedron dofs");

   DenseMatrix T(dof);
   for (int m = 0; m < dof; m++)
   {
      const IntegrationPoint &ip = Nodes.IntPoint(m);
      poly1d.CalcBasis(p, ip.x, shape_x);
      poly1d.CalcBasis(p, ip.y, shape_y);
      poly1d.CalcBasis(p, ip.z, shape_z);
      poly1d.CalcBasis(p, 1. - ip.x - ip.y - ip.z, shape_l);
      const real_t *nm = nk + 3*dof2nk[m];
      o = 0;
      for (int k = 0; k <= p; k++)
         for (int j = 0; j + k <= p; j++)
            for (int i = 0; i + j + k <= p; i++)
            {
               const real_t s = shape_x(i)*shape_y(j)*shape_z(k)*
                                shape_l(p-i-j-k);
               T(o++,m) = s*nm[0];
               T(o++,m) = s*nm[1];
               T(o++,m) = s*nm[2];
            }
   }
   Ti.Factor(T);
}

void BDM_TetrahedronElement::CalcVShape(const IntegrationPoint &ip,
                                      DenseMatrix &shape) const
{
   const int p = order;
#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_l(p + 1);
   DenseMatrix u(dof, dim);
#endif
   poly1d.CalcBasis(p, ip.x, shape_x);
   poly1d.CalcBasis(p, ip.y, shape_y);
   poly1d.CalcBasis(p, ip.z, shape_z);
   poly1d.CalcBasis(p, 1. - ip.x - ip.y - ip.z, shape_l);
   int o = 0;
   for (int k = 0; k <= p; k++)
      for (int j = 0; j + k <= p; j++)
         for (int i = 0; i + j + k <= p; i++)
         {
            const real_t s = shape_x(i)*shape_y(j)*shape_z(k)*shape_l(p-i-j-k);
            u(o,0) = s;  u(o,1) = 0.; u(o,2) = 0.; o++;
            u(o,0) = 0.; u(o,1) = s;  u(o,2) = 0.; o++;
            u(o,0) = 0.; u(o,1) = 0.; u(o,2) = s;  o++;
         }
   Ti.Mult(u, shape);
}

void BDM_TetrahedronElement::CalcDivShape(const IntegrationPoint &ip,
                                        Vector &divshape) const
{
   const int p = order;
#ifdef MFEM_THREAD_SAFE
   Vector shape_x(p + 1), shape_y(p + 1), shape_z(p + 1), shape_l(p + 1);
   Vector dshape_x(p + 1), dshape_y(p + 1), dshape_z(p + 1), dshape_l(p + 1);
   Vector divu(dof);
#endif
   poly1d.CalcBasis(p, ip.x, shape_x, dshape_x);
   poly1d.CalcBasis(p, ip.y, shape_y, dshape_y);
   poly1d.CalcBasis(p, ip.z, shape_z, dshape_z);
   poly1d.CalcBasis(p, 1. - ip.x - ip.y - ip.z, shape_l, dshape_l);
   int o = 0;
   for (int k = 0; k <= p; k++)
      for (int j = 0; j + k <= p; j++)
         for (int i = 0; i + j + k <= p; i++)
         {
            const int l = p-i-j-k;
            divu(o++) = (dshape_x(i)*shape_l(l) - shape_x(i)*dshape_l(l))*
                        shape_y(j)*shape_z(k);
            divu(o++) = (dshape_y(j)*shape_l(l) - shape_y(j)*dshape_l(l))*
                        shape_x(i)*shape_z(k);
            divu(o++) = (dshape_z(k)*shape_l(l) - shape_z(k)*dshape_l(l))*
                        shape_x(i)*shape_y(j);
         }
   Ti.Mult(divu, divshape);
}

} // namespace mfem

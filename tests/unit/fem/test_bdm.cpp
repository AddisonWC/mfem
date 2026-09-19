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

#include "mfem.hpp"
#include "unit_tests.hpp"

using namespace mfem;

namespace
{

Mesh MakeBDMMesh(int dim)
{
   Mesh mesh = dim == 2 ?
               Mesh::MakeCartesian2D(2, 2, Element::TRIANGLE, true) :
               Mesh::MakeCartesian3D(1, 1, 1, Element::TETRAHEDRON);
   // Exercise Piola mappings with shear and nonuniform scaling.
   for (int i = 0; i < mesh.GetNV(); i++)
   {
      real_t *x = mesh.GetVertex(i);
      x[0] = 1.3*x[0] + 0.2*x[1];
      x[1] *= 0.7;
      if (dim == 3)
      {
         x[0] += 0.1*x[2];
         x[2] *= 1.1;
      }
   }
   return mesh;
}

} // namespace

TEST_CASE("BDM polynomial reproduction", "[BDM]")
{
   const int dim = GENERATE(2, 3);
   const int p = GENERATE(1, 2, 3, 4, 5);
   const int basis = GENERATE(BasisType::GaussLegendre, BasisType::OpenUniform,
                              BasisType::OpenHalfUniform);
   CAPTURE(dim, p, basis);
   Mesh mesh = MakeBDMMesh(dim);
   BDM_FECollection fec(p, dim, basis);
   FiniteElementSpace fes(&mesh, &fec);
   GridFunction field(&fes);
   const IntegrationRule &ir = IntRules.Get(mesh.GetElementGeometry(0), 3);

   // Span all of [P_p]^dim, checking every component and its divergence.
   for (int c = 0; c < dim; c++)
      for (int a = 0; a <= p; a++)
         for (int b = 0; a + b <= p; b++)
            for (int k = 0; k <= (dim == 3 ? p-a-b : 0); k++)
            {
               CAPTURE(c, a, b, k);
               VectorFunctionCoefficient coefficient(dim,
               [=](const Vector &x, Vector &v)
               {
                  v = 0.0;
                  v[c] = std::pow(x[0], a)*std::pow(x[1], b);
                  if (dim == 3) { v[c] *= std::pow(x[2], k); }
               });
               field.ProjectCoefficient(coefficient);
               real_t value_error = 0.0, divergence_error = 0.0;
               Vector value(dim), exact(dim), x(dim);
               for (int e = 0; e < mesh.GetNE(); e++)
               {
                  ElementTransformation *T = mesh.GetElementTransformation(e);
                  for (int j = 0; j < ir.GetNPoints(); j++)
                  {
                     const IntegrationPoint &ip = ir.IntPoint(j);
                     T->SetIntPoint(&ip);
                     field.GetVectorValue(*T, ip, value);
                     coefficient.Eval(exact, *T, ip);
                     value -= exact;
                     value_error = std::max(value_error, value.Normlinf());

                     T->Transform(ip, x);
                     const int powers[3] = {a, b, k};
                     real_t divergence = powers[c];
                     if (powers[c] > 0)
                     {
                        for (int d = 0; d < dim; d++)
                        {
                           divergence *= std::pow(x[d], powers[d] - (d == c));
                        }
                     }
                     const real_t div_error =
                        std::abs(field.GetDivergence(*T) - divergence);
                     divergence_error = std::max(divergence_error, div_error);
                  }
               }
               REQUIRE(value_error < 2e-10);
               REQUIRE(divergence_error < 2e-9);
            }
}

TEST_CASE("BDM normal trace compatibility", "[BDM]")
{
   const int dim = GENERATE(2, 3);
   const int p = GENERATE(1, 2, 3, 4, 5);
   const int basis = GENERATE(BasisType::GaussLegendre, BasisType::OpenUniform,
                              BasisType::OpenHalfUniform);
   CAPTURE(dim, p, basis);
   Mesh mesh = MakeBDMMesh(dim);
   BDM_FECollection fec(p, dim, basis);
   FiniteElementSpace fes(&mesh, &fec);
   GridFunction field(&fes);
   for (int i = 0; i < field.Size(); i++) { field[i] = std::sin(real_t(i+1)); }
   std::unique_ptr<FiniteElementCollection> trace(fec.GetTraceCollection());

   // Arbitrary cell-local coefficients must not affect the normal trace.
   // Compare both volume restrictions to the scalar trace on every face.
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      FaceElementTransformations *T = mesh.GetFaceElementTransformations(f);
      const FiniteElement *face_fe =
         trace->FiniteElementForGeometry(T->GetGeometryType());
      Array<int> dofs;
      fes.GetFaceDofs(f, dofs);
      Vector face_values, shape(face_fe->GetDof()), normal(dim), value(dim);
      field.GetSubVector(dofs, face_values);
      const IntegrationRule &ir = IntRules.Get(T->GetGeometryType(), 2*p+2);
      real_t error = 0.0;
      for (int j = 0; j < ir.GetNPoints(); j++)
      {
         const IntegrationPoint &ip = ir.IntPoint(j);
         T->SetAllIntPoints(&ip);
         CalcOrtho(T->Jacobian(), normal);
         face_fe->CalcShape(ip, shape);
         const real_t normal_value = shape*face_values;
         field.GetVectorValue(*T->Elem1, T->GetElement1IntPoint(), value);
         error = std::max(error, std::abs(value*normal - normal_value));
         if (T->Elem2No >= 0)
         {
            field.GetVectorValue(*T->Elem2, T->GetElement2IntPoint(), value);
            error = std::max(error, std::abs(value*normal - normal_value));
         }
      }
      CAPTURE(f);
      REQUIRE(error < 2e-9);
   }
}

TEST_CASE("BDM local transfers", "[BDM]")
{
   const int dim = GENERATE(2, 3);
   const int p = GENERATE(1, 2, 3, 4, 5);
   CAPTURE(dim, p);
   const Geometry::Type geom = dim == 2 ? Geometry::TRIANGLE :
                               Geometry::TETRAHEDRON;
   BDM_FECollection fec(p, dim);
   const FiniteElement &fe = *fec.FiniteElementForGeometry(geom);
   IsoparametricTransformation T;
   T.SetIdentityTransformation(geom);
   DenseMatrix interpolation, restriction(fe.GetDof()), product(fe.GetDof());
   fe.GetLocalInterpolation(T, interpolation);
   for (int i = 0; i < fe.GetDof(); i++) { interpolation(i,i) -= 1.0; }
   REQUIRE(interpolation.MaxMaxNorm() < 2e-10);

   // The image contains the reference element, so every restriction row is
   // defined. The two Piola transfers must be inverses on the full space.
   T.GetPointMat() *= 2.0;
   T.Reset();
   fe.GetLocalInterpolation(T, interpolation);
   fe.GetLocalRestriction(T, restriction);
   Mult(restriction, interpolation, product);
   for (int i = 0; i < fe.GetDof(); i++) { product(i,i) -= 1.0; }
   REQUIRE(product.MaxMaxNorm() < 2e-8);

   // Increasing the polynomial order must preserve a coarse polynomial.
   std::unique_ptr<FiniteElementCollection> high(fec.Clone(p+1));
   const FiniteElement &fine = *high->FiniteElementForGeometry(geom);
   T.SetIdentityTransformation(geom);
   DenseMatrix prolongation, reduction, round_trip(fe.GetDof());
   fine.GetTransferMatrix(fe, T, prolongation);
   fe.GetTransferMatrix(fine, T, reduction);
   Mult(reduction, prolongation, round_trip);
   for (int i = 0; i < fe.GetDof(); i++) { round_trip(i,i) -= 1.0; }
   REQUIRE(round_trip.MaxMaxNorm() < 2e-9);
}

TEST_CASE("BDM refinement transfer", "[BDM]")
{
   const int dim = GENERATE(2, 3);
   const int p = GENERATE(1, 2, 3, 4);
   CAPTURE(dim, p);
   Mesh mesh = MakeBDMMesh(dim);
   BDM_FECollection fec(p, dim);
   FiniteElementSpace fes(&mesh, &fec);
   GridFunction field(&fes);
   VectorFunctionCoefficient coefficient(dim,
   [=](const Vector &x, Vector &v)
   {
      for (int d = 0; d < dim; d++)
      {
         v[d] = std::pow(x[d] + 0.2*x[(d+1)%dim], p);
      }
   });
   field.ProjectCoefficient(coefficient);
   mesh.UniformRefinement();
   fes.Update();
   field.Update();
   REQUIRE(field.ComputeL2Error(coefficient) < 2e-10);
}

#ifdef MFEM_USE_EXCEPTIONS
TEST_CASE("BDM requires a nodal basis", "[BDM]")
{
   for (int basis : {BasisType::IntegratedGLL, BasisType::Positive})
   {
      REQUIRE_THROWS(BDM_TriangleElement(2, basis));
      REQUIRE_THROWS(BDM_TetrahedronElement(2, basis));
      REQUIRE_THROWS(BDM_FECollection(2, 2, basis));
      REQUIRE_THROWS(BDM_FECollection(2, 3, basis));
   }
}
#endif

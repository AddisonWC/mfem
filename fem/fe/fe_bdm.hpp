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

#ifndef MFEM_FE_BDM
#define MFEM_FE_BDM

#include "fe_base.hpp"

namespace mfem
{

/// Arbitrary order nodal Brezzi-Douglas-Marini elements on a triangle.
/** The polynomial space is [P_p]^2. Project() uses point samples and need not
    commute with divergence and scalar L2 projection. */
class BDM_TriangleElement : public VectorFiniteElement
{
   static const real_t nk[14];

#ifndef MFEM_THREAD_SAFE
   mutable Vector shape_x, shape_y, shape_l;
   mutable Vector dshape_x, dshape_y, dshape_l;
   mutable DenseMatrix u;
   mutable Vector divu;
#endif
   Array<int> dof2nk;
   DenseMatrixInverse Ti;

public:
   /** @brief Construct the BDM element of polynomial degree @a p >= 1.
       @a ob_type selects the open nodal basis for normal edge samples.
       Cell-local samples use Gauss-Legendre points. */
   BDM_TriangleElement(const int p,
                       const int ob_type = BasisType::GaussLegendre);
   void CalcVShape(const IntegrationPoint &ip,
                   DenseMatrix &shape) const override;
   void CalcVShape(ElementTransformation &Trans,
                   DenseMatrix &shape) const override
   { CalcVShape_RT(Trans, shape); }
   void CalcDivShape(const IntegrationPoint &ip,
                     Vector &divshape) const override;
   void GetLocalInterpolation(ElementTransformation &Trans,
                              DenseMatrix &I) const override
   { LocalInterpolation_RT(*this, nk, dof2nk, Trans, I); }
   void GetLocalRestriction(ElementTransformation &Trans,
                            DenseMatrix &R) const override
   { LocalRestriction_RT(nk, dof2nk, Trans, R); }
   void GetTransferMatrix(const FiniteElement &fe,
                          ElementTransformation &Trans,
                          DenseMatrix &I) const override
   { LocalInterpolation_RT(CheckVectorFE(fe), nk, dof2nk, Trans, I); }
   using FiniteElement::Project;
   void Project(VectorCoefficient &vc,
                ElementTransformation &Trans, Vector &dofs) const override
   { Project_RT(nk, dof2nk, vc, Trans, dofs); }
   void ProjectFromNodes(Vector &vc, ElementTransformation &Trans,
                         Vector &dofs) const override
   { Project_RT(nk, dof2nk, vc, Trans, dofs); }
   void ProjectMatrixCoefficient(MatrixCoefficient &mc,
                                 ElementTransformation &T,
                                 Vector &dofs) const override
   { ProjectMatrixCoefficient_RT(nk, dof2nk, mc, T, dofs); }
   void Project(const FiniteElement &fe, ElementTransformation &Trans,
                DenseMatrix &I) const override
   { Project_RT(nk, dof2nk, fe, Trans, I); }
   void ProjectCurl(const FiniteElement &fe,
                    ElementTransformation &Trans,
                    DenseMatrix &curl) const override
   { ProjectCurl2D_RT(nk, dof2nk, fe, Trans, curl); }
};

/// Arbitrary order nodal Brezzi-Douglas-Marini elements on a tetrahedron.
/** The polynomial space is [P_p]^3. Project() uses point samples and need not
    commute with divergence and scalar L2 projection. */
class BDM_TetrahedronElement : public VectorFiniteElement
{
   static const real_t nk[30];

#ifndef MFEM_THREAD_SAFE
   mutable Vector shape_x, shape_y, shape_z, shape_l;
   mutable Vector dshape_x, dshape_y, dshape_z, dshape_l;
   mutable DenseMatrix u;
   mutable Vector divu;
#endif
   Array<int> dof2nk;
   DenseMatrixInverse Ti;

public:
   /** @brief Construct the BDM element of polynomial degree @a p >= 1.
       @a ob_type selects the open nodal basis for normal face samples.
       Cell-local samples use Gauss-Legendre points. */
   BDM_TetrahedronElement(const int p,
                          const int ob_type = BasisType::GaussLegendre);
   void CalcVShape(const IntegrationPoint &ip,
                   DenseMatrix &shape) const override;
   void CalcVShape(ElementTransformation &Trans,
                   DenseMatrix &shape) const override
   { CalcVShape_RT(Trans, shape); }
   void CalcDivShape(const IntegrationPoint &ip,
                     Vector &divshape) const override;
   void GetLocalInterpolation(ElementTransformation &Trans,
                              DenseMatrix &I) const override
   { LocalInterpolation_RT(*this, nk, dof2nk, Trans, I); }
   void GetLocalRestriction(ElementTransformation &Trans,
                            DenseMatrix &R) const override
   { LocalRestriction_RT(nk, dof2nk, Trans, R); }
   void GetTransferMatrix(const FiniteElement &fe,
                          ElementTransformation &Trans,
                          DenseMatrix &I) const override
   { LocalInterpolation_RT(CheckVectorFE(fe), nk, dof2nk, Trans, I); }
   using FiniteElement::Project;
   void Project(VectorCoefficient &vc,
                ElementTransformation &Trans, Vector &dofs) const override
   { Project_RT(nk, dof2nk, vc, Trans, dofs); }
   void ProjectFromNodes(Vector &vc, ElementTransformation &Trans,
                         Vector &dofs) const override
   { Project_RT(nk, dof2nk, vc, Trans, dofs); }
   void ProjectMatrixCoefficient(MatrixCoefficient &mc,
                                 ElementTransformation &T,
                                 Vector &dofs) const override
   { ProjectMatrixCoefficient_RT(nk, dof2nk, mc, T, dofs); }
   void Project(const FiniteElement &fe, ElementTransformation &Trans,
                DenseMatrix &I) const override
   { Project_RT(nk, dof2nk, fe, Trans, I); }
   void ProjectCurl(const FiniteElement &fe,
                    ElementTransformation &Trans,
                    DenseMatrix &curl) const override
   { ProjectCurl3D_RT(nk, dof2nk, fe, Trans, curl); }
};

} // namespace mfem

#endif

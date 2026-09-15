// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.

#ifndef MFEM_FE_JM
#define MFEM_FE_JM

#include "fe_base.hpp"

namespace mfem
{

/// Coordinate choices for the same Johnson--Mercier space.
enum class JMBasis { Moments, SplitVertex };

/** The lowest-order, two-dimensional Johnson--Mercier element.

    The element consists of symmetric, piecewise-linear matrix fields on the
    Alfeld split of a triangle. Both bases use four edge and three interior
    DOFs. Moments uses traction moments and reference-pullback cell moments.
    SplitVertex uses endpoint (nn,nt) traction values and physical tensor
    components at the barycenter. The latter basis splits into four vertex
    subspaces of dimensions 4,4,4,3. Matrix fields use the double contravariant
    Piola transformation with a geometry-dependent basis correction.
    Physical mappings assume an affine element transformation. */
class JohnsonMercierTriangleFiniteElement : public FiniteElement
{
private:
   // Rows contain the coefficients of the 15 nodal basis functions in the
   // 27-dimensional broken, symmetric P1 basis on the Alfeld split.
   DenseMatrix basis;
   JMBasis basis_type;
   DenseMatrix reference_change;

   static int GetSubTriangle(const IntegrationPoint &ip);
   static void CalcRawShape(const IntegrationPoint &ip, Vector &raw);
   static void CalcRawDivShape(const IntegrationPoint &ip,
                               DenseMatrix &raw_div);
   void GetFacetTransform(const DenseMatrix &J, DenseMatrix &A) const;
   void GetMomentToSplitVertexMatrix(const DenseMatrix &J,
                                     DenseMatrix &change) const;

   void CalcShape(const IntegrationPoint &, Vector &) const override
   { MFEM_ABORT("Johnson-Mercier shape functions are matrix-valued"); }
   void CalcDShape(const IntegrationPoint &, DenseMatrix &) const override
   { MFEM_ABORT("use CalcDivShape for the Johnson-Mercier element"); }

public:
   explicit JohnsonMercierTriangleFiniteElement(
      JMBasis type = JMBasis::Moments);

   JMBasis GetBasisType() const { return basis_type; }

   /** Map physical moment coefficients to split-vertex coefficients.
       Each edge has (nn,nt) values at its first and second endpoints;
       the interior coefficients are physical (s00,s01,s11) at the barycenter.
       The matrix uses local edge orientations and is independent of the
       selected basis. Only affine triangles are supported. */
   void GetMomentToSplitVertexMatrix(ElementTransformation &Trans,
                                     DenseMatrix &change) const;

   void CalcMShape(const IntegrationPoint &ip,
                   DenseTensor &shape) const override;
   void CalcMShape(ElementTransformation &Trans,
                   DenseTensor &shape) const override;

   void CalcDivShape(const IntegrationPoint &ip,
                     DenseMatrix &divshape) const override;
   void CalcPhysDivShape(ElementTransformation &Trans,
                         DenseMatrix &divshape) const override;

   /** Project a scalar H1 element with three byVDIM components representing
       (s00,s01,s11) using the canonical Johnson--Mercier moment interpolant,
       expressed in the selected basis. In SplitVertex mode this does not
       replace moment interpolation by point interpolation of the source. */
   void Project(const FiniteElement &fe, ElementTransformation &Trans,
                DenseMatrix &I) const override;

   void GetTransferMatrix(const FiniteElement &fe,
                          ElementTransformation &Trans,
                          DenseMatrix &I) const override;

   /// Convert reference-child interpolation to the physical moment bases.
   void GetPhysicalTransferMatrix(const DenseMatrix &reference_transfer,
                                  ElementTransformation &child,
                                  ElementTransformation &fine,
                                  DenseMatrix &I) const;

   void GetLocalInterpolation(ElementTransformation &Trans,
                              DenseMatrix &I) const override
   { GetTransferMatrix(*this, Trans, I); }

   void GetFaceDofs(int face, int **dofs, int *ndofs) const override;
};

} // namespace mfem

#endif

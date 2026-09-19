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

TEST_CASE("H1 Segment Finite Element",
          "[H1_SegmentElement]"
          "[NodalFiniteElement]"
          "[ScalarFiniteElement]"
          "[FiniteElement]")
{
   for (int p = 1; p<=2; p++)
   {
      H1_SegmentElement fe(p);

      if (p == 1)
      {
         SECTION("Attributes")
         {
            REQUIRE( fe.GetDim()            == 1                           );
            REQUIRE( fe.GetGeomType()       == Geometry::SEGMENT           );
            REQUIRE( fe.Space()             == (int) FunctionSpace::Pk     );
            REQUIRE( fe.GetRangeType()      == (int) FiniteElement::SCALAR );
            REQUIRE( fe.GetMapType()        == (int) FiniteElement::VALUE  );
            REQUIRE( fe.GetDerivType()      == (int) FiniteElement::GRAD   );
            REQUIRE( fe.GetDerivRangeType() == (int) FiniteElement::VECTOR );
            REQUIRE( fe.GetDerivMapType()   == (int) FiniteElement::H_CURL );
         }
      }
      SECTION("Sizes for p = " + std::to_string(p))
      {
         REQUIRE( fe.GetDof()   == p+1 );
         REQUIRE( fe.GetOrder() == p   );
      }
   }
}

TEST_CASE("H1 Triangle Finite Element",
          "[H1_TriangleElement]"
          "[NodalFiniteElement]"
          "[ScalarFiniteElement]"
          "[FiniteElement]")
{
   for (int p = 1; p<=3; p++)
   {
      H1_TriangleElement fe(p);

      if (p == 1)
      {
         SECTION("Attributes")
         {
            REQUIRE( fe.GetDim()            == 2                           );
            REQUIRE( fe.GetGeomType()       == Geometry::TRIANGLE          );
            REQUIRE( fe.Space()             == (int) FunctionSpace::Pk     );
            REQUIRE( fe.GetRangeType()      == (int) FiniteElement::SCALAR );
            REQUIRE( fe.GetMapType()        == (int) FiniteElement::VALUE  );
            REQUIRE( fe.GetDerivType()      == (int) FiniteElement::GRAD   );
            REQUIRE( fe.GetDerivRangeType() == (int) FiniteElement::VECTOR );
            REQUIRE( fe.GetDerivMapType()   == (int) FiniteElement::H_CURL );
         }
      }
      SECTION("Sizes for p = " + std::to_string(p))
      {
         REQUIRE( fe.GetDof()   == (p+1)*(p+2)/2 );
         REQUIRE( fe.GetOrder() == p             );
      }
   }
}

TEST_CASE("H1 Quadrilateral Finite Element",
          "[H1_QuadrilateralElement]"
          "[NodalFiniteElement]"
          "[ScalarFiniteElement]"
          "[FiniteElement]")
{
   for (int p = 1; p<=3; p++)
   {
      H1_QuadrilateralElement fe(p);

      if (p == 1)
      {
         SECTION("Attributes")
         {
            REQUIRE( fe.GetDim()            == 2                           );
            REQUIRE( fe.GetGeomType()       == Geometry::SQUARE            );
            REQUIRE( fe.Space()             == (int) FunctionSpace::Qk     );
            REQUIRE( fe.GetRangeType()      == (int) FiniteElement::SCALAR );
            REQUIRE( fe.GetMapType()        == (int) FiniteElement::VALUE  );
            REQUIRE( fe.GetDerivType()      == (int) FiniteElement::GRAD   );
            REQUIRE( fe.GetDerivRangeType() == (int) FiniteElement::VECTOR );
            REQUIRE( fe.GetDerivMapType()   == (int) FiniteElement::H_CURL );
         }
      }
      SECTION("Sizes for p = " + std::to_string(p))
      {
         REQUIRE( fe.GetDof()   == (int)pow(p+1,2) );
         REQUIRE( fe.GetOrder() == p               );
      }
   }
}

TEST_CASE("H1 Tetrahedron Finite Element",
          "[H1_TetrahedronElement]"
          "[NodalFiniteElement]"
          "[ScalarFiniteElement]"
          "[FiniteElement]")
{
   for (int p = 1; p<=4; p++)
   {
      H1_TetrahedronElement fe(p);

      if (p == 1)
      {
         SECTION("Attributes")
         {
            REQUIRE( fe.GetDim()            == 3                           );
            REQUIRE( fe.GetGeomType()       == Geometry::TETRAHEDRON       );
            REQUIRE( fe.Space()             == (int) FunctionSpace::Pk     );
            REQUIRE( fe.GetRangeType()      == (int) FiniteElement::SCALAR );
            REQUIRE( fe.GetMapType()        == (int) FiniteElement::VALUE  );
            REQUIRE( fe.GetDerivType()      == (int) FiniteElement::GRAD   );
            REQUIRE( fe.GetDerivRangeType() == (int) FiniteElement::VECTOR );
            REQUIRE( fe.GetDerivMapType()   == (int) FiniteElement::H_CURL );
         }
      }
      SECTION("Sizes for p = " + std::to_string(p))
      {
         REQUIRE( fe.GetDof()   == (p+1)*(p+2)*(p+3)/6 );
         REQUIRE( fe.GetOrder() == p                   );
      }
   }
}

TEST_CASE("H1 Hexahedron Finite Element",
          "[H1_HexahedronElement]"
          "[NodalFiniteElement]"
          "[ScalarFiniteElement]"
          "[FiniteElement]")
{
   for (int p = 1; p<=4; p++)
   {
      H1_HexahedronElement fe(p);

      if (p == 1)
      {
         SECTION("Attributes")
         {
            REQUIRE( fe.GetDim()            == 3                           );
            REQUIRE( fe.GetGeomType()       == Geometry::CUBE              );
            REQUIRE( fe.Space()             == (int) FunctionSpace::Qk     );
            REQUIRE( fe.GetRangeType()      == (int) FiniteElement::SCALAR );
            REQUIRE( fe.GetMapType()        == (int) FiniteElement::VALUE  );
            REQUIRE( fe.GetDerivType()      == (int) FiniteElement::GRAD   );
            REQUIRE( fe.GetDerivRangeType() == (int) FiniteElement::VECTOR );
            REQUIRE( fe.GetDerivMapType()   == (int) FiniteElement::H_CURL );
         }
      }
      SECTION("Sizes for p = " + std::to_string(p))
      {
         REQUIRE( fe.GetDof()   == (int)pow(p+1,3) );
         REQUIRE( fe.GetOrder() == p               );
      }
   }
}

TEST_CASE("H1 Wedge Finite Element",
          "[H1_WedgeElement]"
          "[NodalFiniteElement]"
          "[ScalarFiniteElement]"
          "[FiniteElement]")
{
   for (int p = 1; p<=4; p++)
   {
      H1_WedgeElement fe(p);

      if (p == 1)
      {
         SECTION("Attributes")
         {
            REQUIRE( fe.GetDim()            == 3                           );
            REQUIRE( fe.GetGeomType()       == Geometry::PRISM             );
            REQUIRE( fe.Space()             == (int) FunctionSpace::Qk     );
            REQUIRE( fe.GetRangeType()      == (int) FiniteElement::SCALAR );
            REQUIRE( fe.GetMapType()        == (int) FiniteElement::VALUE  );
            REQUIRE( fe.GetDerivType()      == (int) FiniteElement::GRAD   );
            REQUIRE( fe.GetDerivRangeType() == (int) FiniteElement::VECTOR );
            REQUIRE( fe.GetDerivMapType()   == (int) FiniteElement::H_CURL );
         }
      }
      SECTION("Sizes for p = " + std::to_string(p))
      {
         REQUIRE( fe.GetDof()   == (p+1)*(p+1)*(p+2)/2 );
         REQUIRE( fe.GetOrder() == p                   );
      }
   }
}

TEST_CASE("H1 Pyramid Finite Element",
          "[H1_PyramidElement]"
          "[NodalFiniteElement]"
          "[ScalarFiniteElement]"
          "[FiniteElement]")
{
   for (int p = 1; p<=4; p++)
   {
      H1_FuentesPyramidElement fe(p);

      if (p == 1)
      {
         SECTION("Attributes")
         {
            REQUIRE( fe.GetDim()            == 3                           );
            REQUIRE( fe.GetGeomType()       == Geometry::PYRAMID           );
            REQUIRE( fe.Space()             == (int) FunctionSpace::Uk     );
            REQUIRE( fe.GetRangeType()      == (int) FiniteElement::SCALAR );
            REQUIRE( fe.GetMapType()        == (int) FiniteElement::VALUE  );
            REQUIRE( fe.GetDerivType()      == (int) FiniteElement::GRAD   );
            REQUIRE( fe.GetDerivRangeType() == (int) FiniteElement::VECTOR );
            REQUIRE( fe.GetDerivMapType()   == (int) FiniteElement::H_CURL );
         }
      }
      SECTION("Sizes for p = " + std::to_string(p))
      {
         REQUIRE( fe.GetDof()   == p*(p*p+3)+1 ); // Fuentes et al
         REQUIRE( fe.GetOrder() == p                   );
      }
   }
   for (int p = 1; p<=4; p++)
   {
      H1_BergotPyramidElement fe(p);

      if (p == 1)
      {
         SECTION("Attributes")
         {
            REQUIRE( fe.GetDim()            == 3                           );
            REQUIRE( fe.GetGeomType()       == Geometry::PYRAMID           );
            REQUIRE( fe.Space()             == (int) FunctionSpace::Uk     );
            REQUIRE( fe.GetRangeType()      == (int) FiniteElement::SCALAR );
            REQUIRE( fe.GetMapType()        == (int) FiniteElement::VALUE  );
            REQUIRE( fe.GetDerivType()      == (int) FiniteElement::GRAD   );
            REQUIRE( fe.GetDerivRangeType() == (int) FiniteElement::VECTOR );
            REQUIRE( fe.GetDerivMapType()   == (int) FiniteElement::H_CURL );
         }
      }
      SECTION("Sizes for p = " + std::to_string(p))
      {
         REQUIRE( fe.GetDof()   == (p+1)*(p+2)*(2*p+3)/6 ); // JSC
         REQUIRE( fe.GetOrder() == p                   );
      }
   }
}

TEST_CASE("Nedelec Segment Finite Element",
          "[ND_SegmentElement]"
          "[VectorFiniteElement]"
          "[FiniteElement]")
{
   for (int p = 1; p<=2; p++)
   {
      ND_SegmentElement fe(p);

      if (p == 1)
      {
         SECTION("Attributes")
         {
            REQUIRE( fe.GetDim()            == 1                            );
            REQUIRE( fe.GetGeomType()       == Geometry::SEGMENT            );
            REQUIRE( fe.Space()             == (int) FunctionSpace::Pk      );
            REQUIRE( fe.GetRangeType()      == (int) FiniteElement::VECTOR  );
            REQUIRE( fe.GetMapType()        == (int) FiniteElement::H_CURL  );
            REQUIRE( fe.GetDerivType()      == (int) FiniteElement::NONE    );
            REQUIRE( fe.GetDerivRangeType() ==
                     (int) FiniteElement::UNKNOWN_RANGE_TYPE);
            REQUIRE( fe.GetDerivMapType()   ==
                     (int) FiniteElement::UNKNOWN_MAP_TYPE);
         }
      }
      SECTION("Sizes for p = " + std::to_string(p))
      {
         REQUIRE( fe.GetDof()   == p   );
         REQUIRE( fe.GetOrder() == p-1 );
      }
   }
}

TEST_CASE("Nedelec Triangular Finite Element",
          "[ND_TriangleElement]"
          "[VectorFiniteElement]"
          "[FiniteElement]")
{
   for (int p = 1; p<=3; p++)
   {
      ND_TriangleElement fe(p);

      if (p == 1)
      {
         SECTION("Attributes")
         {
            REQUIRE( fe.GetDim()            == 2                            );
            REQUIRE( fe.GetGeomType()       == Geometry::TRIANGLE           );
            REQUIRE( fe.Space()             == (int) FunctionSpace::Pk      );
            REQUIRE( fe.GetRangeType()      == (int) FiniteElement::VECTOR  );
            REQUIRE( fe.GetMapType()        == (int) FiniteElement::H_CURL  );
            REQUIRE( fe.GetDerivType()      == (int) FiniteElement::CURL    );
            REQUIRE( fe.GetDerivRangeType() == (int) FiniteElement::SCALAR  );
            REQUIRE( fe.GetDerivMapType()   == (int) FiniteElement::INTEGRAL);
         }
      }
      SECTION("Sizes for p = " + std::to_string(p))
      {
         REQUIRE( fe.GetDof()   == p*(p+2) );
         REQUIRE( fe.GetOrder() == p       );
      }
   }
}

TEST_CASE("Nedelec Quadrilateral Finite Element",
          "[ND_QuadrilateralElement]"
          "[VectorFiniteElement]"
          "[FiniteElement]")
{
   for (int p = 1; p<=3; p++)
   {
      ND_QuadrilateralElement fe(p);

      if (p == 1)
      {
         SECTION("Attributes")
         {
            REQUIRE( fe.GetDim()            == 2                            );
            REQUIRE( fe.GetGeomType()       == Geometry::SQUARE             );
            REQUIRE( fe.Space()             == (int) FunctionSpace::Qk      );
            REQUIRE( fe.GetRangeType()      == (int) FiniteElement::VECTOR  );
            REQUIRE( fe.GetMapType()        == (int) FiniteElement::H_CURL  );
            REQUIRE( fe.GetDerivType()      == (int) FiniteElement::CURL    );
            REQUIRE( fe.GetDerivRangeType() == (int) FiniteElement::SCALAR  );
            REQUIRE( fe.GetDerivMapType()   == (int) FiniteElement::INTEGRAL);
         }
      }
      SECTION("Sizes for p = " + std::to_string(p))
      {
         REQUIRE( fe.GetDof()   == 2*p*(p+1) );
         REQUIRE( fe.GetOrder() == p         );
      }
   }
}

TEST_CASE("Nedelec Tetrahedron Finite Element",
          "[ND_TetrahedronElement]"
          "[VectorFiniteElement]"
          "[FiniteElement]")
{
   for (int p = 1; p<=4; p++)
   {
      ND_TetrahedronElement fe(p);

      if (p == 1)
      {
         SECTION("Attributes")
         {
            REQUIRE( fe.GetDim()            == 3                           );
            REQUIRE( fe.GetGeomType()       == Geometry::TETRAHEDRON       );
            REQUIRE( fe.Space()             == (int) FunctionSpace::Pk     );
            REQUIRE( fe.GetRangeType()      == (int) FiniteElement::VECTOR );
            REQUIRE( fe.GetMapType()        == (int) FiniteElement::H_CURL );
            REQUIRE( fe.GetDerivType()      == (int) FiniteElement::CURL   );
            REQUIRE( fe.GetDerivRangeType() == (int) FiniteElement::VECTOR );
            REQUIRE( fe.GetDerivMapType()   == (int) FiniteElement::H_DIV  );
         }
      }
      SECTION("Sizes for p = " + std::to_string(p))
      {
         REQUIRE( fe.GetDof()   == p*(p+2)*(p+3)/2 );
         REQUIRE( fe.GetOrder() == p               );
      }
   }
}

TEST_CASE("Nedelec Hexahedron Finite Element",
          "[ND_HexahedronElement]"
          "[VectorFiniteElement]"
          "[FiniteElement]")
{
   for (int p = 1; p<=4; p++)
   {
      ND_HexahedronElement fe(p);

      if (p == 1)
      {
         SECTION("Attributes")
         {
            REQUIRE( fe.GetDim()            == 3                           );
            REQUIRE( fe.GetGeomType()       == Geometry::CUBE              );
            REQUIRE( fe.Space()             == (int) FunctionSpace::Qk     );
            REQUIRE( fe.GetRangeType()      == (int) FiniteElement::VECTOR );
            REQUIRE( fe.GetMapType()        == (int) FiniteElement::H_CURL );
            REQUIRE( fe.GetDerivType()      == (int) FiniteElement::CURL   );
            REQUIRE( fe.GetDerivRangeType() == (int) FiniteElement::VECTOR );
            REQUIRE( fe.GetDerivMapType()   == (int) FiniteElement::H_DIV  );
         }
      }
      SECTION("Sizes for p = " + std::to_string(p))
      {
         REQUIRE( fe.GetDof()   == 3*p*(int)pow(p+1,2) );
         REQUIRE( fe.GetOrder() == p                   );
      }
   }
}

TEST_CASE("Raviart-Thomas Triangular Finite Element",
          "[RT_TriangleElement]"
          "[VectorFiniteElement]"
          "[FiniteElement]")
{
   for (int p = 1; p<=3; p++)
   {
      RT_TriangleElement fe(p-1);

      if (p == 1)
      {
         SECTION("Attributes")
         {
            REQUIRE( fe.GetDim()            == 2                             );
            REQUIRE( fe.GetGeomType()       == Geometry::TRIANGLE            );
            REQUIRE( fe.Space()             == (int) FunctionSpace::Pk       );
            REQUIRE( fe.GetRangeType()      == (int) FiniteElement::VECTOR   );
            REQUIRE( fe.GetMapType()        == (int) FiniteElement::H_DIV    );
            REQUIRE( fe.GetDerivType()      == (int) FiniteElement::DIV      );
            REQUIRE( fe.GetDerivRangeType() == (int) FiniteElement::SCALAR   );
            REQUIRE( fe.GetDerivMapType()   == (int) FiniteElement::INTEGRAL );
         }
      }
      SECTION("Sizes for p = " + std::to_string(p))
      {
         REQUIRE( fe.GetDof()   == p*(p+2) );
         REQUIRE( fe.GetOrder() == p       );
      }
   }
}

TEST_CASE("Raviart-Thomas Quadrilateral Finite Element",
          "[RT_QuadrilateralElement]"
          "[VectorFiniteElement]"
          "[FiniteElement]")
{
   for (int p = 1; p<=3; p++)
   {
      RT_QuadrilateralElement fe(p-1);

      if (p == 1)
      {
         SECTION("Attributes")
         {
            REQUIRE( fe.GetDim()            == 2                             );
            REQUIRE( fe.GetGeomType()       == Geometry::SQUARE              );
            REQUIRE( fe.Space()             == (int) FunctionSpace::Qk       );
            REQUIRE( fe.GetRangeType()      == (int) FiniteElement::VECTOR   );
            REQUIRE( fe.GetMapType()        == (int) FiniteElement::H_DIV    );
            REQUIRE( fe.GetDerivType()      == (int) FiniteElement::DIV      );
            REQUIRE( fe.GetDerivRangeType() == (int) FiniteElement::SCALAR   );
            REQUIRE( fe.GetDerivMapType()   == (int) FiniteElement::INTEGRAL );
         }
      }
      SECTION("Sizes for p = " + std::to_string(p))
      {
         REQUIRE( fe.GetDof()   == 2*p*(p+1) );
         REQUIRE( fe.GetOrder() == p         );
      }
   }
}

TEST_CASE("Raviart-Thomas Tetrahedron Finite Element",
          "[RT_TetrahedronElement]"
          "[VectorFiniteElement]"
          "[FiniteElement]")
{
   for (int p = 1; p<=4; p++)
   {
      RT_TetrahedronElement fe(p-1);

      if (p == 1)
      {
         SECTION("Attributes")
         {
            REQUIRE( fe.GetDim()            == 3                             );
            REQUIRE( fe.GetGeomType()       == Geometry::TETRAHEDRON         );
            REQUIRE( fe.Space()             == (int) FunctionSpace::Pk       );
            REQUIRE( fe.GetRangeType()      == (int) FiniteElement::VECTOR   );
            REQUIRE( fe.GetMapType()        == (int) FiniteElement::H_DIV    );
            REQUIRE( fe.GetDerivType()      == (int) FiniteElement::DIV      );
            REQUIRE( fe.GetDerivRangeType() == (int) FiniteElement::SCALAR   );
            REQUIRE( fe.GetDerivMapType()   == (int) FiniteElement::INTEGRAL );
         }
      }
      SECTION("Sizes for p = " + std::to_string(p))
      {
         REQUIRE( fe.GetDof()   == p*(p+1)*(p+3)/2 );
         REQUIRE( fe.GetOrder() == p               );
      }
   }
}

TEST_CASE("Raviart-Thomas Hexahedron Finite Element",
          "[RT_HexahedronElement]"
          "[VectorFiniteElement]"
          "[FiniteElement]")
{
   for (int p = 1; p<=4; p++)
   {
      RT_HexahedronElement fe(p-1);

      if (p == 1)
      {
         SECTION("Attributes")
         {
            REQUIRE( fe.GetDim()            == 3                             );
            REQUIRE( fe.GetGeomType()       == Geometry::CUBE                );
            REQUIRE( fe.Space()             == (int) FunctionSpace::Qk       );
            REQUIRE( fe.GetRangeType()      == (int) FiniteElement::VECTOR   );
            REQUIRE( fe.GetMapType()        == (int) FiniteElement::H_DIV    );
            REQUIRE( fe.GetDerivType()      == (int) FiniteElement::DIV      );
            REQUIRE( fe.GetDerivRangeType() == (int) FiniteElement::SCALAR   );
            REQUIRE( fe.GetDerivMapType()   == (int) FiniteElement::INTEGRAL );
         }
      }
      SECTION("Sizes for p = " + std::to_string(p))
      {
         REQUIRE( fe.GetDof()   == 3*(p+1)*(int)pow(p,2) );
         REQUIRE( fe.GetOrder() == p                     );
      }
   }
}

TEST_CASE("Brezzi-Douglas-Marini Simplex Finite Elements",
          "[BDM]"
          "[BDM_TriangleElement]"
          "[BDM_TetrahedronElement]"
          "[VectorFiniteElement]"
          "[FiniteElement]")
{
   for (int p = 1; p <= 4; p++)
   {
      BDM_TriangleElement tri(p);
      BDM_TetrahedronElement tet(p);

      REQUIRE(tri.GetDim() == 2);
      REQUIRE(tri.GetGeomType() == Geometry::TRIANGLE);
      REQUIRE(tri.Space() == (int) FunctionSpace::Pk);
      REQUIRE(tri.GetRangeType() == (int) FiniteElement::VECTOR);
      REQUIRE(tri.GetMapType() == (int) FiniteElement::H_DIV);
      REQUIRE(tri.GetDerivType() == (int) FiniteElement::DIV);
      REQUIRE(tri.GetDof() == (p + 1)*(p + 2));
      REQUIRE(tri.GetOrder() == p);

      REQUIRE(tet.GetDim() == 3);
      REQUIRE(tet.GetGeomType() == Geometry::TETRAHEDRON);
      REQUIRE(tet.Space() == (int) FunctionSpace::Pk);
      REQUIRE(tet.GetRangeType() == (int) FiniteElement::VECTOR);
      REQUIRE(tet.GetMapType() == (int) FiniteElement::H_DIV);
      REQUIRE(tet.GetDerivType() == (int) FiniteElement::DIV);
      REQUIRE(tet.GetDof() == (p + 1)*(p + 2)*(p + 3)/2);
      REQUIRE(tet.GetOrder() == p);
   }
}

TEST_CASE("Brezzi-Douglas-Marini Finite Element Collection",
          "[BDM]"
          "[BDM_FECollection]"
          "[FiniteElementCollection]")
{
   const int basis = GENERATE(BasisType::GaussLegendre, BasisType::OpenUniform,
                              BasisType::OpenHalfUniform);
   for (int dim = 2; dim <= 3; dim++)
   {
      for (int p = 1; p <= 4; p++)
      {
         CAPTURE(dim, p, basis);
         BDM_FECollection fec(p, dim, basis);
         const FiniteElementCollection *base_fec = &fec;
         REQUIRE(dynamic_cast<const RT_FECollection *>(base_fec) == nullptr);
         REQUIRE(fec.GetOrder() == p);
         REQUIRE(fec.GetConstructorOrder() == p);
         REQUIRE(fec.GetContType() == FiniteElementCollection::NORMAL);
         REQUIRE(fec.GetOpenBasisType() == basis);
         REQUIRE(fec.DofForGeometry(dim == 2 ? Geometry::SEGMENT :
                                    Geometry::TRIANGLE) ==
                 (dim == 2 ? p + 1 : (p + 1)*(p + 2)/2));

         const Geometry::Type volume = dim == 2 ? Geometry::TRIANGLE :
                                       Geometry::TETRAHEDRON;
         REQUIRE(fec.FiniteElementForGeometry(volume) != nullptr);
         REQUIRE(fec.FiniteElementForGeometry(dim == 2 ? Geometry::SQUARE :
                                              Geometry::CUBE) == nullptr);
         if (dim == 3)
         {
            REQUIRE(fec.FiniteElementForGeometry(Geometry::SQUARE) == nullptr);
            REQUIRE(fec.DofForGeometry(Geometry::SQUARE) == 0);
            REQUIRE(fec.DofOrderForOrientation(Geometry::SQUARE, 0) == nullptr);
         }

         std::unique_ptr<FiniteElementCollection> copy(
            FiniteElementCollection::New(fec.Name()));
         REQUIRE(copy != nullptr);
         REQUIRE(std::string(copy->Name()) == fec.Name());
         std::unique_ptr<FiniteElementCollection> higher(fec.Clone(p+1));
         REQUIRE(higher->GetOrder() == p+1);
         REQUIRE(higher->GetConstructorOrder() == p+1);
         REQUIRE(higher->FiniteElementForGeometry(volume)->GetOrder() == p+1);
         REQUIRE(dynamic_cast<BDM_FECollection *>(higher.get()) != nullptr);

         Mesh mesh = dim == 2 ?
                     Mesh::MakeCartesian2D(2, 2, Element::TRIANGLE, true) :
                     Mesh::MakeCartesian3D(1, 1, 1, Element::TETRAHEDRON);
         FiniteElementSpace fespace(&mesh, &fec);
         const int entity_dofs = dim == 2 ? p + 1 : (p + 1)*(p + 2)/2;
         const int interior_dofs = dim == 2 ? (p + 1)*(p - 1) :
                                   (p + 1)*(p + 2)*(p - 1)/2;
         const int entities = dim == 2 ? mesh.GetNEdges() : mesh.GetNFaces();
         REQUIRE(fespace.GetVSize() ==
                 entities*entity_dofs + mesh.GetNE()*interior_dofs);
      }
   }
}

TEST_CASE("L2 Segment Finite Element",
          "[L2_SegmentElement]"
          "[NodalFiniteElement]"
          "[ScalarFiniteElement]"
          "[FiniteElement]")
{
   for (int p = 0; p<=1; p++)
   {
      L2_SegmentElement fe(p);

      if (p == 0)
      {
         SECTION("Attributes")
         {
            REQUIRE( fe.GetDim()            == 1                           );
            REQUIRE( fe.GetGeomType()       == Geometry::SEGMENT           );
            REQUIRE( fe.Space()             == (int) FunctionSpace::Pk     );
            REQUIRE( fe.GetRangeType()      == (int) FiniteElement::SCALAR );
            REQUIRE( fe.GetMapType()        == (int) FiniteElement::VALUE  );
            REQUIRE( fe.GetDerivType()      == (int) FiniteElement::GRAD   );
            REQUIRE( fe.GetDerivRangeType() == (int) FiniteElement::VECTOR );
            REQUIRE( fe.GetDerivMapType()   == (int) FiniteElement::H_CURL );
         }
      }
      SECTION("Sizes for p = " + std::to_string(p))
      {
         REQUIRE( fe.GetDof()   == p+1 );
         REQUIRE( fe.GetOrder() == p   );
      }
   }
}

TEST_CASE("L2 Triangle Finite Element",
          "[L2_TriangleElement]"
          "[NodalFiniteElement]"
          "[ScalarFiniteElement]"
          "[FiniteElement]")
{
   for (int p = 0; p<=2; p++)
   {
      L2_TriangleElement fe(p);

      if (p == 1)
      {
         SECTION("Attributes")
         {
            REQUIRE( fe.GetDim()            == 2                           );
            REQUIRE( fe.GetGeomType()       == Geometry::TRIANGLE          );
            REQUIRE( fe.Space()             == (int) FunctionSpace::Pk     );
            REQUIRE( fe.GetRangeType()      == (int) FiniteElement::SCALAR );
            REQUIRE( fe.GetMapType()        == (int) FiniteElement::VALUE  );
            REQUIRE( fe.GetDerivType()      == (int) FiniteElement::GRAD   );
            REQUIRE( fe.GetDerivRangeType() == (int) FiniteElement::VECTOR );
            REQUIRE( fe.GetDerivMapType()   == (int) FiniteElement::H_CURL );
         }
      }
      SECTION("Sizes for p = " + std::to_string(p))
      {
         REQUIRE( fe.GetDof()   == (p+1)*(p+2)/2 );
         REQUIRE( fe.GetOrder() == p             );
      }
   }
}

TEST_CASE("L2 Quadrilateral Finite Element",
          "[L2_QuadrilateralElement]"
          "[NodalFiniteElement]"
          "[ScalarFiniteElement]"
          "[FiniteElement]")
{
   for (int p = 0; p<=2; p++)
   {
      L2_QuadrilateralElement fe(p);

      if (p == 1)
      {
         SECTION("Attributes")
         {
            REQUIRE( fe.GetDim()            == 2                           );
            REQUIRE( fe.GetGeomType()       == Geometry::SQUARE            );
            REQUIRE( fe.Space()             == (int) FunctionSpace::Qk     );
            REQUIRE( fe.GetRangeType()      == (int) FiniteElement::SCALAR );
            REQUIRE( fe.GetMapType()        == (int) FiniteElement::VALUE  );
            REQUIRE( fe.GetDerivType()      == (int) FiniteElement::GRAD   );
            REQUIRE( fe.GetDerivRangeType() == (int) FiniteElement::VECTOR );
            REQUIRE( fe.GetDerivMapType()   == (int) FiniteElement::H_CURL );
         }
      }
      SECTION("Sizes for p = " + std::to_string(p))
      {
         REQUIRE( fe.GetDof()   == (int)pow(p+1,2) );
         REQUIRE( fe.GetOrder() == p               );
      }
   }
}

TEST_CASE("L2 Tetrahedron Finite Element",
          "[L2_TetrahedronElement]"
          "[NodalFiniteElement]"
          "[ScalarFiniteElement]"
          "[FiniteElement]")
{
   for (int p = 0; p<=3; p++)
   {
      L2_TetrahedronElement fe(p);

      if (p == 1)
      {
         SECTION("Attributes")
         {
            REQUIRE( fe.GetDim()            == 3                           );
            REQUIRE( fe.GetGeomType()       == Geometry::TETRAHEDRON       );
            REQUIRE( fe.Space()             == (int) FunctionSpace::Pk     );
            REQUIRE( fe.GetRangeType()      == (int) FiniteElement::SCALAR );
            REQUIRE( fe.GetMapType()        == (int) FiniteElement::VALUE  );
            REQUIRE( fe.GetDerivType()      == (int) FiniteElement::GRAD   );
            REQUIRE( fe.GetDerivRangeType() == (int) FiniteElement::VECTOR );
            REQUIRE( fe.GetDerivMapType()   == (int) FiniteElement::H_CURL );
         }
      }
      SECTION("Sizes for p = " + std::to_string(p))
      {
         REQUIRE( fe.GetDof()   == (p+1)*(p+2)*(p+3)/6 );
         REQUIRE( fe.GetOrder() == p                   );
      }
   }
}

TEST_CASE("L2 Hexahedron Finite Element",
          "[L2_HexahedronElement]"
          "[NodalFiniteElement]"
          "[ScalarFiniteElement]"
          "[FiniteElement]")
{
   for (int p = 0; p<=3; p++)
   {
      L2_HexahedronElement fe(p);

      if (p == 1)
      {
         SECTION("Attributes")
         {
            REQUIRE( fe.GetDim()            == 3                           );
            REQUIRE( fe.GetGeomType()       == Geometry::CUBE              );
            REQUIRE( fe.Space()             == (int) FunctionSpace::Qk     );
            REQUIRE( fe.GetRangeType()      == (int) FiniteElement::SCALAR );
            REQUIRE( fe.GetMapType()        == (int) FiniteElement::VALUE  );
            REQUIRE( fe.GetDerivType()      == (int) FiniteElement::GRAD   );
            REQUIRE( fe.GetDerivRangeType() == (int) FiniteElement::VECTOR );
            REQUIRE( fe.GetDerivMapType()   == (int) FiniteElement::H_CURL );
         }
      }
      SECTION("Sizes for p = " + std::to_string(p))
      {
         REQUIRE( fe.GetDof()   == (int)pow(p+1,3) );
         REQUIRE( fe.GetOrder() == p               );
      }
   }
}

TEST_CASE("L2 Wedge Finite Element",
          "[L2_WedgeElement]"
          "[NodalFiniteElement]"
          "[ScalarFiniteElement]"
          "[FiniteElement]")
{
   for (int p = 0; p<=3; p++)
   {
      L2_WedgeElement fe(p);

      if (p == 1)
      {
         SECTION("Attributes")
         {
            REQUIRE( fe.GetDim()            == 3                           );
            REQUIRE( fe.GetGeomType()       == Geometry::PRISM             );
            REQUIRE( fe.Space()             == (int) FunctionSpace::Qk     );
            REQUIRE( fe.GetRangeType()      == (int) FiniteElement::SCALAR );
            REQUIRE( fe.GetMapType()        == (int) FiniteElement::VALUE  );
            REQUIRE( fe.GetDerivType()      == (int) FiniteElement::GRAD   );
            REQUIRE( fe.GetDerivRangeType() == (int) FiniteElement::VECTOR );
            REQUIRE( fe.GetDerivMapType()   == (int) FiniteElement::H_CURL );
         }
      }
      SECTION("Sizes for p = " + std::to_string(p))
      {
         REQUIRE( fe.GetDof()   == (p+1)*(p+1)*(p+2)/2 );
         REQUIRE( fe.GetOrder() == p                   );
      }
   }
}

TEST_CASE("Nedelec Wedge Finite Element",
          "[ND_WedgeElement]"
          "[VectorFiniteElement]"
          "[FiniteElement]")
{
   int p = 1;

   ND_WedgeElement fe(p);

   SECTION("Attributes")
   {
      REQUIRE( fe.GetDim()            == 3                     );
      REQUIRE( fe.GetGeomType()       == Geometry::PRISM       );
      REQUIRE( fe.GetDof()            == 3*p*(p+1)*(p+2)/2     );
      REQUIRE( fe.GetOrder()          == p                     );
      REQUIRE( fe.Space()             == (int) FunctionSpace::Qk     );
      REQUIRE( fe.GetRangeType()      == (int) FiniteElement::VECTOR );
      REQUIRE( fe.GetMapType()        == (int) FiniteElement::H_CURL );
      REQUIRE( fe.GetDerivType()      == (int) FiniteElement::CURL   );
      REQUIRE( fe.GetDerivRangeType() == (int) FiniteElement::VECTOR );
      REQUIRE( fe.GetDerivMapType()   == (int) FiniteElement::H_DIV  );
   }
}

TEST_CASE("Raviart-Thomas Wedge Finite Element",
          "[RT_WedgeElement]"
          "[VectorFiniteElement]"
          "[FiniteElement]")
{
   int p = 1;

   RT_WedgeElement fe(p-1);

   SECTION("Attributes")
   {
      REQUIRE( fe.GetDim()            == 3                       );
      REQUIRE( fe.GetGeomType()       == Geometry::PRISM         );
      REQUIRE( fe.GetDof()            == (int)pow(p+1,2)*p/2 +
               /*                     */ (int)pow(p,2)*(p+2)     );
      REQUIRE( fe.GetOrder()          == p                       );
      REQUIRE( fe.Space()             == (int) FunctionSpace::Qk       );
      REQUIRE( fe.GetRangeType()      == (int) FiniteElement::VECTOR   );
      REQUIRE( fe.GetMapType()        == (int) FiniteElement::H_DIV    );
      REQUIRE( fe.GetDerivType()      == (int) FiniteElement::DIV      );
      REQUIRE( fe.GetDerivRangeType() == (int) FiniteElement::SCALAR   );
      REQUIRE( fe.GetDerivMapType()   == (int) FiniteElement::INTEGRAL );
   }
}

# BDM nodal DOF rationale

These discussion notes accompany the simplex BDM implementation. They are not
part of the library documentation or a prerequisite for using the elements.
The library should follow the existing RT/ND implementation and documentation
style; it does not need an embedded unisolvence proof.

## Existing scalar nodes

`L2_TriangleElement` and `L2_TetrahedronElement` in `fem/fe/fe_l2.cpp`
already use normalized barycentric grids constructed from open one-dimensional
points. Gauss–Legendre is the default. For degree m, take m+1 open points a_i
and all nonnegative multi-indices with sum m, then set

    lambda_j = a_{alpha_j} / sum_k a_{alpha_k}.

The RT interior grid is precisely the scalar L2 grid of degree p-1, with
independent component samples at each point. The proposed BDM layout borrows
the point locations and directions of MFEM's ND element of order p-1, with
all those additional DOFs owned by the cell. Normal facet DOFs remain shared.

For BDM triangles, the strictly interior component samples use the scalar
L2 grid of degree p-3. For tetrahedra, the additional face-tangential samples
use triangular L2 grids of degree p-3, and the volume component samples use
the tetrahedral L2 grid of degree p-4. Negative-degree stages are absent.

This reuses MFEM's scalar interpolation convention. We accept unisolvence of
the existing scalar L2 elements as a library premise. The BDM argument below
reduces vector uniqueness to those scalar interpolation problems; no new
scalar-node proof or node family is needed. No claim of optimal conditioning
is made.

## Uniqueness argument for the vector DOFs

Let v belong to [P_p]^d and suppose all its DOFs vanish.

On a triangle:

1. The p+1 normal samples on each edge force its degree-p normal trace to
   vanish. The independent normals at each vertex imply v=0 there.
2. The tangential trace on each edge vanishes at both endpoints and at the
   p-1 additional edge points. It is therefore identically zero.
3. Both components vanish on every edge, so
   v = lambda_0 lambda_1 lambda_2 w, with w in [P_{p-3}]^2.
4. At strictly interior nodes the bubble factor is nonzero. The component
   samples therefore determine w through the existing scalar P_{p-3} grid.
   Hence w=0. For p<3, the boundary vanishing already forces v=0.

On a tetrahedron:

1. The scalar P_p face grids determine the normal trace on each face.
2. Zero face-normal traces force v to be parallel to each edge and zero at
   vertices. The p-1 edge-tangential samples then force v=0 on every edge.
3. On each face v is tangential and zero on its edges. Each of its two
   tangential components contains the triangular bubble factor. The
   remaining degree-p-3 factors are determined by the face-interior samples.
4. Thus v vanishes on every face, and
   v = lambda_0 lambda_1 lambda_2 lambda_3 w, with w in [P_{p-4}]^3.
   The volume-interior component samples determine w. Low-degree stages
   with no possible nonzero bubble polynomial require no samples.

The number of vector DOFs equals dim [P_p]^d, so this uniqueness argument
also establishes existence of interpolation, under the scalar premise above.
Cell-owned tangential samples on the boundary do not impose tangential
continuity between cells.

## Interpolation semantics

An unisolvent nodal basis represents the same BDM polynomial space as a
canonical moment basis. Their interpolants of general nonpolynomial or
higher-degree fields need not agree. In particular, the point-based Project
operation does not automatically commute with divergence and scalar L2
projection. This distinction also matters when comparing to moment-based
implementations; it does not invalidate the conforming approximation space.

## Implementation and verification

The implementation rejects IntegratedGLL: nodal volume samples cannot share
coefficients with an integrated segment trace basis. The open nodal options
GaussLegendre, OpenUniform, and OpenHalfUniform control the normal facet
samples. Cell-local samples retain MFEM's default Gauss-Legendre layout.

`tests/unit/fem/test_bdm.cpp` covers full-degree vector polynomial reproduction
and divergence on affine meshes with shear, arbitrary-field normal traces on
shared faces, local interpolation/restriction, transfer between polynomial
orders, and mesh refinement. Collection metadata, factory names, and cloning
are covered in `tests/unit/fem/test_fe.cpp`. The common test tag is `[BDM]`.
Invalid-basis exception tests are enabled when MFEM_USE_EXCEPTIONS is enabled.

Run through the normal MFEM unit-test target, for example:

```sh
cmake --build build-hx --target unit_tests -j 6
build-hx/tests/unit/unit_tests '[BDM]'
```

Mathematical discussion remains here, with concise implementation comments
and API documentation in the library.

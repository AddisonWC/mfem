// Additive HX for symmetric H(div) stresses with sigma n = 0 on part of the
// boundary. See reports/hx_bc.md for the DOF analysis and results.
//
// Every space is restricted by a sparse matrix P whose columns span the
// constrained subspace: A_r = P^T A P, and the auxiliary maps become
// P_s^T Pi P_h1 and P_s^T Airy P_phi. The stress P_s has orthonormal
// columns, so P_s^T gives coordinates of anything already in its range.
#include "mfem.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

extern "C" void dstev_(const char*, const int*, double*, double*, double*,
                       const int*, double*, int*);
extern "C" void dsyev_(const char*, const char*, const int*, double*,
                       const int*, double*, double*, const int*, int*);
extern "C" void dsygv_(const int*, const char*, const char*, const int*,
                       double*, const int*, double*, const int*, double*,
                       double*, const int*, int*);

using namespace mfem;
using namespace std;

namespace
{
using Direction = array<real_t,2>;

/// Traction-free part Gamma_s of the boundary, as mesh edges and vertices.
struct TractionBoundary
{
   vector<int> edges, edge_component;
   vector<vector<Direction>> tangents; // unit tangents of incident Gamma_s edges
   vector<int> vertex_component;       // -1 off Gamma_s
   int components = 0;
};

// Sides of the unit square: b,r,t,l; h is the lower half of the left side;
// a is the whole boundary. Classification uses the unrotated coordinates.
TractionBoundary MarkBoundary(Mesh &mesh, const string &spec)
{
   TractionBoundary bc;
   const int nv = mesh.GetNV();
   bc.tangents.resize(nv);
   bc.vertex_component.assign(nv, -1);
   auto on = [&](char side) { return spec == "a" || spec.find(side) != string::npos; };
   Array<int> v;
   vector<int> parent(nv);
   for (int i = 0; i < nv; i++) { parent[i] = i; }
   function<int(int)> root = [&](int i) { return parent[i] == i ? i : parent[i] = root(parent[i]); };
   for (int b = 0; b < mesh.GetNBE(); b++)
   {
      mesh.GetBdrElementVertices(b, v);
      const real_t *p = mesh.GetVertex(v[0]), *q = mesh.GetVertex(v[1]);
      const real_t x = (p[0]+q[0])/2, y = (p[1]+q[1])/2, eps = 1e-10;
      const bool marked = (on('b') && y < eps) || (on('t') && y > 1-eps) ||
                          (on('r') && x > 1-eps) || (on('l') && x < eps) ||
                          (on('h') && x < eps && y < 0.5);
      if (!marked) { continue; }
      bc.edges.push_back(mesh.GetBdrElementFaceIndex(b));
      parent[root(v[0])] = root(v[1]);
      bc.vertex_component[v[0]] = bc.vertex_component[v[1]] = 0;
   }
   vector<int> label(nv, -1);
   for (int i = 0; i < nv; i++)
   {
      if (bc.vertex_component[i] < 0) { continue; }
      int &l = label[root(i)];
      if (l < 0) { l = bc.components++; }
      bc.vertex_component[i] = l;
   }
   for (int e : bc.edges)
   {
      mesh.GetEdgeVertices(e, v);
      bc.edge_component.push_back(bc.vertex_component[v[0]]);
   }
   return bc;
}

void SetTangents(Mesh &mesh, TractionBoundary &bc)
{
   Array<int> v;
   for (int e : bc.edges)
   {
      mesh.GetEdgeVertices(e, v);
      const real_t *p = mesh.GetVertex(v[0]), *q = mesh.GetVertex(v[1]);
      const real_t l = hypot(q[0]-p[0], q[1]-p[1]);
      const Direction t = {(q[0]-p[0])/l, (q[1]-p[1])/l};
      for (int j = 0; j < 2; j++) { bc.tangents[v[j]].push_back(t); }
   }
}

/** Orthonormal basis of the symmetric h = (xx,xy,yy) with h d = 0 for every
    direction d. One direction leaves a line; two distinct ones leave zero. */
DenseMatrix SymmetricNullSpace(const vector<Direction> &directions)
{
   double G[9] = {0}, w[3], work[16];
   for (const auto &d : directions)
   {
      const double rows[2][3] = {{d[0], d[1], 0}, {0, d[0], d[1]}};
      for (const auto &r : rows)
         for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) { G[i+3*j] += r[i]*r[j]; }
   }
   const int n = 3, lwork = 16; int info;
   dsyev_("V", "U", &n, G, &n, w, work, &lwork, &info);
   MFEM_VERIFY(info == 0, "dsyev failed");
   int k = 0;
   while (k < 3 && w[k] < 1e-10) { k++; }
   DenseMatrix N(3, k);
   for (int j = 0; j < k; j++) for (int i = 0; i < 3; i++) { N(i,j) = G[i+3*j]; }
   return N;
}

Direction Normal(const Direction &t) { return {t[1], -t[0]}; }

/// Free combinations (columns) of the vertex DOFs at a vertex on Gamma_s.
using VertexRule = function<DenseMatrix(const vector<Direction> &, int)>;

// Cartesian tensor values: sigma n = 0 on each incident Gamma_s edge.
DenseMatrix StressRule(const vector<Direction> &t, int ndof)
{
   if (ndof == 0) { return DenseMatrix(0, 0); }
   MFEM_VERIFY(ndof == 3, "expected (xx,xy,yy) vertex DOFs");
   vector<Direction> n;
   for (const auto &d : t) { n.push_back(Normal(d)); }
   return SymmetricNullSpace(n);
}

DenseMatrix ZeroRule(const vector<Direction> &, int ndof) { return DenseMatrix(ndof, 0); }

/** Clamped potential: value and gradient vanish; Hessian H t = 0 for each
    incident Gamma_s tangent, so that Airy(phi) n = 0 at the vertex. The
    deliberately wrong variants leave the Hessian free, or set it to zero. */
enum class HessianBC { Tangent, Free, Zero };
DenseMatrix PotentialRule(HessianBC mode, const vector<Direction> &t, int ndof)
{
   if (ndof == 3) { return DenseMatrix(3, 0); } // HCT: value and gradient
   MFEM_VERIFY(ndof == 6, "expected (v,x,y,xx,xy,yy) vertex DOFs");
   DenseMatrix N;
   if (mode == HessianBC::Tangent) { N = SymmetricNullSpace(t); }
   else if (mode == HessianBC::Free) { N.Diag(1.0, 3); }
   else { N.SetSize(3, 0); }
   DenseMatrix F(6, N.Width());
   F = 0.0;
   for (int j = 0; j < N.Width(); j++) for (int i = 0; i < 3; i++) { F(3+i,j) = N(i,j); }
   return F;
}

/// DOFs of the affine function c0 + c1 x + c2 y, fitted element by element.
Vector AffineDofs(FiniteElementSpace &fes, const real_t c[3])
{
   Vector g(fes.GetVSize());
   g = 0.0;
   Array<int> dofs;
   const IntegrationRule &rule = IntRules.Get(Geometry::TRIANGLE, 12);
   for (int e = 0; e < fes.GetNE(); e++)
   {
      const FiniteElement &fe = *fes.GetFE(e);
      ElementTransformation &T = *fes.GetElementTransformation(e);
      const int nd = fe.GetDof();
      DenseMatrix S(rule.GetNPoints(), nd), StS(nd);
      Vector f(rule.GetNPoints()), shape(nd), rhs(nd), local(nd), x(2);
      for (int q = 0; q < rule.GetNPoints(); q++)
      {
         T.SetIntPoint(&rule.IntPoint(q));
         T.Transform(rule.IntPoint(q), x);
         fe.CalcPhysShape(T, shape);
         S.SetRow(q, shape);
         f(q) = c[0] + c[1]*x(0) + c[2]*x(1);
      }
      // Derivative DOFs scale like powers of h; equilibrate before solving.
      Vector scale(nd);
      for (int j = 0; j < nd; j++)
      {
         real_t norm = 0.0;
         for (int q = 0; q < S.Height(); q++) { norm += S(q,j)*S(q,j); }
         scale(j) = 1.0/sqrt(norm);
      }
      S.RightScaling(scale);
      MultAtB(S, S, StS);
      S.MultTranspose(f, rhs);
      DenseMatrixInverse(StS).Mult(rhs, local);
      for (int j = 0; j < nd; j++) { local(j) *= scale(j); }
      fes.GetElementDofs(e, dofs);
      for (int i = 0; i < nd; i++)
      { g(UnsignIndex(dofs[i])) = (dofs[i] < 0 ? -1.0 : 1.0)*local(i); }
   }
   return g;
}

/** Build P. Interior edge DOFs on Gamma_s are removed. With @a gauge and an
    empty Gamma_s, the value and gradient at vertex 0 are removed. With
    @a cohomology, each Gamma_s component after the first gets three columns
    making the constrained DOFs there those of a free affine function. */
unique_ptr<SparseMatrix> Restriction(FiniteElementSpace &fes,
                                     const TractionBoundary &bc, VertexRule rule,
                                     bool gauge, bool cohomology)
{
   const int n = fes.GetVSize();
   vector<char> fixed(n, 0);
   vector<vector<pair<int,real_t>>> cols;
   vector<DenseMatrix> free(fes.GetMesh()->GetNV());
   Array<int> dofs;
   for (int v = 0; v < fes.GetMesh()->GetNV(); v++)
   {
      if (bc.vertex_component[v] < 0) { continue; }
      fes.GetVertexVDofs(v, dofs);
      free[v] = rule(bc.tangents[v], dofs.Size());
      for (int i = 0; i < dofs.Size(); i++) { fixed[dofs[i]] = 1; }
      for (int j = 0; j < free[v].Width(); j++)
      {
         cols.emplace_back();
         for (int i = 0; i < dofs.Size(); i++)
            if (std::abs(free[v](i,j)) > 1e-14) { cols.back().push_back({dofs[i], free[v](i,j)}); }
      }
   }
   for (int e : bc.edges)
   {
      fes.GetEdgeInteriorVDofs(e, dofs);
      for (int d : dofs) { fixed[UnsignIndex(d)] = 1; }
   }
   if (gauge && bc.edges.empty())
   {
      fes.GetVertexDofs(0, dofs);
      for (int i = 0; i < 3; i++) { fixed[dofs[i]] = 1; }
   }
   for (int i = 0; i < n; i++) if (!fixed[i]) { cols.push_back({{i, 1.0}}); }
   for (int c = 1; cohomology && c < bc.components; c++)
   {
      for (int a = 0; a < 3; a++)
      {
         const real_t coeff[3] = {real_t(a == 0), real_t(a == 1), real_t(a == 2)};
         const Vector g = AffineDofs(fes, coeff);
         cols.emplace_back();
         for (int v = 0; v < fes.GetMesh()->GetNV(); v++)
         {
            if (bc.vertex_component[v] != c) { continue; }
            fes.GetVertexVDofs(v, dofs);
            // Constrained part of the vertex jet: (I - F F^T) g.
            // The vertex jet of an affine function is exact: value, gradient,
            // and a zero Hessian. Only edge DOFs come from the fit.
            Vector gv(dofs.Size()), Ftg(free[v].Width()), proj(dofs.Size());
            const real_t *x = fes.GetMesh()->GetVertex(v);
            gv = 0.0;
            gv(0) = coeff[0] + coeff[1]*x[0] + coeff[2]*x[1];
            gv(1) = coeff[1]; gv(2) = coeff[2];
            proj = 0.0;
            if (free[v].Width() > 0)
            {
               free[v].MultTranspose(gv, Ftg);
               free[v].Mult(Ftg, proj);
            }
            for (int i = 0; i < dofs.Size(); i++)
               if (std::abs(gv(i)-proj(i)) > 1e-13) { cols.back().push_back({dofs[i], gv(i)-proj(i)}); }
         }
         for (size_t k = 0; k < bc.edges.size(); k++)
         {
            if (bc.edge_component[k] != c) { continue; }
            fes.GetEdgeInteriorVDofs(bc.edges[k], dofs);
            for (int d : dofs)
               if (std::abs(g(UnsignIndex(d))) > 1e-13)
               { cols.back().push_back({UnsignIndex(d), g(UnsignIndex(d))}); }
         }
      }
   }
   auto P = unique_ptr<SparseMatrix>(new SparseMatrix(n, int(cols.size())));
   for (size_t j = 0; j < cols.size(); j++)
      for (const auto &entry : cols[j]) { P->Set(entry.first, int(j), entry.second); }
   P->Finalize();
   P->SortColumnIndices();
   return P;
}

SparseMatrix *Restrict(const SparseMatrix &A, const SparseMatrix &P)
{
   SparseMatrix *Ar = RAP(P, A, P);
   Ar->SortColumnIndices();
   return Ar;
}

/// max |(I - Ps Ps^T) M| / max |M|: does the image of M lie in range(Ps)?
real_t Leak(const SparseMatrix &Ps, const SparseMatrix &M)
{
   unique_ptr<SparseMatrix> PsT(Transpose(Ps));
   unique_ptr<SparseMatrix> coords(mfem::Mult(*PsT, M));
   unique_ptr<SparseMatrix> back(mfem::Mult(Ps, *coords));
   unique_ptr<SparseMatrix> diff(Add(1.0, M, -1.0, *back));
   return diff->MaxNorm()/std::max(M.MaxNorm(), real_t(1e-300));
}

/// Vertex patches (as in ex43.hpp) expressed in restricted coordinates.
Table RestrictedPatches(FiniteElementSpace &fes, const SparseMatrix &P)
{
   Mesh &mesh = *fes.GetMesh();
   unique_ptr<Table> v2e(mesh.GetVertexToElementTable());
   Table patches;
   patches.MakeI(mesh.GetNV());
   vector<vector<int>> rows(mesh.GetNV());
   Array<int> dofs, edges, orient, ev;
   vector<int> marker(P.Width(), -1);
   for (int v = 0; v < mesh.GetNV(); v++)
   {
      vector<int> full;
      fes.GetVertexDofs(v, dofs);
      for (int d : dofs) { full.push_back(UnsignIndex(d)); }
      for (int k = 0; k < v2e->RowSize(v); k++)
      {
         const int e = v2e->GetRow(v)[k];
         fes.GetElementInteriorDofs(e, dofs);
         for (int d : dofs) { full.push_back(UnsignIndex(d)); }
         mesh.GetElementEdges(e, edges, orient);
         for (int edge : edges)
         {
            mesh.GetEdgeVertices(edge, ev);
            if (ev.Find(v) < 0) { continue; }
            fes.GetEdgeInteriorDofs(edge, dofs);
            for (int d : dofs) { full.push_back(UnsignIndex(d)); }
         }
      }
      for (int d : full)
      {
         for (int k = P.GetI()[d]; k < P.GetI()[d+1]; k++)
         {
            const int c = P.GetJ()[k];
            if (marker[c] != v) { marker[c] = v; rows[v].push_back(c); }
         }
      }
      patches.AddColumnsInRow(v, int(rows[v].size()));
   }
   patches.MakeJ();
   for (int v = 0; v < mesh.GetNV(); v++)
      for (int c : rows[v]) { patches.AddConnection(v, c); }
   patches.ShiftUpI();
   return patches;
}

class PatchSmoother : public Solver
{
   vector<Array<int>> patches;
   vector<DenseMatrix> inverses;
   real_t scale;
   mutable Vector rhs, sol;
public:
   PatchSmoother(const SparseMatrix &A, const Table &table, real_t s)
      : Solver(A.Height()), scale(s)
   {
      for (int p = 0; p < table.Size(); p++)
      {
         Array<int> dofs;
         table.GetRow(p, dofs);
         if (dofs.Size() == 0) { continue; }
         DenseMatrix local(dofs.Size());
         A.GetSubMatrix(dofs, dofs, local);
         inverses.emplace_back();
         DenseMatrixInverse(local, true).GetInverseMatrix(inverses.back());
         patches.push_back(dofs);
      }
   }
   void Mult(const Vector &x, Vector &y) const override
   {
      y = 0.0;
      for (size_t p = 0; p < patches.size(); p++)
      {
         x.GetSubVector(patches[p], rhs);
         sol.SetSize(rhs.Size());
         inverses[p].Mult(rhs, sol);
         y.AddElementVector(patches[p], sol);
      }
      y *= scale;
   }
   void SetOperator(const Operator &) override { MFEM_ABORT("fixed smoother"); }
};

/// R + sum_k M_k L_k^{-1} M_k^T, in restricted coordinates.
class AdditiveHX : public Solver
{
   const Solver &R;
   vector<const SparseMatrix *> maps;
   vector<const Solver *> inverses;
   mutable Vector rhs, sol, cor;
public:
   AdditiveHX(const Solver &smoother) : Solver(smoother.Height()), R(smoother) { }
   void AddSpace(const SparseMatrix &M, const Solver &inv)
   { maps.push_back(&M); inverses.push_back(&inv); }
   void Mult(const Vector &x, Vector &y) const override
   {
      R.Mult(x, y);
      for (size_t k = 0; k < maps.size(); k++)
      {
         rhs.SetSize(maps[k]->Width()); sol.SetSize(rhs.Size()); cor.SetSize(height);
         maps[k]->MultTranspose(x, rhs);
         inverses[k]->Mult(rhs, sol);
         maps[k]->Mult(sol, cor);
         y += cor;
      }
   }
   void SetOperator(const Operator &) override { MFEM_ABORT("fixed HX"); }
};

struct Spectrum { real_t lo = 0, hi = 0; int steps = 0; };

/// Ritz extremes of BA in the A inner product, full reorthogonalization.
Spectrum Lanczos(const SparseMatrix &A, const Solver &B, int maxit, int seed)
{
   const int n = A.Height();
   maxit = std::min(n, maxit);
   vector<Vector> Q, AQ;
   Vector q(n), aq(n), w(n), aw(n);
   std::mt19937 gen(seed);
   std::normal_distribution<double> normal;
   for (int i = 0; i < n; i++) { q(i) = normal(gen); }
   A.Mult(q, aq);
   real_t norm = sqrt(q*aq);
   q /= norm; aq /= norm;
   vector<double> diag, off;
   Spectrum s, previous;
   for (int k = 0; k < maxit; k++)
   {
      Q.push_back(q); AQ.push_back(aq);
      B.Mult(aq, w);
      diag.push_back(aq*w);
      for (int pass = 0; pass < 2; pass++)
         for (int j = 0; j <= k; j++) { w.Add(-(AQ[j]*w), Q[j]); }
      A.Mult(w, aw);
      const real_t beta = sqrt(std::max(real_t(0), w*aw));
      const int m = k+1;
      if (m % 5 == 0 || m == maxit || beta < 1e-12)
      {
         vector<double> d = diag, e = off;
         e.resize(m);
         vector<double> work(std::max(1, 2*m-2));
         int info;
         dstev_("N", &m, d.data(), e.data(), nullptr, &m, work.data(), &info);
         MFEM_VERIFY(info == 0, "dstev failed");
         s.lo = d.front(); s.hi = d.back(); s.steps = m;
         const bool settled = previous.steps &&
                              std::abs(s.lo-previous.lo) < 1e-6*s.lo &&
                              std::abs(s.hi-previous.hi) < 1e-6*s.hi;
         if (settled || beta < 1e-12) { return s; }
         previous = s;
      }
      off.push_back(beta);
      q = w; q /= beta;
      A.Mult(q, aq);
   }
   return s;
}

DenseMatrix ToDense(const SparseMatrix &A)
{
   DenseMatrix D;
   A.ToDenseMatrix(D);
   return D;
}

/// Number of generalized eigenvalues of (K, M) below tol (M SPD).
int GeneralizedNullity(const SparseMatrix &K, const SparseMatrix &M, real_t tol)
{
   DenseMatrix k = ToDense(K), m = ToDense(M);
   const int n = k.Height(), itype = 1;
   vector<double> w(n), work(1);
   int lwork = -1, info;
   dsygv_(&itype, "N", "U", &n, k.Data(), &n, m.Data(), &n, w.data(), work.data(), &lwork, &info);
   lwork = int(work[0]); work.resize(lwork);
   dsygv_(&itype, "N", "U", &n, k.Data(), &n, m.Data(), &n, w.data(), work.data(), &lwork, &info);
   MFEM_VERIFY(info == 0, "dsygv failed");
   return int(std::count_if(w.begin(), w.end(), [&](double x) { return x < tol; }));
}

/// Numerical rank of a rectangular sparse matrix, via eigenvalues of M^T M.
int Rank(const SparseMatrix &M)
{
   unique_ptr<SparseMatrix> MtM(TransposeMult(M, M));
   DenseMatrix g = ToDense(*MtM);
   const int n = g.Height();
   vector<double> w(n), work(1);
   int lwork = -1, info;
   dsyev_("N", "U", &n, g.Data(), &n, w.data(), work.data(), &lwork, &info);
   lwork = int(work[0]); work.resize(lwork);
   dsyev_("N", "U", &n, g.Data(), &n, w.data(), work.data(), &lwork, &info);
   MFEM_VERIFY(info == 0, "dsyev failed");
   const double top = n ? w.back() : 0.0;
   return int(std::count_if(w.begin(), w.end(), [&](double x) { return x > 1e-11*top; }));
}

DenseMatrix MatrixH1Weight()
{
   DenseMatrix w(3);
   w = 0.0;
   w(0,0) = w(2,2) = 1.0; w(1,1) = 2.0;
   return w;
}

} // namespace

int main(int argc, char *argv[])
{
   string element = "aw", bc_spec = "b", h1bc = "normal", smoother = "patch",
          hessian = "tangent";
   int n0 = 2, rmin = 0, rmax = 4, lanczos = 300, exact_max = 2500;
   bool cohomology = true, header = true, check = false;
   int failures = 0;
   real_t rotation = 0.0, jitter = 0.0, tol = 1e-8;
   OptionsParser args(argc, argv);
   args.AddOption(&element, "-e", "--element", "jm, aw, hz, hzzz, or hzzzr.");
   args.AddOption(&bc_spec, "-bc", "--traction-boundary",
                  "Sides with sigma n = 0: any of b,r,t,l; h = lower half of l; "
                  "a = all; none.");
   args.AddOption(&h1bc, "-h1bc", "--h1-boundary",
                  "Constraint on the P1 tensors at Gamma_s vertices: normal "
                  "(sigma n = 0) or full (sigma = 0).");
   args.AddOption(&cohomology, "-cohom", "--cohomology", "-no-cohom",
                  "--no-cohomology", "Allow affine potentials on extra Gamma_s components.");
   args.AddOption(&hessian, "-hess", "--potential-hessian",
                  "Argyris/Bell vertex Hessian on Gamma_s: tangent (H t = 0, "
                  "correct), free, or zero.");
   args.AddOption(&smoother, "-s", "--smoother", "patch, jacobi, or gs.");
   args.AddOption(&n0, "-n", "--coarse", "Coarse mesh is n x n squares, two triangles each.");
   args.AddOption(&rmin, "-rmin", "--min-refinement", "First refinement level.");
   args.AddOption(&rmax, "-rmax", "--max-refinement", "Last refinement level.");
   args.AddOption(&rotation, "-rot", "--rotation", "Rotate the square (degrees).");
   args.AddOption(&jitter, "-jitter", "--jitter",
                  "Random interior vertex displacement, as a fraction of h.");
   args.AddOption(&lanczos, "-lanczos", "--lanczos", "Maximum Lanczos steps (0 = off).");
   args.AddOption(&exact_max, "-exact-max", "--exact-max",
                  "Largest restricted stress size for dense exactness checks.");
   args.AddOption(&tol, "-tol", "--tolerance", "PCG relative tolerance.");
   args.AddOption(&header, "-header", "--header", "-no-header", "--no-header", "Print CSV header.");
   args.AddOption(&check, "-check", "--check", "-no-check", "--no-check",
                  "Exit with status 3 unless both maps land in the restricted "
                  "space, the checked complex is exact, and PCG converges.");
   args.ParseCheck(cerr);

   string stress_name, potential_name;
   if (element == "jm") { stress_name = "JM_2D_P1"; potential_name = "HCT_2D_P3"; }
   else if (element == "aw") { stress_name = "AW_2D_P3"; potential_name = "Argyris_2D_P5"; }
   else if (element == "hz") { stress_name = "HZ_2D_P3"; potential_name = "Argyris_2D_P5"; }
   else if (element == "hzzz") { stress_name = "HZZZ_2D_P3"; potential_name = "Bell_2D_P5"; }
   else if (element == "hzzzr") { stress_name = "HZZZr_2D_P3"; potential_name = "Bell_2D_P5"; }
   else { MFEM_ABORT("unknown element " << element); }
   MFEM_VERIFY(h1bc == "normal" || h1bc == "full", "-h1bc is normal or full");
   if (bc_spec == "none") { bc_spec = ""; }
   MFEM_VERIFY(hessian == "tangent" || hessian == "free" || hessian == "zero",
               "-hess is tangent, free, or zero");
   const HessianBC hessian_bc = hessian == "tangent" ? HessianBC::Tangent :
                                hessian == "free" ? HessianBC::Free : HessianBC::Zero;
   using namespace std::placeholders;
   const VertexRule potential_rule = std::bind(PotentialRule, hessian_bc, _1, _2);

   if (header)
   {
      cout << "element,bc,h1bc,hess,cohom,smoother,rot,jitter,level,h,components,"
           "n_stress,n_h1,n_potential,leak_pi,leak_airy,ker_div,rank_airy,"
           "pcg_its,lambda_min,lambda_max,kappa,lanczos_steps\n";
   }
   for (int level = rmin; level <= rmax; level++)
   {
      const int N = n0 << level;
      Mesh mesh = Mesh::MakeCartesian2D(N, N, Element::TRIANGLE, true);
      TractionBoundary bc = MarkBoundary(mesh, bc_spec);
      std::mt19937 gen(level+17);
      std::uniform_real_distribution<double> unif(-1.0, 1.0);
      const real_t c = cos(rotation*M_PI/180), s = sin(rotation*M_PI/180);
      for (int v = 0; v < mesh.GetNV(); v++)
      {
         real_t *x = mesh.GetVertex(v);
         const bool interior = x[0] > 1e-10 && x[0] < 1-1e-10 &&
                               x[1] > 1e-10 && x[1] < 1-1e-10;
         if (interior && jitter > 0)
         {
            x[0] += jitter*unif(gen)/N;
            x[1] += jitter*unif(gen)/N;
         }
         const real_t x0 = x[0], x1 = x[1];
         x[0] = c*x0 - s*x1; x[1] = s*x0 + c*x1;
      }
      SetTangents(mesh, bc);

      unique_ptr<FiniteElementCollection> stress_fec(FiniteElementCollection::New(stress_name.c_str()));
      unique_ptr<FiniteElementCollection> potential_fec(FiniteElementCollection::New(potential_name.c_str()));
      H1_FECollection h1_fec(1, 2);
      FiniteElementSpace stress(&mesh, stress_fec.get());
      FiniteElementSpace potential(&mesh, potential_fec.get());
      FiniteElementSpace h1(&mesh, &h1_fec, 3, Ordering::byVDIM);

      auto Ps = Restriction(stress, bc, StressRule, false, false);
      auto Ph = Restriction(h1, bc, h1bc == "full" ? VertexRule(ZeroRule) : VertexRule(StressRule),
                            false, false);
      auto Pp = Restriction(potential, bc, potential_rule, true, cohomology);

      BilinearForm mass(&stress), divdiv(&stress);
      mass.AddDomainIntegrator(new MatrixFEMassIntegrator);
      divdiv.AddDomainIntegrator(new MatrixDivDivIntegrator);
      mass.Assemble(); mass.Finalize();
      divdiv.Assemble(); divdiv.Finalize();
      unique_ptr<SparseMatrix> A(Add(mass.SpMat(), divdiv.SpMat()));
      unique_ptr<SparseMatrix> Ar(Restrict(*A, *Ps));

      MatrixConstantCoefficient weight(MatrixH1Weight());
      BilinearForm h1_form(&h1), hess(&potential);
      h1_form.AddDomainIntegrator(new VectorMassIntegrator(weight));
      h1_form.AddDomainIntegrator(new VectorDiffusionIntegrator(weight));
      hess.AddDomainIntegrator(new HessianIntegrator);
      h1_form.Assemble(); h1_form.Finalize();
      hess.Assemble(); hess.Finalize();
      unique_ptr<SparseMatrix> Lh(Restrict(h1_form.SpMat(), *Ph));
      unique_ptr<SparseMatrix> Lp(Restrict(hess.SpMat(), *Pp));

      DiscreteLinearOperator pi(&h1, &stress), airy(&potential, &stress);
      pi.AddDomainInterpolator(new IdentityInterpolator);
      airy.AddDomainInterpolator(new AiryInterpolator);
      pi.Assemble(); pi.Finalize();
      airy.Assemble(); airy.Finalize();
      unique_ptr<SparseMatrix> PiP(mfem::Mult(pi.SpMat(), *Ph));
      unique_ptr<SparseMatrix> AiryP(mfem::Mult(airy.SpMat(), *Pp));
      const real_t leak_pi = Leak(*Ps, *PiP), leak_airy = Leak(*Ps, *AiryP);
      unique_ptr<SparseMatrix> PsT(Transpose(*Ps));
      unique_ptr<SparseMatrix> Pir(mfem::Mult(*PsT, *PiP));
      unique_ptr<SparseMatrix> Airyr(mfem::Mult(*PsT, *AiryP));

      int ker_div = -1, rank_airy = -1;
      if (Ar->Height() <= exact_max)
      {
         unique_ptr<SparseMatrix> Mr(Restrict(mass.SpMat(), *Ps));
         unique_ptr<SparseMatrix> Dr(Restrict(divdiv.SpMat(), *Ps));
         ker_div = GeneralizedNullity(*Dr, *Mr, 1e-8);
         rank_airy = Rank(*Airyr);
      }

      unique_ptr<Solver> R;
      Table patches;
      if (smoother == "patch")
      {
         patches = RestrictedPatches(stress, *Ps);
         R.reset(new PatchSmoother(*Ar, patches, 0.33));
      }
      else if (smoother == "jacobi") { R.reset(new DSmoother(*Ar, DSmoother::JACOBI)); }
      // Forward/backward sweeps keep the preconditioner symmetric for CG.
      else if (smoother == "gs") { R.reset(new GSSmoother(*Ar, GSSmoother::SYMMETRIC)); }
      else { MFEM_ABORT("unknown smoother " << smoother); }
      UMFPackSolver h1_inv(*Lh), pot_inv(*Lp);
      AdditiveHX hx(*R);
      hx.AddSpace(*Pir, h1_inv);
      hx.AddSpace(*Airyr, pot_inv);

      Vector b(Ar->Height()), x(Ar->Height());
      b.Randomize(1);
      x = 0.0;
      CGSolver cg;
      cg.SetOperator(*Ar); cg.SetPreconditioner(hx);
      cg.SetRelTol(tol); cg.SetAbsTol(0.0); cg.SetMaxIter(2000); cg.SetPrintLevel(-1);
      cg.Mult(b, x);

      Spectrum sp;
      if (lanczos > 0) { sp = Lanczos(*Ar, hx, lanczos, 7); }
      cout << element << ',' << (bc_spec.empty() ? "none" : bc_spec) << ',' << h1bc << ','
           << hessian << ',' << cohomology << ',' << smoother << ',' << rotation << ',' << jitter << ','
           << level << ',' << 1.0/N << ',' << bc.components << ','
           << Ar->Height() << ',' << Lh->Height() << ',' << Lp->Height() << ','
           << leak_pi << ',' << leak_airy << ',' << ker_div << ',' << rank_airy << ','
           << (cg.GetConverged() ? cg.GetNumIterations() : -cg.GetNumIterations()) << ','
           << sp.lo << ',' << sp.hi << ',' << (sp.lo > 0 ? sp.hi/sp.lo : 0.0) << ','
           << sp.steps << endl;
      if (check && (leak_pi > 1e-8 || leak_airy > 1e-8 || ker_div != rank_airy ||
                    !cg.GetConverged()))
      {
         cerr << "check failed at level " << level << '\n';
         failures++;
      }
   }
   return failures ? 3 : 0;
}

// Weighted residual-correction compositions for the HX spectral experiments.
#ifndef ADDISONS_EXPERIMENTS_HX_COMPOSITION_HPP
#define ADDISONS_EXPERIMENTS_HX_COMPOSITION_HPP
#include "common.hpp"
#include <string>
namespace mfem
{
namespace jm_hx
{
struct SpectrumConfig { std::string name, order; double s,h,c; };
class SpectrumComposition : public Solver
{
   const SparseMatrix &A;
   const Solver &S;
   const AuxiliaryCorrection &aux;
   SpectrumConfig cfg;
   mutable Vector r,z,az,other;
public:
   SpectrumComposition(const SparseMatrix &a, const Solver &s,
                       const AuxiliaryCorrection &q, SpectrumConfig f)
      : Solver(a.Height()), A(a), S(s), aux(q), cfg(f), r(height),z(height),
        az(height),other(height) {}
   void Term(char t, const Vector &x, Vector &y) const
   {
      if (t=='S') { S.Mult(x,y); y*=cfg.s; }
      if (t=='H') { aux.H1Mult(x,y); y*=cfg.h; }
      if (t=='C') { aux.AiryMult(x,y); y*=cfg.c; }
   }
   void Mult(const Vector &x, Vector &y) const override
   {
      y=0.;
      if (cfg.order=="add")
      { for (char t : std::string("SHC")) { Term(t,x,z); y+=z; } }
      else if (cfg.order=="H(C+S)H")
      {
         Term('H',x,y); A.Mult(y,az); r=x; r-=az;
         // Both summands see the SAME residual, with independent weights.
         Term('C',r,z); Term('S',r,other); z+=other;
         y+=z; A.Mult(z,az); r-=az;
         Term('H',r,z); y+=z;
      }
      else
      {
         r=x;
         for (char t : cfg.order)
         { Term(t,r,z); y+=z; A.Mult(z,az); r-=az; }
      }
   }
   void MultTranspose(const Vector &x, Vector &y) const override { Mult(x,y); }
   void SetOperator(const Operator&) override { MFEM_ABORT("fixed composition"); }
};
}
}
#endif

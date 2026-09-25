// Energy-inner-product Lanczos for BA: restarted ARPACK or full reorthogonalization.
// See reports/hx_spectrum.md for the pencil, projection identity, and reproduction.
#include "hx_composition.hpp"
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <limits>
#include <string>
#include <fstream>
#include <sstream>
#include <set>
extern "C" void dstev_(const char*, const int*, double*, double*, double*,
                       const int*, double*, int*);
extern "C" void dgeev_(const char*, const char*, const int*, double*,
                       const int*,
                       double*, double*, double*, const int*, double*, const int*, double*, const int*,
                       int*);
// ARPACK-ng LP64 C interface; generalized symmetric mode 2.
extern "C" void dsaupd_c(int*,const char*,int,const char*,int,double,double*,
                         int,double*,int,int*,int*,double*,double*,int,int*);
extern "C" void dseupd_c(int,const char*,const int*,double*,double*,int,double,
                         const char*,int,const char*,int,double,double*,int,double*,int,int*,int*,
                         double*,double*,int,int*);
using namespace mfem;
using namespace mfem::jm_hx;
struct Estimate { double lo,hi,rl,rh; int steps; bool ok; bool bounded; double upper_gap; };
Estimate Lanczos(const SparseMatrix &A, const Operator &B, int seed, int maxit,
                 double tol, double top_bound=0., double bound_tol=1e-7)
{
   int n=A.Height(); maxit=std::min(n,maxit);
   std::vector<Vector> Q,AQ;
   Vector q(n),aq(n),w(n),aw(n),u(n),au(n),bu(n);
   std::mt19937 gen(seed); std::normal_distribution<double> normal;
   for (int i=0; i<n; i++) { q(i)=normal(gen); }
   A.Mult(q,aq); double norm=std::sqrt(q*aq); q/=norm; aq/=norm;
   std::vector<double> diag,off;
   Estimate result{};
   for (int k=0; k<maxit; k++)
   {
      Q.emplace_back(q); AQ.emplace_back(aq);
      B.Mult(aq,w); double alpha=aq*w;
      diag.push_back(alpha);
      // Two passes avoid ghost Ritz values. Store A*q to avoid one SpMV per dot.
      for (int pass=0; pass<2; pass++)
         for (int j=0; j<=k; j++) { double h=AQ[j]*w; w.Add(-h,Q[j]); }
      A.Mult(w,aw); double beta=std::sqrt(std::max(0.,w*aw));
      int m=k+1;
      if (m%10==0 || m==maxit || beta<1e-13)
      {
         auto d=diag; auto e=off; e.resize(m);
         std::vector<double> v(m*m),work(std::max(1,2*m-2)); int info;
         const char job='V';
         dstev_(&job,&m,d.data(),e.data(),v.data(),&m,work.data(),&info);
         MFEM_VERIFY(info==0,"tridiagonal eigensolve failed");
         result.lo=d.front(); result.hi=d.back(); result.steps=m;
         double residuals[2];
         for (int t=0; t<2; t++)
         {
            int col=t ? m-1:0; u=0.;
            for (int j=0; j<m; j++) { u.Add(v[j+col*m],Q[j]); }
            A.Mult(u,au); B.Mult(au,bu);
            // The true Rayleigh quotient is a variational lower bound on lmax.
            if (t==1 && top_bound>0.) { result.hi=(au*bu)/(u*au); }
            bu.Add(-d[col],u);
            A.Mult(bu,aw);
            residuals[t]=std::sqrt(std::max(0.,bu*aw)/(u*au))/std::abs(d[col]);
         }
         result.rl=residuals[0]; result.rh=residuals[1];
         MFEM_VERIFY(top_bound<=0. ||
                     result.hi<=top_bound*(1.+1e-8),"smoother-center upper bound violated");
         result.upper_gap=top_bound>0. ? std::max(0.,top_bound-result.hi) : 0.;
         result.bounded=top_bound>0. && result.upper_gap<=bound_tol*top_bound;
         result.ok=result.rl<tol && (result.rh<tol || result.bounded);
         if (result.ok || beta<1e-13) { return result; }
      }
      off.push_back(beta); q=w; q/=beta; A.Mult(q,aq);
   }
   return result;
}
Estimate RestartedLanczos(const SparseMatrix &A, const Operator &B, int seed,
                          int maxit, double tol, int ncv, bool projector_center)
{
   int n=A.Height(), nev=projector_center?1:2, ido=0, info=1;
   const char *which=projector_center ? "SA" : "BE";
   ncv=std::min(n,ncv); int lwork=ncv*(ncv+8);
   std::vector<double> resid(n),v(n*ncv),wd(3*n),wl(lwork),d(nev),z(n*nev);
   std::vector<int> select(ncv); int iparam[11]= {},ipntr[11]= {};
   iparam[0]=1; iparam[2]=maxit; iparam[6]=2;
   std::mt19937 gen(seed); std::normal_distribution<double> normal;
   for (double &x:resid) { x=normal(gen); }
   Vector ax(n),tmp(n);
   // Generalized pencil (A B A, A): OP=B A, metric=A. No A inverse is needed.
   do
   {
      dsaupd_c(&ido,"G",n,which,nev,tol,resid.data(),ncv,v.data(),n,
               iparam,ipntr,wd.data(),wl.data(),lwork,&info);
      if (ido==-1 || ido==1)
      {
         Vector x(wd.data()+ipntr[0]-1,n),y(wd.data()+ipntr[1]-1,n);
         A.Mult(x,ax); B.Mult(ax,y);
         if (ido==1) { A.Mult(y,tmp); x=tmp; }
      }
      else if (ido==2)
      {
         Vector x(wd.data()+ipntr[0]-1,n),y(wd.data()+ipntr[1]-1,n);
         A.Mult(x,y);
      }
      else { MFEM_VERIFY(ido==99,"unexpected ARPACK request"); }
   }
   while (ido!=99);
   if (info!=0)
   {
      std::cerr << "ARPACK did not converge: " << info << ", converged " << iparam[4]
                << std::endl;
      double nan=std::numeric_limits<double>::quiet_NaN();
      return {nan,nan,nan,nan,iparam[8],false};
   }
   dseupd_c(1,"A",select.data(),d.data(),z.data(),n,0.,"G",n,which,nev,tol,
            resid.data(),ncv,v.data(),n,iparam,ipntr,wd.data(),wl.data(),lwork,&info);
   MFEM_VERIFY(info==0,"ARPACK extraction failed");
   Estimate e{d[0],projector_center?1.:d[1],0.,0.,iparam[8],false};
   for (int t=0; t<(projector_center?1:2); t++)
   {
      Vector u(z.data()+t*n,n); A.Mult(u,ax); double norm=u*ax;
      B.Mult(ax,tmp); tmp.Add(-d[t],u); A.Mult(tmp,ax);
      double r=std::sqrt(std::max(0.,tmp*ax)/norm)/std::abs(d[t]);
      (t?e.rh:e.rl)=r;
   }
   // For a palindrome with exact Airy projection P in the center,
   // I-BA = F^*(I-P)F in the A metric. It is PSD and rank deficient;
   // hence lambda_max(BA)=1, even if lambda_min is negative.
   // res_max=0 is an identity placeholder, not a computed Ritz residual.
   if (projector_center) { e.hi=1.; e.rh=0.; }
   e.ok=std::max(e.rl,e.rh)<std::max(10*tol,1e-7);
   return e;
}
int main(int argc,char**argv)
{
   const char *meshfile="data/inline-tri.mesh", *selection="all";
   const char *patch="all", *h1="all", *order="add", *skip_file="";
   double ws=1.,wh=1.,wc=1.;
   int ref=0,levels=1,maxit=2000,seed=17,ncv=80; double tol=1e-9;
   bool restarted=true, bound_center=true;
   int full_cap=600; double bound_tol=1e-7;
   OptionsParser args(argc,argv);
   args.AddOption(&meshfile,"-m","--mesh","Square triangle mesh.");
   args.AddOption(&ref,"-r","--refinements","Initial refinement.");
   args.AddOption(&levels,"-levels","--levels","Refinement levels.");
   args.AddOption(&maxit,"-max-it","--max-iterations",
                  "Maximum steps (full) or restart cycles (restarted).");
   args.AddOption(&restarted,"-restart","--restart","-full","--full",
                  "Use implicitly restarted Lanczos.");
   args.AddOption(&bound_center,"-bound-center","--bound-center",
                  "-no-bound-center","--no-bound-center",
                  "Use full Lanczos and the 3*wS upper bound for weight-one Airy, smoother-centered sweeps.");
   args.AddOption(&full_cap,"-full-cap","--full-cap",
                  "Maximum Lanczos vectors in automatic bounded-center mode.");
   args.AddOption(&bound_tol,"-bound-tol","--bound-tolerance",
                  "Relative width allowed for the upper-eigenvalue bracket.");
   args.AddOption(&ncv,"-ncv","--krylov-size","Restarted Krylov dimension.");
   args.AddOption(&seed,"-seed","--seed","Random seed.");
   args.AddOption(&tol,"-tol","--tolerance","Relative A-norm Ritz residual.");
   args.AddOption(&selection,"-select","--selection",
                  "all, core, fine, additive_limit, chosen, tune, middle, middle_tune, middle_refine, middle_chosen, middle_best, custom, or exact config name.");
   args.AddOption(&patch,"-patch","--patch","all, small, large.");
   args.AddOption(&h1,"-h1","--h1","all, split, unsplit.");
   args.AddOption(&order,"-order","--order",
                  "Custom composition: add, H(C+S)H, or a palindrome of S,H,C.");
   args.AddOption(&ws,"-ws","--weight-s","Custom patch weight.");
   args.AddOption(&wh,"-wh","--weight-h","Custom H1 weight.");
   args.AddOption(&wc,"-wc","--weight-c","Custom Airy weight.");
   args.AddOption(&skip_file,"-skip","--skip-completed",
                  "Skip configurations already recorded in this CSV (same seed and method). Output only new rows.");
   args.ParseCheck(std::cerr);
   MFEM_VERIFY(ref>=0 && levels>0 && maxit>0 && tol>0 &&
               ncv>2,"invalid numerical options");
   MFEM_VERIFY(std::string(patch)=="all" || std::string(patch)=="small" ||
               std::string(patch)=="large","invalid patch");
   MFEM_VERIFY(std::string(h1)=="all" || std::string(h1)=="split" ||
               std::string(h1)=="unsplit","invalid H1 space");
   std::string ord(order),reverse(ord.rbegin(),ord.rend());
   MFEM_VERIFY(ord=="add" || ord=="H(C+S)H" || (!ord.empty() && ord==reverse &&
                                                ord.find_first_not_of("SHC")==std::string::npos),
               "multiplicative order must be a palindrome of S,H,C");
   MFEM_VERIFY(ws>0 && wh>0 && wc>0,"weights must be positive");
   auto key=[](int r, const std::string &p, const std::string &h,
               const std::string &name, const std::string &composition,
               double a, double b, double c, int seed,
               const std::string &method)
   {
      std::ostringstream out;
      out << std::setprecision(12) << r << ',' << p << ',' << h << ',' << name << ','
          << composition << ','
          << a << ',' << b << ',' << c << ',' << seed << ',' << method;
      return out.str();
   };
   std::set<std::string> completed;
   if (std::string(skip_file)!="")
   {
      std::ifstream input(skip_file); MFEM_VERIFY(input,"cannot open skip CSV");
      std::string line; std::getline(input,line);
      while (std::getline(input,line))
      {
         std::istringstream in(line); std::vector<std::string> f; std::string field;
         while (std::getline(in,field,',')) { f.push_back(field); }
         MFEM_VERIFY(f.size()>=23,"invalid skip CSV");
         if (f[21]=="1") completed.insert(key(std::stoi(f[0]),f[4],f[5],f[6],f[7],
                                                 std::stod(f[8]),std::stod(f[9]),std::stod(f[10]),std::stoi(f[11]),f[12]));
      }
   }
   Mesh mesh(meshfile);
   for (int i=0; i<ref; i++) { mesh.UniformRefinement(); }
   std::cout<<std::setprecision(12)
            <<"refinement,elements,dofs,h1_dofs,patch,h1,config,composition,wS,wH,wC,seed,method,ncv,tolerance,steps,lambda_min,lambda_max,kappa,res_min,res_max,converged,seconds,max_method,dense_min,dense_max,dense_imag,max_upper_bound,max_upper_gap\n";
   bool ok=true;
   int selected=0;
   for (int lev=0; lev<levels; lev++)
   {
      if (lev) { mesh.UniformRefinement(); }
      JohnsonMercierFECollection mfec,vfec(JMBasis::SplitVertex);
      FiniteElementSpace m(&mesh,&mfec),v(&mesh,&vfec);
      BilinearForm am(&m),av(&v);
      for (auto a : {&am,&av})
      { a->AddDomainIntegrator(new MatrixFEMassIntegrator); a->AddDomainIntegrator(new MatrixDivDivIntegrator); a->Assemble(); a->Finalize(); }
      auto map=BasisMatrix(m,v);
      for (bool split : {false,true})
      {
         if (std::string(h1)!="all" && std::string(h1)!=(split?"split":"unsplit")) { continue; }
         AuxiliaryCorrection aux(m,v,*map,split);
         for (bool small : {false,true})
         {
            if (std::string(patch)!="all" && std::string(patch)!=(small?"small":"large")) { continue; }
            PatchSolver local(av.SpMat(),v,small); MappedSolver S(*map,local);
            double safe=small ? .5 : .33;
            std::vector<SpectrumConfig> configs= {{"base","add",1,1,1},
               {"s01","add",.1,1,1},{"s033","add",.33,1,1},
               {"s3","add",3,1,1},{"s10","add",10,1,1},
               {"h033","add",1,.33,1},{"h3","add",1,3,1},
               {"c033","add",1,1,.33},{"c3","add",1,1,3},
               {"balanced","add",safe,1,1},
               {"mult_SHCHS","SHCHS",safe,1,1},
               {"mult_HSCSH","HSCSH",safe,1,1},
               {"mult_SCHCS","SCHCS",safe,1,1},
               {"mult_raw","SHCHS",1,1,1},
               {"chosen_add","add",small?(split?.75:1.):.75,small?1.:2.,1.},
               {"chosen_mult","SHCHS",small?.6:.4,1.,1.}
            };
            if (std::string(selection)=="custom") configs= {{"custom",order,ws,wh,wc}};
            if (std::string(selection)=="tune")
            {
               configs.clear();
               for (double s : {.25,.33,.5,.75,1.,1.5,2.})
                  for (double h : {.5,1.,2.,3.})
                     configs.push_back({"tune_add","add",s,h,1.});
               for (double s : {.2,.33,.4,.5,.6})
                  for (double h : {.5,1.})
                     configs.push_back({"tune_mult","SHCHS",s,h,1.});
            }
            if (std::string(selection)=="middle")
            {
               double matched=small?.6:.4;
               configs= {{"middle_HCSCH","HCSCH",1,1,1},
                  {"middle_CHSHC","CHSHC",1,1,1},
                  {"middle_HCSCH_matched","HCSCH",matched,1,1},
                  {"middle_CHSHC_matched","CHSHC",matched,1,1},
                  {"hybrid_unit","H(C+S)H",1,1,1},
                  {"hybrid_matched","H(C+S)H",matched,1,1},
                  {"hybrid_c05","H(C+S)H",1,1,.5},
                  {"hybrid_s2","H(C+S)H",2,1,1},
                  {"hybrid_h05","H(C+S)H",1,.5,1},
                  {"hybrid_s2h05","H(C+S)H",2,.5,1},
                  {"reference_SHCHS","SHCHS",matched,1,1}
               };
            }
            if (std::string(selection)=="middle_tune")
            {
               configs.clear();
               for (const std::string ord : {"HCSCH","CHSHC"})
                  for (double s : {.5,1.,2.,4.})
                     for (double h : {.5,1.})
                        configs.push_back({"tune_middle",ord,s,h,1.});
               for (double s : {.5,1.,2.})
                  for (double c : {.5,1.,2.})
                     for (double h : {.5,1.})
                        configs.push_back({"tune_hybrid","H(C+S)H",s,h,c});
            }
            if (std::string(selection)=="middle_refine")
            {
               configs.clear();
               for (double s : {.25,1./3,.4,.5,.6,.75,1.})
               {
                  for (const std::string ord : {"HCSCH","CHSHC"})
                     configs.push_back({"refine_middle",ord,s,1.,1.});
                  for (double c : {.25,.4,.5,.6,.75,1.})
                     configs.push_back({"refine_hybrid","H(C+S)H",s,1.,c});
               }
            }
            if (std::string(selection)=="middle_chosen" ||
                std::string(selection)=="middle_best")
            {
               configs.clear();
               const double s=small ? (split?(std::string(selection)=="middle_best"?.6:.5):.4)
                              : 1./3;
               configs.push_back({"chosen_middle","HCSCH",s,1.,1.});
               configs.push_back({"chosen_middle","CHSHC",s,1.,1.});
               configs.push_back({"chosen_hybrid","H(C+S)H",s,1.,.5});
            }
            for (auto cfg:configs)
            {
               std::string sel(selection);
               // Indefiniteness is already established by the coarse dense check.
               if (sel=="all" && restarted && cfg.name=="mult_raw" && ref+lev>0) { continue; }
               if (sel=="all" && cfg.name.find("chosen_")==0) { continue; }
               if (sel=="chosen" && cfg.name.find("chosen_")!=0) { continue; }
               if (sel=="fine" && cfg.name!="base" && cfg.name!="s033" &&
                   cfg.name!="mult_SHCHS") { continue; }
               if (sel=="additive_limit" && cfg.name!="base" && cfg.name!="s033") { continue; }
               if (sel=="core" && cfg.name!="base" && cfg.name!="s033" && cfg.name!="s3" &&
                   cfg.name!="h3" && cfg.name!="mult_SHCHS" && cfg.name!="mult_HSCSH") { continue; }
               if (sel!="all" && sel!="core" && sel!="fine" && sel!="additive_limit" &&
                   sel!="chosen" && sel!="tune" && sel!="middle" && sel!="middle_tune" &&
                   sel!="middle_refine" && sel!="middle_chosen" && sel!="middle_best" &&
                   sel!=cfg.name) { continue; }
               selected++;
               double top_bound=0.;
               if (bound_center && (cfg.order=="HCSCH" || cfg.order=="CHSHC") && cfg.c==1. &&
                   cfg.h<=1.)
               {
                  // Patch support overlap <=3 gives lambda_max(SA)<=3.
                  // H1 inclusion gives lambda_max(HA)<=2; I-hHA is a contraction.
                  // With weight-one Airy, I-CA is an orthogonal projection.
                  top_bound=std::max(1.,3*cfg.s);
               }
               const std::string method=top_bound>1.?"full_bound":(restarted
                                                                   ?"restarted":"full");
               if (completed.count(key(ref+lev,small?"small":"large",split?"split":"unsplit",
                                       cfg.name,cfg.order,cfg.s,cfg.h,cfg.c,seed,method))) { continue; }
               std::cerr << "Estimating r=" << ref+lev << ' ' << (small?"small":"large")
                         << ' ' << (split?"split":"unsplit") << ' ' << cfg.name << std::endl;
               SpectrumComposition B(am.SpMat(),S,aux,cfg);
               auto start=std::chrono::steady_clock::now();
               bool projector_center=cfg.order!="add" && cfg.order!="H(C+S)H" &&
                                     cfg.order[cfg.order.size()/2]=='C' &&
                                     cfg.c==1.;
               // A C factor makes the error rank deficient. With top_bound=1,
               // the upper endpoint is exactly one even when S is central.
               projector_center=projector_center || top_bound==1.;
               auto e=top_bound>1. ? Lanczos(am.SpMat(),B,seed,full_cap,tol,top_bound,
                                             bound_tol) :
                      (restarted ? RestartedLanczos(am.SpMat(),B,seed,maxit,tol,ncv,
                                                    projector_center) : Lanczos(am.SpMat(),B,seed,maxit,tol));
               double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()
                                                        -start).count();
               double dense_lo=0.,dense_hi=0.,dense_imag=0.;
               if (m.GetVSize()<=320 && e.ok)
               {
                  int n=m.GetVSize(),one=1,info,lwork=8*n;
                  DenseMatrix ba(n); Vector unit(n),aunit(n),bunit(n);
                  for (int j=0; j<n; j++)
                  {
                     unit=0.; unit(j)=1.; am.SpMat().Mult(unit,aunit); B.Mult(aunit,bunit);
                     for (int i=0; i<n; i++) { ba(i,j)=bunit(i); }
                  }
                  std::vector<double> wr(n),wi(n),work(lwork); double dummy;
                  const char no='N';
                  dgeev_(&no,&no,&n,ba.Data(),&n,wr.data(),wi.data(),&dummy,&one,&dummy,&one,
                         work.data(),&lwork,&info);
                  MFEM_VERIFY(info==0,"dense eigenvalue validation failed");
                  dense_lo=*std::min_element(wr.begin(),wr.end());
                  dense_hi=*std::max_element(wr.begin(),wr.end());
                  for (double x:wi) { dense_imag=std::max(dense_imag,std::abs(x)); }
                  MFEM_VERIFY(std::abs(e.lo-dense_lo)<1e-6*std::abs(dense_lo) &&
                              std::abs(e.hi-dense_hi)<1e-6*std::abs(dense_hi) &&
                              dense_imag<1e-8,"Lanczos disagrees with dense spectrum");
               }
               ok=ok&&e.ok;
               std::cout<<ref+lev<<','<<mesh.GetNE()<<','<<m.GetVSize()<<','<<aux.H1Size()
                        <<','<<(small?"small":"large")<<','<<(split?"split":"unsplit")
                        <<','<<cfg.name<<','<<cfg.order<<','<<cfg.s<<','<<cfg.h<<','<<cfg.c<<','<<seed<<','<<
                        method
                        <<','<<ncv<<','<<tol<<','<<e.steps<<','<<e.lo<<','<<e.hi<<','<<
                        (e.lo>0?e.hi/e.lo:std::numeric_limits<double>::quiet_NaN())
                        <<','<<e.rl<<','<<e.rh<<','<<e.ok<<','<<sec<<','<<(restarted&&
                                                                           projector_center?"projector_identity":(e.bounded?"lanczos_bracket":"lanczos"))
                        <<','<<dense_lo<<','<<dense_hi<<','<<dense_imag<<','<<top_bound<<','<<e.upper_gap<<std::endl;
            }
         }
      }
   }
   MFEM_VERIFY(selected>0,"unknown or empty configuration selection");
   return ok?0:3;
}

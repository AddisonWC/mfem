#!/bin/sh
# Reproduce results/hx_bc.csv. Run from the experiment build directory:
#   sh ../johnson_mercier/hx_bc_sweep.sh [jobs] > hx_bc.csv
# Level 6 runs need about 3.5 GB each and are run one at a time.
jobs=${1:-6}
export OPENBLAS_NUM_THREADS=1
./hx_bc -rmax -1 2>/dev/null # header only
{
for e in jm aw hz hzzz hzzzr; do
   for bc in none b h bl blt lr a; do echo "-e $e -bc $bc"; done
done
for e in jm aw; do
   for bc in b h bl a; do echo "-e $e -bc $bc -h1bc full"; done
   for bc in none bl a; do
      for s in jacobi gs; do echo "-e $e -bc $bc -s $s"; done
   done
   for bc in bl a; do echo "-e $e -bc $bc -jitter 0.2"; done
done
for e in jm aw hz; do echo "-e $e -bc lr -no-cohom"; done
for e in jm aw hzzz; do
   for bc in h bl a; do echo "-e $e -bc $bc -rot 30"; done
done
for e in aw hzzz; do
   for bc in h bl a; do
      for hs in free zero; do echo "-e $e -bc $bc -hess $hs"; done
   done
done
} | xargs -P "$jobs" -I{} sh -c './hx_bc -rmax 5 {} -no-header 2>/dev/null'
for o in "-bc none" "-bc a" "-bc bl" "-bc blt" "-bc lr" "-bc lr -no-cohom" "-bc bl -h1bc full"; do
   ./hx_bc -e aw $o -rmin 6 -rmax 6 -lanczos 200 -no-header 2>/dev/null
done

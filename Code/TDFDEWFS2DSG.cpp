#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define PI 3.141592654
using namespace std;

const int nOrder = 3;
const int nTimePreSnap = 100;

typedef struct {
  int nx, nz;
  int Nx, Nz;
  int sx, sz;
  int npx, npz;
  float dx, dz;
} dim;
typedef struct {
  float *vp, *vs, *rho;
} media;
typedef struct {
  float *vxx, *vxz, *vzx, *vzz,
        *txxx, *txxz, *tzzx, *tzzz,
        *txzx, *txzz;
  float *vxt, *vzt, *txxt, *tzzt, *txzt;
} wave;
typedef struct {
  int nt;
  float dt;
  float ampl, f0, t0;
  float d0x, d0z;
  float C[nOrder];
} coeff;

void wave_exp(dim D, char *filename, float *P) {
//
  FILE *fp = fopen(filename, "wb");
    fwrite(&D.nx, sizeof(float), 1, fp);
    fwrite(&D.nz, sizeof(float), 1, fp);
    for(int i = 0; i < D.nz; i++) {
      fwrite(&P[(i + D.npz + nOrder)*D.Nx + D.npx + nOrder], sizeof(float), D.nx, fp);
    }
  fclose(fp);
/*
  FILE *fp = fopen(filename, "wt");
    for(int i = 0; i < D.nz; i++) {
      for(int j = 0; j < D.nx; j++)
        fprintf(fp, "%lf, ", (double)P[(i + D.npz + nOrder)*D.Nx + D.npx + nOrder + j]);
      fprintf(fp, "\n");
    }
  fclose(fp);
*/
}

void wave_exe(wave W, media M, dim D, coeff C) {
  int ix, iz, idx;
  int sidx = D.npx + D.sx + nOrder - 1 + (D.npz + D.sz + nOrder - 1)*D.Nx;
  float srclet;
  float *dpmlx, *dpmlz, *lambda, *mu;
  float *factor[10];

  int i;
  int memSize = D.Nx*D.Nz*sizeof(float);
  dpmlx  = (float*) malloc(memSize);
  dpmlz  = (float*) malloc(memSize);
  lambda = (float*) malloc(memSize);
  mu     = (float*) malloc(memSize);
  for(i = 0; i < 10; i++)
    factor[i] = (float*) malloc(memSize);

  for(idx = 0; idx < D.Nx*D.Nz; idx++) {
    mu[idx] = M.rho[idx]*M.vs[idx]*M.vs[idx];
    lambda[idx] = M.rho[idx]*M.vp[idx]*M.vp[idx] - 2*mu[idx];
  }
  for(iz = 0; iz < D.Nz; iz++)
    for(ix = 0; ix < D.Nx; ix++) {
      idx = iz*D.Nx + ix;
      if(ix < D.npx + nOrder && ix >= nOrder) 
        dpmlx[idx] = C.d0x*pow(1.0*(D.npx + nOrder - ix)/D.npx, 2);
      else if(ix >= D.Nx - D.npx - nOrder && ix < D.Nx - nOrder) 
        dpmlx[idx] = C.d0x*pow(1.0*(ix + D.npx + nOrder + 1 - D.Nx)/D.npx, 2);
      else
        dpmlx[idx] = 0.0;

      if(iz < D.npz + nOrder && iz >= nOrder) 
        dpmlz[idx] = C.d0z*pow(1.0*(D.npz + nOrder - iz)/D.npz, 2);
      else if(iz >= D.Nz - D.npz - nOrder && iz < D.Nz - nOrder) 
        dpmlz[idx] = C.d0z*pow(1.0*(iz + D.npz + nOrder + 1 - D.Nz)/D.npz, 2);
      else
        dpmlz[idx] = 0.0;
    }

  for(idx = 0; idx < D.Nx*D.Nz; idx++) {
    factor[0][idx] = (2 - C.dt*dpmlx[idx])/(2 + C.dt*dpmlx[idx]);
    factor[1][idx] = (2 - C.dt*dpmlz[idx])/(2 + C.dt*dpmlz[idx]);
    factor[2][idx] = 2*C.dt/(2 + C.dt*dpmlx[idx])/M.rho[idx]/D.dx;
    factor[3][idx] = 2*C.dt/(2 + C.dt*dpmlz[idx])/M.rho[idx]/D.dz;
    factor[4][idx] = 2*C.dt/(2 + C.dt*dpmlx[idx])*(lambda[idx] + 2*mu[idx])/D.dx;
    factor[5][idx] = 2*C.dt/(2 + C.dt*dpmlz[idx])*lambda[idx]/D.dz;
    factor[6][idx] = 2*C.dt/(2 + C.dt*dpmlx[idx])*lambda[idx]/D.dx;
    factor[7][idx] = 2*C.dt/(2 + C.dt*dpmlz[idx])*(lambda[idx] + 2*mu[idx])/D.dz;
    factor[8][idx] = 2*C.dt/(2 + C.dt*dpmlx[idx])*mu[idx]/D.dx;
    factor[9][idx] = 2*C.dt/(2 + C.dt*dpmlz[idx])*mu[idx]/D.dz;
  }

  for(idx = 0; idx < D.Nx*D.Nz; idx++) {
    W.vxx [idx] = 0.0; W.vxz [idx] = 0.0; W.vxt [idx] = 0.0;
    W.vzx [idx] = 0.0; W.vzz [idx] = 0.0; W.vzt [idx] = 0.0;
    W.txxx[idx] = 0.0; W.txxz[idx] = 0.0; W.txxt[idx] = 0.0;
    W.tzzx[idx] = 0.0; W.tzzz[idx] = 0.0; W.tzzt[idx] = 0.0;
    W.txzx[idx] = 0.0; W.txzz[idx] = 0.0; W.txzt[idx] = 0.0;
  }

  int it;
  float Psum;
  char file[200];
  for(it = 0; it < C.nt; it++) {
    if(it%nTimePreSnap == 0) printf("calculating and exporting for step %5d ...\n", it);

    for(iz = nOrder; iz < D.Nz - nOrder; iz++)
      for(ix = nOrder; ix < D.Nx - nOrder; ix++) {
        idx = iz*D.Nx + ix;
        Psum = 0.0;
        for(i = 0; i < nOrder; i++)
          Psum += C.C[i]*(W.vxt[idx + i + 1] - W.vxt[idx - i]);
        W.txxx[idx] = factor[0][idx]*W.txxx[idx] + factor[4][idx]*Psum;
        W.tzzx[idx] = factor[0][idx]*W.tzzx[idx] + factor[6][idx]*Psum;
        Psum = 0.0;
        for(i = 0; i < nOrder; i++)
          Psum += C.C[i]*(W.vzt[idx + i*D.Nx] - W.vzt[idx - (i + 1)*D.Nx]);
        W.txxz[idx] = factor[1][idx]*W.txxz[idx] + factor[5][idx]*Psum;
        W.tzzz[idx] = factor[1][idx]*W.tzzz[idx] + factor[7][idx]*Psum;
      }

    for(iz = nOrder; iz < D.Nz - nOrder; iz++)
      for(ix = nOrder; ix < D.Nx - nOrder; ix++) {
        idx = iz*D.Nx + ix;
        Psum = 0.0;
        for(i = 0; i < nOrder; i++)
          Psum += C.C[i]*(W.vzt[idx + i] - W.vzt[idx - i - 1]);
        W.txzx[idx] = factor[0][idx]*W.txzx[idx] + factor[8][idx]*Psum;
        Psum = 0.0;
        for(i = 0; i < nOrder; i++)
          Psum += C.C[i]*(W.vxt[idx + (i + 1)*D.Nx] - W.vxt[idx - i*D.Nx]);
        W.txzz[idx] = factor[1][idx]*W.txzz[idx] + factor[9][idx]*Psum;
      }
    for(idx = 0; idx < D.Nx*D.Nz; idx++)
      W.txzt[idx] = W.txzx[idx] + W.txzz[idx];

    srclet = (1 - 2*pow(PI*C.f0*((C.dt*it) - C.t0), 2))*exp(- pow(PI*C.f0*(C.dt*it - C.t0), 2));
    W.txxx[sidx] += C.ampl*srclet/4;
    W.txxz[sidx] += C.ampl*srclet/4;
    W.tzzx[sidx] += C.ampl*srclet/4;
    W.tzzz[sidx] += C.ampl*srclet/4;
    for(idx = 0; idx < D.Nx*D.Nz; idx++) {
      W.txxt[idx] = W.txxx[idx] + W.txxz[idx];
      W.tzzt[idx] = W.tzzx[idx] + W.tzzz[idx];
    }
//
    for(iz = nOrder; iz < D.Nz - nOrder; iz++)
      for(ix = nOrder; ix < D.Nx - nOrder; ix++) {
        idx = iz*D.Nx + ix;
        Psum = 0.0;
        for(i = 0; i < nOrder; i++)
          Psum += C.C[i]*(W.txxt[idx + i] - W.txxt[idx - i - 1]);
        W.vxx[idx] = factor[0][idx]*W.vxx[idx] + factor[2][idx]*Psum;
        Psum = 0.0;
        for(i = 0; i < nOrder; i++)
          Psum += C.C[i]*(W.txzt[idx + i*D.Nx] - W.txzt[idx - (i + 1)*D.Nx]);
        W.vxz[idx] = factor[1][idx]*W.vxz[idx] + factor[3][idx]*Psum;
      }
    for(idx = 0; idx < D.Nx*D.Nz; idx++)
      W.vxt[idx] = W.vxx[idx] + W.vxz[idx];

    for(iz = nOrder; iz < D.Nz - nOrder; iz++)
      for(ix = nOrder; ix < D.Nx - nOrder; ix++) {
        idx = iz*D.Nx + ix;
        Psum = 0.0;
        for(i = 0; i < nOrder; i++)
          Psum += C.C[i]*(W.txzt[idx + i + 1] - W.txzt[idx - i]);
        W.vzx[idx] = factor[0][idx]*W.vzx[idx] + factor[2][idx]*Psum;
        Psum = 0.0;
        for(i = 0; i < nOrder; i++)
          Psum += C.C[i]*(W.tzzt[idx + (i + 1)*D.Nx] - W.tzzt[idx - i*D.Nx]);
        W.vzz[idx] = factor[1][idx]*W.vzz[idx] + factor[3][idx]*Psum;
      }
    for(idx = 0; idx < D.Nx*D.Nz; idx++)
      W.vzt[idx] = W.vzx[idx] + W.vzz[idx];

    if(it%nTimePreSnap == 0) {
      sprintf(file, "P%05d.bin", it);
      wave_exp(D, file, W.txxt);
    }
  }

  free(dpmlx); free(dpmlz);
  free(lambda); free(mu);
  for(i = 0; i < 10; i++)
    free(factor[i]);
}

int main(int argc, char *argv[]) {

  int nx = 500, nz = 600;
  int npmlx = 20, npmlz = 20;
  int sx = 80, sz = 80;
  float dx = 5.0, dz = 5.0;
  int nt = 1000;
  float dt = 1.0e-3;
  int nppw = 12;
  float ampl = 1.0e0;

  wave W; media M; dim D; coeff C;

  int Nx, Nz;
  int memSize;
  int i;
  
  cout << "Input nt = ";
  cin >> nt;

  int prod1, prod2;
  for(int m = 1; m < nOrder + 1; m++) {
    prod1 = 1;
    for(i = 1; i <= nOrder; i++)
      if(i != m) prod1 *= (2*i - 1)*(2*i - 1);
    prod2 = 1;
    for(i = 1; i <= nOrder; i++)
      if(i != m) prod2 *= abs((2*m - 1)*(2*m - 1) - (2*i - 1)*(2*i - 1));
    C.C[m - 1] = pow(-1.0, m + 1)*prod1/(2*m - 1)/prod2;
  }

  Nx = nx + 2*npmlx + 2*nOrder;
  Nz = nz + 2*npmlz + 2*nOrder;
  memSize = Nx*Nz*sizeof(float);

  D.nx = nx; D.nz = nz;
  D.Nx = Nx; D.Nz = Nz;
  D.sx = sx; D.sz = sz;
  D.npx = npmlx; D.npz = npmlz;
  D.dx = dx; D.dz= dz;
  C.nt = nt; C.dt = dt;
  C.ampl = ampl;

  M.vp  = (float*) malloc(memSize);
  M.vs  = (float*) malloc(memSize);
  M.rho = (float*) malloc(memSize);
  for(i = 0; i < Nx*Nz; i++) {
    M.vp [i] = 2000.0;
    M.vs [i] = 1000.0;
    M.rho[i] = 1000.0;
  }

  float *vsmin = M.vs, *vsmax = M.vs;
  for(i = 1; i < Nx*Nz; i++) {
    if(*vsmin > M.vs[i]) vsmin = &M.vs[i];
    if(*vsmax < M.vs[i]) vsmax = &M.vs[i];
  }

  C.f0 = (*vsmin)/(max(dx, dz)*nppw);
  C.t0 = 1.0/C.f0;
  C.d0x = 3*(*vsmax)/dx*(8.0/15 - 3.0/100*npmlx + 1.0/1500*npmlx*npmlx);
  C.d0z = 3*(*vsmax)/dz*(8.0/15 - 3.0/100*npmlz + 1.0/1500*npmlz*npmlz);

  W.vxx  = (float*) malloc(memSize);
  W.vxz  = (float*) malloc(memSize);
  W.vzx  = (float*) malloc(memSize);
  W.vzz  = (float*) malloc(memSize);
  W.txxx = (float*) malloc(memSize);
  W.txxz = (float*) malloc(memSize);
  W.tzzx = (float*) malloc(memSize);
  W.tzzz = (float*) malloc(memSize);
  W.txzx = (float*) malloc(memSize);
  W.txzz = (float*) malloc(memSize);
  W.vxt  = (float*) malloc(memSize);
  W.vzt  = (float*) malloc(memSize);
  W.txxt = (float*) malloc(memSize);
  W.tzzt = (float*) malloc(memSize);
  W.txzt = (float*) malloc(memSize);

  wave_exe(W, M, D, C);

  free(M.vp); free(M.vs); free(M.rho);

  free(W.vxx ); free(W.vxz ); free(W.vxt );
  free(W.vzx ); free(W.vzz ); free(W.vzt );
  free(W.txxx); free(W.txxz); free(W.txxt);
  free(W.tzzx); free(W.tzzz); free(W.tzzt);
  free(W.txzx); free(W.txzz); free(W.txzt);

  return 0;
}

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
  float dt;
  float d0x, d0z;
  float C[nOrder];
} coeff;
typedef struct {
  float *f0, *f1, *f2, *f3, *f4, *f5, *f6, *f7, *f8, *f9;
} factor;

__global__ void pre_eval(wave W, media M, dim D, coeff C, factor F) {
  int ix = threadIdx.x + blockIdx.x*blockDim.x;
  int iz = threadIdx.y + blockIdx.y*blockDim.y;
  int idx = iz*D.Nx + ix;
  float dpmlx = 0.0, dpmlz = 0.0;
  float lambda, mu;

  if(ix < D.Nx && iz < D.Nz) {
    mu = M.rho[idx]*M.vs[idx]*M.vs[idx];
    lambda = M.rho[idx]*M.vp[idx]*M.vp[idx] - 2*mu;
    if(ix < D.npx + nOrder && ix >= nOrder)
      dpmlx = C.d0x*pow(1.0*(D.npx + nOrder - ix)/D.npx, 2);
    if(ix >= D.Nx - D.npx - nOrder && ix < D.Nx - nOrder)
      dpmlx = C.d0x*pow(1.0*(ix + D.npx + nOrder + 1 - D.Nx)/D.npx, 2);
    if(iz < D.npz + nOrder && iz >= nOrder)
      dpmlz = C.d0z*pow(1.0*(D.npz + nOrder - iz)/D.npz, 2);
    if(iz >= D.Nz - D.npz - nOrder && iz < D.Nz - nOrder)
      dpmlz = C.d0z*pow(1.0*(iz + D.npz + nOrder + 1 - D.Nz)/D.npz, 2);

    F.f0[idx] = (2 - C.dt*dpmlx)/(2 + C.dt*dpmlx);
    F.f1[idx] = (2 - C.dt*dpmlz)/(2 + C.dt*dpmlz);
    F.f2[idx] = 2*C.dt/(2 + C.dt*dpmlx)/M.rho[idx]/D.dx;
    F.f3[idx] = 2*C.dt/(2 + C.dt*dpmlz)/M.rho[idx]/D.dz;
    F.f4[idx] = 2*C.dt/(2 + C.dt*dpmlx)*(lambda + 2*mu)/D.dx;
    F.f5[idx] = 2*C.dt/(2 + C.dt*dpmlz)*lambda/D.dz;
    F.f6[idx] = 2*C.dt/(2 + C.dt*dpmlx)*lambda/D.dx;
    F.f7[idx] = 2*C.dt/(2 + C.dt*dpmlz)*(lambda + 2*mu)/D.dz;
    F.f8[idx] = 2*C.dt/(2 + C.dt*dpmlx)*mu/D.dx;
    F.f9[idx] = 2*C.dt/(2 + C.dt*dpmlz)*mu/D.dz;

    W.vxx [idx] = 0.0; W.vxz [idx] = 0.0; W.vxt [idx] = 0.0;
    W.vzx [idx] = 0.0; W.vzz [idx] = 0.0; W.vzt [idx] = 0.0;
    W.txxx[idx] = 0.0; W.txxz[idx] = 0.0; W.txxt[idx] = 0.0;
    W.tzzx[idx] = 0.0; W.tzzz[idx] = 0.0; W.tzzt[idx] = 0.0;
    W.txzx[idx] = 0.0; W.txzz[idx] = 0.0; W.txzt[idx] = 0.0;
  }
}

__global__ void vel_eval(wave W, dim D, coeff C, factor F, int sidx) {
  int ix = threadIdx.x + blockIdx.x*blockDim.x;
  int iz = threadIdx.y + blockIdx.y*blockDim.y;
  int idx = iz*D.Nx + ix;
  int i;
  float Psum;
  
  if(ix >= nOrder && ix < D.Nx - nOrder && iz >= nOrder && iz < D.Nz - nOrder) {
    Psum = 0.0;
    for(i = 0; i < nOrder; i++)
      Psum += C.C[i]*(W.txxt[idx + i] - W.txxt[idx - i - 1]);
    W.vxx[idx] = F.f0[idx]*W.vxx[idx] + F.f2[idx]*Psum;
    Psum = 0.0;
    for(i = 0; i < nOrder; i++)
      Psum += C.C[i]*(W.txzt[idx + i*D.Nx] - W.txzt[idx - (i + 1)*D.Nx]);
    W.vxz[idx] = F.f1[idx]*W.vxz[idx] + F.f3[idx]*Psum;
    W.vxt[idx] = W.vxx[idx] + W.vxz[idx];

    Psum = 0.0;
    for(i = 0; i < nOrder; i++)
      Psum += C.C[i]*(W.txzt[idx + i + 1] - W.txzt[idx - i]);
    W.vzx[idx] = F.f0[idx]*W.vzx[idx] + F.f2[idx]*Psum;
    Psum = 0.0;
    for(i = 0; i < nOrder; i++)
      Psum += C.C[i]*(W.tzzt[idx + (i + 1)*D.Nx] - W.tzzt[idx - i*D.Nx]);
    W.vzz[idx] = F.f1[idx]*W.vzz[idx] + F.f3[idx]*Psum;
    W.vzt[idx] = W.vzx[idx] + W.vzz[idx];
  }
}

__global__ void str_eval(wave W, dim D, coeff C, int sidx, float srclet, factor F) {
  int ix = threadIdx.x + blockIdx.x*blockDim.x;
  int iz = threadIdx.y + blockIdx.y*blockDim.y;
  int idx = iz*D.Nx + ix;
  int i;
  float Psum;

  if(ix >= nOrder && ix < D.Nx - nOrder && iz >= nOrder && iz < D.Nz - nOrder) {
    Psum = 0.0;
    for(i = 0; i < nOrder; i++)
      Psum += C.C[i]*(W.vxt[idx + i + 1] - W.vxt[idx - i]);
    W.txxx[idx] = F.f0[idx]*W.txxx[idx] + F.f4[idx]*Psum;
    W.tzzx[idx] = F.f0[idx]*W.tzzx[idx] + F.f6[idx]*Psum;
    Psum = 0.0;
    for(i = 0; i < nOrder; i++)
      Psum += C.C[i]*(W.vzt[idx + i*D.Nx] - W.vzt[idx - (i + 1)*D.Nx]);
    W.txxz[idx] = F.f1[idx]*W.txxz[idx] + F.f5[idx]*Psum;
    W.tzzz[idx] = F.f1[idx]*W.tzzz[idx] + F.f7[idx]*Psum;

    Psum = 0.0;
    for(i = 0; i < nOrder; i++)
      Psum += C.C[i]*(W.vzt[idx + i] - W.vzt[idx - i - 1]);
    W.txzx[idx] = F.f0[idx]*W.txzx[idx] + F.f8[idx]*Psum;
    Psum = 0.0;
    for(i = 0; i < nOrder; i++)
      Psum += C.C[i]*(W.vxt[idx + (i + 1)*D.Nx] - W.vxt[idx - i*D.Nx]);
    W.txzz[idx] = F.f1[idx]*W.txzz[idx] + F.f9[idx]*Psum;
    W.txzt[idx] = W.txzx[idx] + W.txzz[idx];

    if(idx == sidx) {
      W.txxx[idx] += srclet/4;
      W.txxz[idx] += srclet/4;
      W.tzzx[idx] += srclet/4;
      W.tzzz[idx] += srclet/4;
    }
    W.txxt[idx] = W.txxx[idx] + W.txxz[idx];
    W.tzzt[idx] = W.tzzx[idx] + W.tzzz[idx];
  }
}

void exp_wave(dim D, char *filename, float *P) {
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

int main(int argc, char *argv[]) {

  int nx = 500, nz = 600;
  int npmlx = 20, npmlz = 20;
  int sx = 80, sz = 80;
  float dx = 5.0, dz = 5.0;
  int nt = 1000;
  float dt = 1.0e-3;
  int nppw = 12;
  float ampl = 1.0e0;

  wave W; media M; dim D; coeff C; factor F;

  int Nx, Nz;
  size_t memSize;
  int i, j;

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

  float *Vp  = (float*) malloc(memSize);
  float *Vs  = (float*) malloc(memSize);
  float *Rho = (float*) malloc(memSize);
  for(i = 0; i < Nx*Nz; i++) {
    Vp [i] = 2000.0;
    Vs [i] = 1000.0;
    Rho[i] = 1000.0;
  }

  cudaMalloc((float**) &M.vp , memSize);
  cudaMalloc((float**) &M.vs , memSize);
  cudaMalloc((float**) &M.rho, memSize);
  cudaMemcpy(M.vp , Vp , memSize, cudaMemcpyHostToDevice);
  cudaMemcpy(M.vs , Vs , memSize, cudaMemcpyHostToDevice);
  cudaMemcpy(M.rho, Rho, memSize, cudaMemcpyHostToDevice);

  float *vsmin = Vs, *vsmax = Vs;
  for(i = 1; i < Nx*Nz; i++) {
    if(*vsmin > Vs[i]) vsmin = &Vs[i];
    if(*vsmax < Vs[i]) vsmax = &Vs[i];
  }

  float f0, t0;
  f0 = (*vsmin)/(max(dx, dz)*nppw);
  t0 = 1.0/f0;
  C.dt = dt;
  C.d0x = 3*(*vsmax)/dx*(8.0/15 - 3.0/100*npmlx + 1.0/1500*npmlx*npmlx);
  C.d0z = 3*(*vsmax)/dz*(8.0/15 - 3.0/100*npmlz + 1.0/1500*npmlz*npmlz);

  cudaMalloc((float**) &F.f0, memSize);
  cudaMalloc((float**) &F.f1, memSize);
  cudaMalloc((float**) &F.f2, memSize);
  cudaMalloc((float**) &F.f3, memSize);
  cudaMalloc((float**) &F.f4, memSize);
  cudaMalloc((float**) &F.f5, memSize);
  cudaMalloc((float**) &F.f6, memSize);
  cudaMalloc((float**) &F.f7, memSize);
  cudaMalloc((float**) &F.f8, memSize);
  cudaMalloc((float**) &F.f9, memSize);

  cudaMalloc((float**) &W.vxx , memSize);
  cudaMalloc((float**) &W.vxz , memSize);
  cudaMalloc((float**) &W.vzx , memSize);
  cudaMalloc((float**) &W.vzz , memSize);
  cudaMalloc((float**) &W.txxx, memSize);
  cudaMalloc((float**) &W.txxz, memSize);
  cudaMalloc((float**) &W.tzzx, memSize);
  cudaMalloc((float**) &W.tzzz, memSize);
  cudaMalloc((float**) &W.txzx, memSize);
  cudaMalloc((float**) &W.txzz, memSize);
  cudaMalloc((float**) &W.vxt , memSize);
  cudaMalloc((float**) &W.vzt , memSize);
  cudaMalloc((float**) &W.txxt, memSize);
  cudaMalloc((float**) &W.tzzt, memSize);
  cudaMalloc((float**) &W.txzt, memSize);

  dim3 Block(32, 16);
  dim3 Grid(ceil(1.0*Nx/Block.x), ceil(1.0*Nz/Block.y));
  
  cout << "Block = " << Block.x << " " << Block.y << endl;
  cout << "Grid = " << Grid.x << " " << Grid.y << endl;

  float *P;
  P = (float*) malloc(memSize);

  float srclet;
  int sidx = npmlx + sx + nOrder - 1 + (npmlz + sz + nOrder - 1)*D.Nx;
  pre_eval <<< Grid, Block >>> (W, M, D, C, F);

  int it = 0;
  char file[200];
  for(i = 0; i*nTimePreSnap < nt; i++) {
    printf("calculating and exporting for step %5d ...\n", it);
    for(j = 0; j < nTimePreSnap; j++, it++) {
      if(it > nt) break;
      srclet = ampl*(1 - 2*pow(PI*f0*(dt*it - t0), 2))*exp( - pow(PI*f0*(dt*it - t0), 2));

      str_eval <<< Grid, Block >>> (W, D, C, sidx, srclet, F);
      vel_eval <<< Grid, Block >>> (W, D, C, F, sidx);

      if((it - 1)%nTimePreSnap == 0) {
        cudaDeviceSynchronize();
        sprintf(file, "P%05d.bin", it - 1);
        cudaMemcpy(P, W.txxt, memSize, cudaMemcpyDeviceToHost);
        exp_wave(D, file, P);
      }
    }
  }

  free(Vp); free(Vs); free(Rho);
  cudaFree(M.vp); cudaFree(M.vs); cudaFree(M.rho);

  free(P);

  cudaFree(W.vxx ); cudaFree(W.vxz ); cudaFree(W.vxt );
  cudaFree(W.vzx ); cudaFree(W.vzz ); cudaFree(W.vzt );
  cudaFree(W.txxx); cudaFree(W.txxz); cudaFree(W.txxt);
  cudaFree(W.tzzx); cudaFree(W.tzzz); cudaFree(W.tzzt);
  cudaFree(W.txzx); cudaFree(W.txzz); cudaFree(W.txzt);

  return 0;
}

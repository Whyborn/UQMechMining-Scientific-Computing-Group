/*
Vortex visualisation using the Triple Decomposition Method

@author: Nick Gibbons
*/

#include <stdio.h>
#include <math.h>
#include "nmsimplex.h"
#include "lapacke.h"

#define pi 3.141592653589793

// Using global things here so that refscalar can have only 1 argument
double Sn[3][3];
double Wn[3][3];
double scratch[1000]; // Also I don't want this reallocated every iteration
static int lenscratch = 1000;
static double guess[8][3] = {{ 045.0/180.0*pi, 045.0/180.0*pi, 005.0/180.0*pi}, 
                             { 135.0/180.0*pi, 045.0/180.0*pi, 005.0/180.0*pi}, 
                             { 225.0/180.0*pi, 045.0/180.0*pi, 005.0/180.0*pi}, 
                             { 315.0/180.0*pi, 045.0/180.0*pi, 005.0/180.0*pi}, 
                             { 045.0/180.0*pi,-045.0/180.0*pi, 005.0/180.0*pi}, 
                             { 135.0/180.0*pi,-045.0/180.0*pi, 005.0/180.0*pi}, 
                             { 225.0/180.0*pi,-045.0/180.0*pi, 005.0/180.0*pi}, 
                             { 315.0/180.0*pi,-045.0/180.0*pi, 005.0/180.0*pi}};

void ComputeQ(double a, double b, double c, double Q[3][3]){
    // Compute the rotations matrix from kolar_vort07 appendix A (Checked)
    Q[0][0] = cos(a)*cos(b)*cos(c) - sin(a)*sin(c);
    Q[0][1] = sin(a)*cos(b)*cos(c) + cos(a)*sin(c);
    Q[0][2] =                      - sin(b)*cos(c);
    Q[1][0] =-cos(a)*cos(b)*sin(c) - sin(a)*cos(c);
    Q[1][1] =-sin(a)*cos(b)*sin(c) + cos(a)*cos(c);
    Q[1][2] =                        sin(b)*sin(c);
    Q[2][0] = cos(a)*sin(b);
    Q[2][1] = sin(a)*sin(b);
    Q[2][2] = cos(b);
    return;
}

void ComputeT(double a[3][3], double aT[3][3]){
    // Compute matrix transverse by copying
    int i,j;
    for (i=0; i<3; i++){
        for (j=0; j<3; j++){
            aT[j][i] = a[i][j];
        }
    }
    return;
}

void matmul(double a[3][3], double b[3][3], double c[3][3]){
    // Unrolled 3x3 only matrix multiplication (Checked)
    
    c[0][0] = a[0][0]*b[0][0] + a[0][1]*b[1][0] + a[0][2]*b[2][0];
    c[0][1] = a[0][0]*b[0][1] + a[0][1]*b[1][1] + a[0][2]*b[2][1];
    c[0][2] = a[0][0]*b[0][2] + a[0][1]*b[1][2] + a[0][2]*b[2][2];

    c[1][0] = a[1][0]*b[0][0] + a[1][1]*b[1][0] + a[1][2]*b[2][0];
    c[1][1] = a[1][0]*b[0][1] + a[1][1]*b[1][1] + a[1][2]*b[2][1];
    c[1][2] = a[1][0]*b[0][2] + a[1][1]*b[1][2] + a[1][2]*b[2][2];
    
    c[2][0] = a[2][0]*b[0][0] + a[2][1]*b[1][0] + a[2][2]*b[2][0];
    c[2][1] = a[2][0]*b[0][1] + a[2][1]*b[1][1] + a[2][2]*b[2][1];
    c[2][2] = a[2][0]*b[0][2] + a[2][1]*b[1][2] + a[2][2]*b[2][2];
    return;
}

double refscalar(double angles[3]){ 
    // Compute the Reference Frame Scalar (kolar_vort07 equation 10) (Checked)
    double a,b,c, REF, temp[3][3];
    double Q[3][3], Sdash[3][3], Wdash[3][3], QT[3][3];

    a = angles[0];
    b = angles[1];
    c = angles[2];
    ComputeQ(a,b,c,Q);
    ComputeT(Q,QT);

    matmul(Q,Sn,temp);
    matmul(temp,QT,Sdash);

    matmul(Q,Wn,temp);
    matmul(temp,QT,Wdash);

    REF = fabs(Sdash[0][1]*Wdash[0][1]);
    REF+= fabs(Sdash[1][2]*Wdash[1][2]);
    REF+= fabs(Sdash[2][0]*Wdash[2][0]);
    return -REF/1e6;// Divide to keep numbers nicer? Waste?
}

void eigvals3x3(double* A, double* reig, double* ieig, int* ier){
    // Wrapper function for LAPACK dgeev
    char doL, doR;
    int lenL, lenR, N;
    double vL[1], vR[1];

    doL = 'N'; doR = 'N';
    lenL = 1; lenR = 1;
    N=3;

    dgeev_(&doL, &doR, &N, A, &N, reig, ieig, vL, &lenL, vR, &lenR, scratch, &lenscratch, ier);
    return;
}

void sort3(double* a){
    // In place sort of an array of 3 numbers (GPU branching trouble ??)
    double temp;
    if (a[0] > a[1]){
        temp = a[0];
        a[0] = a[1];
        a[1] = temp;
	}
    if (a[0] > a[2]){
        temp = a[0];
        a[0] = a[2];
        a[2] = temp;
	}
    if (a[1] > a[2]){
        temp = a[1];
        a[1] = a[2];
        a[2] = temp;
	}
    return;
}

void ComputeL2(double* du, int Nx, double* L, int verbose){
    /* Interface to Python Function: Do everything here
    du: Double pointer to numpy array of velocity grads (nx, neq, neq) (IN)
    Nx: Int of nx (number of cells) (IN)
    L : Double to numpy array of lambda 2 values (nx) (OUT)

    Notes: You should probably return the angles instead maybe????
    */
    int i,j,ier,n,p;
    double LL[3][3],ievals[3],revals[3],temp[3][3];

    p=Nx/100;
    if (p==0) p=1; // Check to make sure that the percent thing works

    for (n=0; n<Nx; n++){
        if ((n%p==0)&&(verbose==1)) {
            printf("@ %f percent\n", n*100.0/Nx);
        }

        // Copy data to local arrays for caching purposes
        for (i=0; i<3; i++){
            for (j=0; j<3; j++){
                Sn[i][j] = 0.5*(du[9*n + 3*i + j] + du[9*n + 3*j + i]);
                Wn[i][j] = 0.5*(du[9*n + 3*i + j] - du[9*n + 3*j + i]);
            }
        }

        // Now compute the lambda 2 criterion
        matmul(Sn,Sn,temp);
        matmul(Wn,Wn,LL);

        for (i=0; i<3; i++) for (j=0; j<3; j++) LL[i][j] += temp[i][j]; 

        eigvals3x3(&LL[0][0], revals, ievals,&ier); // A little trick here to pass &LL
        if (ier!=0){
            printf("Error computing eigenvalues: %i\n", ier);
            exit(1);
		}
        sort3(revals);
        L[n] = revals[1];
    }
    return;
}

void ComputeL2TDM(double* du, int Nx, double* L, int verbose){
    /* Interface to Python Function: Do everything here
    S : Double pointer to numpy array of strains (nx, neq, neq) (IN)
    W : Double pointer to numpy array of Vorticities (nx, neq, neq) (IN)
    du: Double pointer to numpy array of velocity grads (nx, neq, neq) (IN)
    Nx: Int of nx (number of cells) (IN)
    L : Double to numpy array of lambda 2 values (nx) (OUT)

    Notes: You should probably return the angles instead maybe????
    */
    int i,j,n,p,ier,nc;
    double result,a,b,c;
    double dun[3][3], temp[3][3], resdash[3][3], LL[3][3];
    double dudash[3][3], el[3][3], rr[3][3], Q[3][3], QT[3][3];
    double angles[3], revals[3], ievals[3], res[3][3];

    p=Nx/100;
    if (p==0) p=1; // Check to make sure that the percent thing works

    for (n=0; n<Nx; n++){
        if ((n%p==0)&&(verbose==1)) {
            printf("@ %f percent\n", n*100.0/Nx);
        }

        // Copy data to local arrays for caching purposes
        for (i=0; i<3; i++){
            for (j=0; j<3; j++){
                dun[i][j] = du[9*n + 3*i + j];
                Sn[i][j] = 0.5*(du[9*n + 3*i + j] + du[9*n + 3*j + i]);
                Wn[i][j] = 0.5*(du[9*n + 3*i + j] - du[9*n + 3*j + i]);
            }
        }

        nc=1;
        while(nc>0){
            angles[0] = guess[nc-1][0]; angles[1] = guess[nc-1][1]; angles[2] = guess[nc-1][2];
            // Yeah this is going to be a problem... There are multiple solutions...
            // you may need to work out a test for this, or use a repeated sequence
            // of guessing and take the lowest value

            // Use Nelder-Mead to find the Basic reference frame angles (BRF)
            result=simplex(refscalar, angles, 3, 1e-7, 1.0, NULL);
            //printf("result:%f a: %f b: %f c: %f\n",result, angles[0], angles[1], angles[2]);

            if (isnan(result)){ // it didn't work
                if (nc>9) { // Don't try again
                    printf("Fault in N %i, solver not converged\n",n);
                    printf("du[9] = {%.17g,\n %.17g,\n %.17g,\n",dun[0][0],dun[0][1],dun[0][2]);
                    printf("%.17g,\n %.17g,\n %.17g,\n",dun[1][0],dun[1][1],dun[1][2]);
                    printf("%.17g,\n %.17g,\n %.17g};",dun[2][0],dun[2][1],dun[2][2]);
                    exit(1);
                }
                else{ //Try again
                    nc++;
                }
            }
            else { // It worked!
                nc=0;
            }
        }
       
        a = angles[0]; b = angles[1]; c = angles[2];

        // Now compute the Shear tensor and return to lab frame
        ComputeQ(a,b,c,Q);
        ComputeT(Q,QT);
        matmul(Q,dun,temp);
        matmul(temp,QT,dudash);  

        // Construct the residual tensor in the BRF (kolar_vort07 equation 8b)
        for (i=0; i<3; i++){
            for (j=0; j<3; j++){
                if (i==j) {
                  resdash[i][j] = dudash[i][j];
                }
                else{
                  result = fmin(fabs(dudash[i][j]),fabs(dudash[j][i]));
                  resdash[i][j]=copysign(result, dudash[i][j]);
                }
            }
        }

        // Subtract it from the velocity gradients to get the sheartensor in the BRF
        //for (i=0; i<3; i++){ // Just to check?
        //    for (i=0; j<3; j++){
        //        sheardash[i][j] = dudash[i][j] - resdash[i][j];
        //    }
        //}
        
        // Now transform the residual tensor back into the lab frame
        matmul(QT,resdash,temp);
        matmul(temp,Q,res); // Since inverse(Q) == QT;
        
        for (i=0; i<3; i++){
            for (j=0; j<3; j++){
                el[i][j] = 0.5*(res[i][j] + res[j][i]);
                rr[i][j] = 0.5*(res[i][j] - res[j][i]);
            }
        }
        
        // Now compute the lambda 2 criterion
        matmul(el,el,temp);
        matmul(rr,rr,LL);

        for (i=0; i<3; i++){
            for (j=0; j<3; j++){
                LL[i][j] += temp[i][j]; 
            }
        }
        eigvals3x3(&LL[0][0], revals, ievals,&ier); // A little trick here to pass &LL
        if (ier!=0){
            printf("Error computing eigenvalues: %i\n", ier);
            exit(1);
		}
        sort3(revals);
        L[n] = revals[1];
    }
    return;
}


int main(int argc, char* argv []){
    //  test
    int i,j;
    double L[1], S[9], W[9];

    //double du[9] ={-1764368.18259121,  842278.96366864, 3872913.7410053 ,
    //                30442.6142948 ,  629865.44086344,  410334.15700657,
    //               -467807.16789635, -181200.90148082,  985354.4824681 };
    //double S[9] = {-1764368.18259121,  436360.78898172, 1702553.28655447,
    //                436360.78898172,  629865.44086344,  114566.62776288,
    //                1702553.28655447,  114566.62776288,  985354.4824681 };

    //double W[9] = {       0.        ,  405918.17468692, 2170360.45445083,
    //               -405918.17468692,       0.        ,  295767.5292437 ,
    //               -2170360.45445083, -295767.5292437 ,       0.        };

    //double du[9] = { -1.31361036e+05,  -9.48324458e+06,  -3.80935454e+07,
    //           1.61793813e+05,  -1.58325629e+07,  -6.28565086e+07,
    //          -4.77490594e+04,   2.59194235e+06,   1.01567858e+07};
    double du[9] = {-131361.03641460655,
     -9483244.5785234664,
      -38093545.411565095,
      161793.81333122187,
       -15832562.879114239,
        -62856508.621331587,
        -47749.059393181145,
         2591942.3547066413,
          10156785.756732743};

    for (i=0; i<3; i++){
        for (j=0; j<3; j++){
            S[3*i+j] = 0.5*(du[3*i+j] + du[3*j+i]);
            W[3*i+j] = 0.5*(du[3*i+j] - du[3*j+i]);
        }
    }

    ComputeL2TDM(du, 1, L,0);

    double answer = 2.39385944e+09; // From serial_vortid.py

    if ((L[0]-answer)/answer*100>5e-5){
        printf("Unit test failed!\n");
        printf("Returned L2 of: %f  Answer: %f\n", L[0], answer);
    }
    else{
        printf("Test Passed!\n");
        printf("Returned L2 of: %f  Answer: %f\n", L[0], answer);
    }

    return 0;
}

/*
Simple integrator in c

Compile with
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void integrate_particle(double* x, double* xd, double* xdd, double* xs, int N, int dim, double dt){
    // Simple fixed timestep integrator
    for(int i=0; i<N; i++) {
        for(int d=0; d<dim; d++) {
            x[d] += xd[d]*dt + 0.5*xdd[d]*dt*dt;
            xd[d] += xdd[d]*dt;
            xs[i*dim + d] = x[d];
        }
    }
    return;
}


#include <omp.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef int64_t Int;
typedef double Real;
double alpha   = 0.15;
double epsilon = 1e-6;
Int n_seeds    = 512;

void ppr(double* p ,Int seed, double alpha, double epsilon, Int n_nodes, Int* indptr, Int* indices, Int* degrees) {
    double* r       = (double*)malloc(n_nodes * sizeof(double));
    double* r_prime = (double*)malloc(n_nodes * sizeof(double));
    for(Int i = 0; i < n_nodes; i++) {
        p[i]       = 0;
        r[i]       = 0;
        r_prime[i] = 0;
    }
    
    r[seed] = 1;
    
    // --
    // Run
    
    Int* frontier = (Int*)malloc(n_nodes * sizeof(Int));
    for(Int i = 0; i < n_nodes; i++) 
        frontier[i] = -1;
    
    frontier[0] = seed;
    Int frontier_size = 1;
    
    while(1) {
        if(frontier_size == 0) break;
        
        memcpy(r_prime, r, n_nodes * sizeof(Int));
        
        for(Int i = 0; i < frontier_size; i++) {
            Int node_idx = frontier[i];
            p[node_idx] += (2 * alpha) / (1 + alpha) * r[node_idx];
            r_prime[node_idx] = 0;
        }
        
        for(Int i = 0; i < frontier_size; i++) {
            Int src_idx = frontier[i];
            Int deg     = degrees[src_idx];
            Int offset  = indptr[src_idx];
            for(Int j = 0; j < deg; j++) {
                Int dst_idx   = indices[offset + j];
                double update = ((1 - alpha) / (1 + alpha)) * r[src_idx] / deg;
                r_prime[dst_idx] += update;
            }
        }
        
        memcpy(r, r_prime, n_nodes * sizeof(Int));
        
        frontier_size = 0;
        for(Int i = 0; i < n_nodes; i++) {
            if(
                (degrees[i] > 0) &&
                (r[i] >= degrees[i] * epsilon)
            ) {
                frontier[frontier_size] = i;
                frontier_size++;    
            }
        }
    }
}

int main(int argc, char *argv[]) {
    
    // -- 
    // IO
    
    FILE *fptr = NULL; 
    fptr = fopen("data/jhu.bin", "rb");    
    Int n_nodes;
    fread(&n_nodes, sizeof(Int), 1, fptr);
    
    Int n_edges;
    fread(&n_edges, sizeof(Int), 1, fptr);
        
    Int* indptr = (Int*)malloc((n_nodes + 1) * sizeof(Int));
    fread(indptr, sizeof(Int), n_nodes + 1, fptr);
    
    Int* indices = (Int*)malloc(n_edges * sizeof(Int));
    fread(indices, sizeof(Int), n_edges, fptr);

    Int* degrees = (Int*)malloc(n_nodes * sizeof(Int));
    for(Int i = 0; i < n_nodes; i++)
        degrees[i] = indptr[i + 1] - indptr[i];
    
    // --
    // Run
    
    double* P = (double*)malloc(n_nodes * n_seeds * sizeof(double)); // Array to n_seeds * n_nodes dense output
    
    #pragma omp parallel for
    for(Int seed = 0; seed < n_seeds; seed++) {
        ppr(
            &P[seed * n_nodes],
            seed,
            alpha,
            epsilon,
            n_nodes,
            indptr,
            indices,
            degrees
        );
    }

    for(Int seed = 0; seed < n_seeds; seed++) {
        double pp = 0;
        for(Int i = 0; i < n_nodes; i++) {
            pp += P[seed * n_nodes + i];
        }
        printf("pp %ld %0.16f\n", seed, pp);
    }
    
    return 0;
}

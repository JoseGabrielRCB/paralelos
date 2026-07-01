#ifndef GRAFO_H
#define GRAFO_H

#include <stdint.h>
#include <stddef.h>
#include <stdio.h>

typedef struct Aresta {
    uint32_t u;
    uint32_t v;
    double   peso;
} Aresta;

typedef struct DSU {
    uint32_t *pai;
    uint32_t *rank;
    uint32_t  n;
} DSU;

int dsu_init(DSU *dsu, uint32_t n);

void dsu_free(DSU *dsu);

uint32_t dsu_find(DSU *dsu, uint32_t x);

int dsu_union(DSU *dsu, uint32_t a, uint32_t b);

uint32_t dsu_num_componentes(DSU *dsu);

Aresta *ler_grafo_texto(const char *caminho, uint32_t *V, uint64_t *E);

Aresta *ler_grafo_binario(const char *caminho, uint32_t *V, uint64_t *E);

FILE *grafo_bin_abrir(const char *caminho);

uint64_t grafo_bin_ler_bloco(FILE *f, Aresta *buf, uint64_t cap);

int grafo_bin_info(const char *caminho, uint32_t *V, uint64_t *E);

Aresta *ler_grafo_mpiio(const char *caminho, uint32_t *V, uint64_t *E_total,
                        int rank, int nprocs, uint64_t *E_local,
                        uint64_t *byte_ini, uint64_t *byte_fim);

Aresta *ler_grafo_binario_mpiio(const char *caminho, uint32_t *V, uint64_t *E_total,
                                int rank, int nprocs, uint64_t *E_local,
                                uint64_t *rec_ini, uint64_t *rec_fim);

Aresta *ler_grafo_binario_dist(const char *caminho, uint32_t *V, uint64_t *E_total,
                               int rank, int nprocs, uint64_t *E_local,
                               uint64_t *rec_ini, uint64_t *rec_fim);

#endif

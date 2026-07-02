/* Base do projeto AGM (Boruvka): tipos, DSU e leitores de I/O. */
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

/* DSU / Union-Find: raiz(x) = componente de x. */
typedef struct DSU {
    uint32_t *pai;
    uint32_t *rank;
    double   *peso;
    uint32_t  n;
} DSU;

int dsu_init(DSU *dsu, uint32_t n);

void dsu_free(DSU *dsu);

uint32_t dsu_find(DSU *dsu, uint32_t x);

int dsu_union(DSU *dsu, uint32_t a, uint32_t b, double peso_aresta);

/* [SOMA PARCIAL] Peso acumulado do componente ao qual 'x' pertence. */
double dsu_peso_componente(DSU *dsu, uint32_t x);

uint32_t dsu_num_componentes(DSU *dsu);

/* --- Leitores sequenciais (carregam o grafo inteiro na RAM) --- */
/* Texto: linha1=V, linha2=E, demais "u v peso". */
Aresta *ler_grafo_texto(const char *caminho, uint32_t *V, uint64_t *E);
/* Binario: sem cabecalho; E=tamanho/16, V=(maior indice)+1. */
Aresta *ler_grafo_binario(const char *caminho, uint32_t *V, uint64_t *E);

/* --- Leitura binaria em STREAMING (memoria O(bloco), p/ .bin gigante) --- */
FILE *grafo_bin_abrir(const char *caminho);
uint64_t grafo_bin_ler_bloco(FILE *f, Aresta *buf, uint64_t cap);
/* Passe unico so p/ descobrir V e E sem carregar o grafo. */
int grafo_bin_info(const char *caminho, uint32_t *V, uint64_t *E);

/* --- Leitores paralelos (implementados so com -DHAVE_MPI; senao stub NULL) --- */
/* Texto via MPI-IO: cada rank le sua fatia de bytes com alinhamento de linha. */
Aresta *ler_grafo_mpiio(const char *caminho, uint32_t *V, uint64_t *E_total,
                        int rank, int nprocs, uint64_t *E_local,
                        uint64_t *byte_ini, uint64_t *byte_fim);

/* Binario via MPI-IO COLETIVO: exige o arquivo visivel em todos os nos (NFS). */
Aresta *ler_grafo_binario_mpiio(const char *caminho, uint32_t *V, uint64_t *E_total,
                                int rank, int nprocs, uint64_t *E_local,
                                uint64_t *rec_ini, uint64_t *rec_fim);

/* Binario SEM disco compartilhado: rank 0 le e envia as fatias por MPI_Send.
 * Mesma particao do coletivo -> resultado identico. E o leitor padrao. */
Aresta *ler_grafo_binario_dist(const char *caminho, uint32_t *V, uint64_t *E_total,
                               int rank, int nprocs, uint64_t *E_local,
                               uint64_t *rec_ini, uint64_t *rec_fim);

#endif

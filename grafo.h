/*
 * grafo.h — Estruturas base do projeto de Árvore Geradora Mínima (Borůvka).
 *
 * FASE 0: apenas a BASE (tipos, DSU e protótipos de I/O). O algoritmo de
 * Borůvka em si NÃO é implementado aqui.
 *
 * Desvios autorizados refletidos nesta base (ver ALGORTIMO_REFERENCIA.md):
 *   - D1: o grafo é representado por LISTA DE ARESTAS (struct Aresta), e não
 *         pela matriz de pesos W(u,v) do pseudocódigo literal — a matriz é
 *         inviável (~1,96M vértices => ~16 TB de RAM).
 *   - D2: peso da aresta em uint32_t; a soma da MST é acumulada em uint64_t.
 *
 * Convenções:
 *   - Índices de vértice em BASE 0 (0 .. V-1). V é lido da 1ª linha dos dados.
 */
#ifndef GRAFO_H
#define GRAFO_H

#include <stdint.h>
#include <stddef.h>

/* ------------------------------------------------------------------------- *
 * Aresta — uma aresta não direcionada ponderada.
 *   u, v : vértices (base 0).
 *   peso : uint32_t (D2 — pesos passam de 4×10^9? não; passam de 4×10^9 só na
 *          soma. O peso individual cabe em uint32_t pois vai além de int com
 *          sinal mas dentro de 2^32-1).
 * ------------------------------------------------------------------------- */
typedef struct Aresta {
    uint32_t u;
    uint32_t v;
    uint32_t peso;
} Aresta;

/* ------------------------------------------------------------------------- *
 * DSU (Disjoint Set Union / Union-Find) — implementa a função `raiz` do
 * pseudocódigo: dá a raiz (componente) à qual cada vértice pertence.
 *   pai  : floresta de representantes (pai[i] == i  =>  i é raiz).
 *   rank : limite superior da altura da árvore (union by rank).
 *   n    : número de vértices (tamanho dos vetores).
 * ------------------------------------------------------------------------- */
typedef struct DSU {
    uint32_t *pai;
    uint32_t *rank;
    uint32_t  n;
} DSU;

/* Inicializa o DSU com n conjuntos singleton (pai[i]=i, rank[i]=0).
 * Retorna 1 em sucesso, 0 em falha de alocação. */
int dsu_init(DSU *dsu, uint32_t n);

/* Libera a memória interna do DSU (pode ser chamada com dsu já zerado). */
void dsu_free(DSU *dsu);

/* Encontra a raiz do conjunto de x, com COMPRESSÃO DE CAMINHO. */
uint32_t dsu_find(DSU *dsu, uint32_t x);

/* Une os conjuntos de a e b por UNION BY RANK.
 * Retorna 1 se houve fusão (estavam em conjuntos distintos),
 *         0 se já pertenciam ao mesmo conjunto (fusão redundante). */
int dsu_union(DSU *dsu, uint32_t a, uint32_t b);

/* Conta quantos componentes (conjuntos distintos) existem no DSU. */
uint32_t dsu_num_componentes(DSU *dsu);

/* ------------------------------------------------------------------------- *
 * Protótipos de I/O (STUBS na Fase 0 — implementados em fases posteriores).
 *
 * Formato do arquivo de dados:
 *   Linha 1: V (número de vértices)
 *   Linha 2: E (número de arestas)
 *   Linhas seguintes: "u v peso"
 * ------------------------------------------------------------------------- */

/* Lê o grafo de um arquivo texto (versão sequencial).
 * Aloca e retorna o vetor de arestas; preenche *V e *E.
 * Retorna NULL em caso de erro. (STUB — Fase 0) */
Aresta *ler_grafo_texto(const char *caminho, uint32_t *V, uint64_t *E);

/* Lê/distribui o grafo via MPI-IO (versão paralela, Estilo B).
 *
 * Comportamento:
 *   - rank 0 lê o CABEÇALHO (V na 1ª linha, E na 2ª) e o difunde (MPI_Bcast),
 *     junto com o offset de byte onde começam as arestas (início da linha 3);
 *   - cada rank lê APENAS sua fatia contígua de BYTES com MPI_File_read_at_all,
 *     com alinhamento de bordas para nunca cortar uma linha ao meio;
 *   - cada rank faz o parsing da sua fatia para um vetor local de Aresta,
 *     ignorando linhas em branco/mal-formadas (há 1 linha vazia no fim).
 *
 * Parâmetros de saída:
 *   *V, *E_total : cabeçalho (confiamos no E do cabeçalho);
 *   *E_local     : nº de arestas efetivamente carregadas por ESTE rank;
 *   *byte_ini/*byte_fim : intervalo de bytes lido por este rank (para o log).
 *
 * Retorno: vetor local de arestas (NULL em erro).
 *
 * OBS: a assinatura é deliberadamente livre de tipos MPI para que grafo.h possa
 * ser incluído também pelos alvos compilados com gcc (teste/sequencial). A
 * IMPLEMENTAÇÃO só existe quando grafo.c é compilado com -DHAVE_MPI (mpicc);
 * caso contrário, é um stub que retorna NULL. Internamente usa MPI_COMM_WORLD. */
Aresta *ler_grafo_mpiio(const char *caminho, uint32_t *V, uint64_t *E_total,
                        int rank, int nprocs, uint64_t *E_local,
                        uint64_t *byte_ini, uint64_t *byte_fim);

#endif /* GRAFO_H */

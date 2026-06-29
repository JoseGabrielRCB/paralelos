/*
 * grafo.h — Estruturas base do projeto de Árvore Geradora Mínima (Borůvka).
 *
 * FASE 0: apenas a BASE (tipos, DSU e protótipos de I/O). O algoritmo de
 * Borůvka em si NÃO é implementado aqui.
 *
 * Desvios autorizados refletidos nesta base (ver ALGORTIMO_REFERENCIA.md):
 *   - D1: o grafo é representado por LISTA DE ARESTAS (struct Aresta), e não
 *         pela matriz de pesos W(u,v) do pseudocódigo literal — a matriz é
 *         inviável (~10M vértices => ~800 TB de RAM).
 *   - D2: peso da aresta em double (ponto flutuante); a soma da MST também
 *         é acumulada em double. (Antes era uint32_t/uint64_t, mas o arquivo
 *         de execução real `graph.bin` traz pesos em ponto flutuante.)
 *
 * Convenções:
 *   - Índices de vértice em BASE 0 (0 .. V-1). V é lido da 1ª linha dos dados.
 */
#ifndef GRAFO_H
#define GRAFO_H

#include <stdint.h>
#include <stddef.h>
#include <stdio.h>

/* ------------------------------------------------------------------------- *
 * Aresta — uma aresta não direcionada ponderada.
 *   u, v : vértices (base 0), uint32_t.
 *   peso : double (D2 — o arquivo de execução `graph.bin` armazena o peso como
 *          ponto flutuante de 64 bits).
 *
 * IMPORTANTE — layout casado com o disco: o arquivo binário grava cada aresta
 * como u(uint32) + v(uint32) + peso(double) = 16 bytes contíguos. Esta struct
 * tem EXATAMENTE esse layout (u@0, v@4, peso@8, sizeof==16, sem padding), então
 * a leitura binária é uma cópia direta de bytes, sem parsing campo a campo.
 * ------------------------------------------------------------------------- */
typedef struct Aresta {
    uint32_t u;
    uint32_t v;
    double   peso;
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
 * Formato texto: linha 1 = V, linha 2 = E, demais "u v peso".
 * Retorna NULL em caso de erro. */
Aresta *ler_grafo_texto(const char *caminho, uint32_t *V, uint64_t *E);

/* Lê o grafo de um arquivo BINÁRIO (versão sequencial).
 *
 * Formato binário (o de `graph.bin`): SEM cabeçalho; o arquivo é uma sequência
 * de registros contíguos de 16 bytes, cada um:
 *     u    : uint32 little-endian
 *     v    : uint32 little-endian
 *     peso : double (8 bytes) little-endian
 * Logo  E = tamanho_do_arquivo / 16. Como não há V no cabeçalho, V é deduzido
 * como (maior índice de vértice visto) + 1.
 *
 * Aloca e retorna o vetor de arestas (cópia direta dos bytes); preenche *V e *E.
 * Retorna NULL em caso de erro. */
Aresta *ler_grafo_binario(const char *caminho, uint32_t *V, uint64_t *E);

/* ------------------------------------------------------------------------- *
 * Leitura BINÁRIA em STREAMING (memória O(bloco)).
 *
 * Para quando o arquivo NÃO cabe na RAM (ex.: graph.bin, 12,8 GB numa máquina
 * com 3,5 GB). Em vez de carregar todas as arestas, o algoritmo re-lê o arquivo
 * em BLOCOS a cada fase de Borůvka. O uso de memória passa a ser O(tamanho do
 * bloco) + as estruturas dimensionadas por V (DSU, melhor-aresta-por-componente).
 * ------------------------------------------------------------------------- */

/* Abre o arquivo binário para leitura sequencial (em blocos). NULL em erro. */
FILE *grafo_bin_abrir(const char *caminho);

/* Lê o próximo bloco de até 'cap' arestas para 'buf'. Retorna o nº de arestas
 * lidas (0 indica fim de arquivo). */
uint64_t grafo_bin_ler_bloco(FILE *f, Aresta *buf, uint64_t cap);

/* Faz UM passe em streaming para descobrir E (= tamanho/16) e V (= maior índice
 * de vértice + 1), sem carregar o grafo. Retorna 1 em sucesso, 0 em erro. */
int grafo_bin_info(const char *caminho, uint32_t *V, uint64_t *E);

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
 *   byte_ini, byte_fim : intervalo de bytes lido por este rank (para o log).
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

/* Lê/distribui o grafo BINÁRIO via MPI-IO (versão paralela, Estilo B).
 *
 * Muito mais simples que a versão texto: como cada registro tem 16 bytes FIXOS,
 * a partição é aritmética exata, sem alinhamento de bordas de linha. O rank i
 * é dono dos registros [i*E/np, (i+1)*E/np); cada rank lê só essa fatia com
 * MPI_File_read_at_all e reinterpreta os bytes como vetor de Aresta.
 *
 * Como o binário não tem cabeçalho:
 *   - E_total = tamanho_do_arquivo / 16 (calculado via MPI_File_get_size);
 *   - V = (maior índice de vértice global) + 1, obtido com MPI_Allreduce(MAX)
 *     do maior índice local de cada rank.
 *
 * Parâmetros de saída:
 *   *V, *E_total : deduzidos do arquivo (ver acima);
 *   *E_local     : nº de arestas (registros) carregadas por ESTE rank;
 *   rec_ini, rec_fim : intervalo de REGISTROS [rec_ini, rec_fim) deste rank.
 *
 * Retorno: vetor local de arestas (NULL em erro). Stub (retorna NULL) quando
 * compilado sem -DHAVE_MPI. */
Aresta *ler_grafo_binario_mpiio(const char *caminho, uint32_t *V, uint64_t *E_total,
                                int rank, int nprocs, uint64_t *E_local,
                                uint64_t *rec_ini, uint64_t *rec_fim);

#endif /* GRAFO_H */

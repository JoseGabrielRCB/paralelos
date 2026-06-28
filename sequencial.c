/*
 * sequencial.c — Árvore Geradora Mínima por BORŮVKA (forma LITERAL da Parte I
 *                do ALGORITMO_REFERENCIA.md).
 *
 * Transcrição do pseudocódigo (Parte I):
 *     Inicialize floresta F = (V, E') com E' = {}
 *     concluído := false
 *     enquanto não concluído faça
 *         Encontre os componentes conectados de F e atribua a cada vértice seu componente
 *         Inicialize a aresta de menor peso de cada componente como "Nenhuma"
 *         para cada aresta uv em E, com u e v em componentes diferentes:
 *             se for-preferido-sobre(uv, menor[comp(u)]): menor[comp(u)] := uv
 *             se for-preferido-sobre(uv, menor[comp(v)]): menor[comp(v)] := uv
 *         se todos os componentes têm menor == "Nenhuma":  concluído := true
 *         senão: para cada componente com menor != "Nenhuma": adicione menor em E'
 *
 * Desvios autorizados (anotados no código):
 *   - D1: usamos LISTA DE ARESTAS (struct Aresta) no lugar da matriz W(u,v),
 *         inviável para ~1,96M vértices (~16 TB). A semântica de W(u,v) e do
 *         MIN por componente é preservada: o "menor por componente" é calculado
 *         pela varredura da lista de arestas.
 *   - D2: peso em uint32_t; soma total da MST em uint64_t.
 *   - D3 (parada): a Parte I termina LITERALMENTE quando "todos os componentes
 *         têm menor == Nenhuma" (nenhuma árvore pode ser mesclada). Para o grafo
 *         DESCONEXO isso ocorre antes de virar 1 só componente — o algoritmo
 *         então encerra retornando a FLORESTA geradora mínima. Detectamos esse
 *         caso (componentes > 1) e avisamos que a condição "1 componente" (forma
 *         citada na Parte II) NÃO foi atendida por o grafo ser desconexo.
 *         GUARDA adicional: se uma fase não realizar nenhuma fusão, encerramos.
 *
 * Uso:  ./sequencial <arquivo_dados> [arquivo_saida]
 *   - <arquivo_dados>: caminho passado por argv[1] (NÃO hardcodado).
 *   - [arquivo_saida]: opcional (padrão "mst_sequencial.txt").
 *
 * Saídas:
 *   - stdout: "Peso total da MST: <soma_uint64>"
 *   - arquivo de saída: as arestas escolhidas, uma "u v peso" por linha.
 *   - stderr: tempo de I/O e tempo de cálculo (clock_gettime CLOCK_MONOTONIC).
 */
#define _POSIX_C_SOURCE 199309L

#include "grafo.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <inttypes.h>

/* Sentinela para "aresta de menor peso = Nenhuma".
 *
 * IMPORTANTE: este é um sentinela de ÍNDICE de aresta (um valor de índice
 * inválido), NÃO um sentinela de PESO. Os índices válidos vão de 0 a E-1
 * (E ~ milhões), nunca UINT32_MAX, então não há colisão possível com nenhum
 * peso de aresta — inclusive entradas de teste com peso == UINT32_MAX. Por
 * isso o sequencial é robusto e dá resultado IDÊNTICO ao paralelo (que, por
 * trafegar a aresta na rede, usa um campo de validade explícito em vez de
 * índice). Logo, não é preciso alterar a lógica aqui. */
#define ARESTA_NENHUMA UINT32_MAX

/* Diferença em segundos entre dois instantes do CLOCK_MONOTONIC. */
static double seg(struct timespec a, struct timespec b)
{
    return (double) (b.tv_sec - a.tv_sec) + (double) (b.tv_nsec - a.tv_nsec) / 1e9;
}

/* ------------------------------------------------------------------------- *
 * for-preferido-sobre(cand, atual)  — pseudocódigo, Parte I.
 *
 *   return (atual is "None")
 *          ou (peso(cand) < peso(atual))
 *          ou (peso(cand) == peso(atual) e regra-desempate(cand, atual))
 *
 * regra-desempate: ordem TOTAL determinística — em empate de peso, prefere o
 * MENOR par (u,v), tomado de forma normalizada (min(u,v), depois max(u,v)),
 * de modo que a orientação da aresta não importe e o resultado seja reproduzível.
 *
 * 'cand' e 'atual' são índices na lista de arestas; 'atual' pode ser
 * ARESTA_NENHUMA (== "None"). Retorna 1 se 'cand' é preferida sobre 'atual'.
 * ------------------------------------------------------------------------- */
static int for_preferido_sobre(const Aresta *arestas, uint32_t cand, uint32_t atual)
{
    if (atual == ARESTA_NENHUMA)
        return 1; /* atual is "None" */

    uint32_t pc = arestas[cand].peso;
    uint32_t pa = arestas[atual].peso;
    if (pc < pa) return 1; /* peso estritamente menor */
    if (pc > pa) return 0;

    /* Empate de peso -> regra-desempate por menor par (min,max). */
    uint32_t cu = arestas[cand].u, cv = arestas[cand].v;
    uint32_t au = arestas[atual].u, av = arestas[atual].v;
    uint32_t cmin = cu < cv ? cu : cv, cmax = cu < cv ? cv : cu;
    uint32_t amin = au < av ? au : av, amax = au < av ? av : au;

    if (cmin != amin) return cmin < amin;
    return cmax < amax; /* se (min,max) iguais => não preferida (retorna 0) */
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "uso: %s <arquivo_dados> [arquivo_saida]\n", argv[0]);
        return 1;
    }
    const char *entrada = argv[1];
    const char *saida   = (argc >= 3) ? argv[2] : "mst_sequencial.txt";

    struct timespec t0, t1;

    /* ===================== 1) LEITURA (I/O) ===================== */
    uint32_t V = 0;
    uint64_t E = 0;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    Aresta *arestas = ler_grafo_texto(entrada, &V, &E);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double tempo_io = seg(t0, t1);

    if (arestas == NULL || V == 0) {
        fprintf(stderr, "Erro: falha ao ler o grafo.\n");
        free(arestas);
        return 1;
    }

    /* ===================== 2) CÁLCULO (Borůvka) ===================== */
    clock_gettime(CLOCK_MONOTONIC, &t0);

    DSU dsu;
    if (!dsu_init(&dsu, V)) {
        fprintf(stderr, "Erro: falha ao inicializar DSU (V=%" PRIu32 ").\n", V);
        free(arestas);
        return 1;
    }

    /* comp[i]     = componente (raiz) do vértice i, recalculado a cada fase.
     * menor[c]    = índice da aresta de menor peso do componente c ("Nenhuma" se ARESTA_NENHUMA).
     * mst[]       = arestas adicionadas a E' (a floresta resultante). */
    uint32_t *comp  = (uint32_t *) malloc((size_t) V * sizeof(uint32_t));
    uint32_t *menor = (uint32_t *) malloc((size_t) V * sizeof(uint32_t));
    Aresta   *mst   = (Aresta *)   malloc((size_t) V * sizeof(Aresta)); /* <= V-1 arestas */
    if (comp == NULL || menor == NULL || mst == NULL) {
        fprintf(stderr, "Erro: falha de alocacao no calculo.\n");
        free(comp); free(menor); free(mst);
        dsu_free(&dsu); free(arestas);
        return 1;
    }

    uint64_t soma   = 0;  /* D2: soma dos pesos em uint64_t */
    uint64_t mst_n  = 0;  /* número de arestas em E' */
    uint32_t fases  = 0;
    int desconexo_sem_fusao = 0; /* sinaliza a guarda D3 (fase sem fusão) */

    int concluido = 0;
    while (!concluido) {
        fases++;

        /* (a) Componentes conectados de F: cada vértice recebe sua raiz. */
        for (uint32_t i = 0; i < V; i++)
            comp[i] = dsu_find(&dsu, i);

        /* (b) Aresta de menor peso de cada componente := "Nenhuma". */
        for (uint32_t i = 0; i < V; i++)
            menor[i] = ARESTA_NENHUMA;

        /* (c) Para cada aresta uv com u,v em componentes diferentes,
         *     atualiza o menor do componente de u e o do componente de v. */
        for (uint64_t i = 0; i < E; i++) {
            uint32_t cu = comp[arestas[i].u];
            uint32_t cv = comp[arestas[i].v];
            if (cu == cv)
                continue; /* mesma componente: não conecta árvores distintas */

            if (for_preferido_sobre(arestas, (uint32_t) i, menor[cu]))
                menor[cu] = (uint32_t) i;
            if (for_preferido_sobre(arestas, (uint32_t) i, menor[cv]))
                menor[cv] = (uint32_t) i;
        }

        /* (d) "se todos os componentes têm menor == Nenhuma" -> concluído. */
        int algum = 0;
        for (uint32_t c = 0; c < V; c++) {
            if (menor[c] != ARESTA_NENHUMA) { algum = 1; break; }
        }

        if (!algum) {
            /* Forma LITERAL de parada da Parte I: nenhuma árvore pode ser
             * mesclada (não há aresta entre componentes diferentes). */
            concluido = 1;
        } else {
            /* (e) Adiciona a aresta de menor peso de cada componente em E'.
             * O dsu_union evita ciclos: uma aresta escolhida simultaneamente
             * pelos dois extremos (comp(u) e comp(v)), ou que fecharia ciclo
             * dentro da fase, retorna 0 e não é contada — preservando a
             * propriedade de FLORESTA exigida pelo pseudocódigo. */
            uint64_t fusoes = 0;
            for (uint32_t c = 0; c < V; c++) {
                if (menor[c] == ARESTA_NENHUMA)
                    continue;
                Aresta e = arestas[menor[c]];
                if (dsu_union(&dsu, e.u, e.v)) {
                    mst[mst_n++] = e;
                    soma += (uint64_t) e.peso; /* D2 */
                    fusoes++;
                }
            }

            /* GUARDA D3: se havia arestas candidatas mas nenhuma fusão ocorreu,
             * não há progresso possível — encerra retornando a floresta.
             * (Com DSU válido isto não deveria ocorrer; mantido por segurança.) */
            if (fusoes == 0) {
                desconexo_sem_fusao = 1;
                concluido = 1;
            }
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double tempo_calc = seg(t0, t1);

    uint32_t n_comp = dsu_num_componentes(&dsu);

    /* ===================== 3) SAÍDA ===================== */

    /* 3.1 stdout: peso total. */
    printf("Peso total da MST: %" PRIu64 "\n", soma);

    /* 3.2 arquivo: MST como LISTA DE ARESTAS, uma por linha, no formato do
     *     enunciado "<u> → (<peso>) → <v>" (mesmo formato do paralelo.c, para
     *     saídas idênticas). Podem ser milhões -> NÃO vão ao stdout. */
    FILE *fo = fopen(saida, "w");
    if (fo == NULL) {
        fprintf(stderr, "Aviso: nao foi possivel abrir '%s' para escrita.\n", saida);
    } else {
        for (uint64_t i = 0; i < mst_n; i++)
            fprintf(fo, "%" PRIu32 " → (%" PRIu32 ") → %" PRIu32 "\n",
                    mst[i].u, mst[i].peso, mst[i].v);
        fclose(fo);
    }

    /* 3.3 stderr: tempos e diagnóstico. */
    fprintf(stderr, "----------------------------------------\n");
    fprintf(stderr, "Vertices (V).......: %" PRIu32 "\n", V);
    fprintf(stderr, "Arestas  (E).......: %" PRIu64 "\n", E);
    fprintf(stderr, "Fases de Boruvka...: %" PRIu32 "\n", fases);
    fprintf(stderr, "Arestas na floresta: %" PRIu64 "\n", mst_n);
    fprintf(stderr, "Componentes finais.: %" PRIu32 "\n", n_comp);
    fprintf(stderr, "Arestas gravadas em: %s\n", saida);
    fprintf(stderr, "Tempo de I/O.......: %.6f s\n", tempo_io);
    fprintf(stderr, "Tempo de calculo...: %.6f s\n", tempo_calc);
    if (n_comp > 1) {
        /* D3: condição literal "1 componente" (Parte II) NÃO atingida. */
        fprintf(stderr,
            "[D3] Condicao literal '1 componente' NAO atendida: grafo DESCONEXO "
            "(%" PRIu32 " componentes). Retornada a FLORESTA geradora minima.\n", n_comp);
        if (desconexo_sem_fusao)
            fprintf(stderr,
                "[D3] (encerrado pela guarda: fase sem nenhuma fusao)\n");
    }
    fprintf(stderr, "----------------------------------------\n");

    /* ===================== 4) LIBERAÇÃO ===================== */
    free(comp); free(menor); free(mst);
    dsu_free(&dsu);
    free(arestas);
    return 0;
}

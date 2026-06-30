/*
 * paralelo.c — Árvore Geradora Mínima por BORŮVKA PARALELO (forma LITERAL da
 *              Parte II do ALGORITMO_REFERENCIA.md), no ESTILO B (coletivas MPI).
 *
 * -------------------------------------------------------------------------
 * Pseudocódigo de referência (Parte II) e o mapeamento usado aqui:
 *
 *   para cada u ∈ V faça em paralelo: raiz(i) := i        -> dsu_init (DSU replicado)
 *   enquanto not(finalizado):
 *     para cada u: vertice+proximo(u) := v com raiz(u)!=raiz(v) e W(u,v) mínimo
 *                                                          -> "melhor aresta por componente"
 *     para cada componente k: escolha o vértice de menor W(u,vertice+proximo(u))
 *                                                          -> reduz ao MIN por componente (best[k])
 *     combine os novos vértices e crie a floresta Fi+1     -> aplica as uniões no DSU
 *     se raiz(u)==raiz(v) para todo u,v: finalizado:=True  -> parada (1 componente) + GUARDA D3
 *
 * "raiz" é Union-Find/DSU; "vertice+proximo(u)" é o vértice de OUTRA árvore mais
 * próximo de u (aresta de menor peso que sai do componente de u).
 *
 * -------------------------------------------------------------------------
 * Desvios autorizados (todos anotados):
 *   - D1: lista de arestas (struct Aresta) no lugar da matriz de pesos W(u,v),
 *         inviável para ~1,96M vértices (~16 TB). O "MIN{W(u,w)}" por componente
 *         vira a varredura das arestas locais buscando a de menor peso que liga
 *         o componente a outro componente.
 *   - D2: peso em uint32_t; soma total da MST em uint64_t.
 *   - D3 (parada): o critério literal "raiz(u)==raiz(v) para todos" (1 só
 *         componente) é implementado; GUARDA: se uma fase não fizer NENHUMA
 *         fusão e ainda houver >1 componente (grafo desconexo), encerra
 *         retornando a FLORESTA geradora mínima — a condição literal (1
 *         componente) NÃO é atingida e isso é avisado.
 *   - Estilo B (coletivas): o "faça em paralelo" é feito por processos MPI que
 *         particionam as arestas; a combinação dos mínimos por componente usa
 *         MPI_Allreduce + MPI_Op CUSTOMIZADA (mesmo desempate determinístico),
 *         de modo que TODOS os ranks fiquem com o vetor global idêntico e
 *         apliquem as MESMAS uniões na MESMA ordem (resultado reproduzível,
 *         independente do nº de processos). Leitura via MPI-IO (ler_grafo_mpiio).
 *
 * O DSU é REPLICADO e idêntico em todos os ranks; ele NÃO trafega pela rede.
 * Só trafegam: o cabeçalho (Bcast) e o vetor "melhor aresta por componente"
 * (Allreduce) a cada fase.
 *
 * Sem threads. Uso:  mpirun -np <N> ./paralelo <arquivo_dados> [arquivo_saida]
 *   - <arquivo_dados> por argv[1] (NÃO hardcodado);
 *   - [arquivo_saida] opcional (padrão "mst_paralelo.txt").
 *
 * Saídas:
 *   - stdout (rank 0): "Peso total da MST: <soma_uint64>"
 *   - arquivo (rank 0): arestas da MST e o "caminho" no formato do enunciado
 *   - stderr (rank 0): tempo TOTAL e tempo SÓ de cálculo (MPI_Wtime + barreiras)
 *   - log_rank<rank>.txt: diagnóstico por rank (arestas locais, bytes, fases...)
 */
#include <mpi.h>
#include "grafo.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <stddef.h>   /* offsetof */
#include <string.h>
#include <ctype.h>

/* Estrutura reduzida que trafega no Allreduce: a "melhor aresta" de um componente.
 *
 * "Nenhuma" (componente sem aresta de saída) é representado por um CAMPO DE
 * VALIDADE explícito (valido == 0), e NÃO por um peso sentinela mágico. Isso é
 * robusto: o maior peso real (4.294.967.269) fica a apenas 26 de UINT32_MAX, e
 * entradas de teste desconhecidas podem conter peso == UINT32_MAX — um sentinela
 * de peso confundiria essa aresta com "Nenhuma". Com 'valido' não há ambiguidade.
 *
 * Layout (atenção a padding): 'valido' (uint8) + 'peso' (double) + 2x uint32.
 * O double força alinhamento de 8 bytes, logo há padding após 'valido'. O tipo
 * MPI é montado com offsetof de cada campo e MPI_Type_create_resized para fixar
 * a extensão em sizeof(ArestaMin), garantindo o avanço correto entre elementos. */
typedef struct {
    uint8_t  valido; /* 0 = "Nenhuma" (sem aresta de saída); 1 = aresta válida */
    double   peso;   /* D2 — peso em ponto flutuante (como em graph.bin) */
    uint32_t u;
    uint32_t v;
} ArestaMin;

/* ------------------------------------------------------------------------- *
 * Desempate determinístico (== for-preferido-sobre / regra-desempate):
 * 'x' é preferida sobre 'y' sse:
 *     x é válida e (y é inválida
 *                   ou peso(x) < peso(y)
 *                   ou (peso igual e par (min,max) de x < par (min,max) de y)).
 * Ordem TOTAL entre arestas válidas -> resultado reproduzível. Uma aresta
 * inválida (valido==0, "Nenhuma") nunca é preferida; e qualquer aresta válida
 * é preferida sobre uma inválida.
 * ------------------------------------------------------------------------- */
/* Verdadeiro se 'nome' termina em 'suf' (case-insensitive): escolhe o leitor
 * (".bin" => binário/graph.bin; caso contrário => texto via MPI-IO de linhas). */
static int tem_sufixo(const char *nome, const char *suf)
{
    size_t ln = strlen(nome), ls = strlen(suf);
    if (ln < ls) return 0;
    const char *p = nome + (ln - ls);
    for (size_t i = 0; i < ls; i++)
        if (tolower((unsigned char) p[i]) != tolower((unsigned char) suf[i]))
            return 0;
    return 1;
}

static inline int aresta_preferida(const ArestaMin *x, const ArestaMin *y)
{
    if (!x->valido) return 0; /* x é "Nenhuma": nunca preferida */
    if (!y->valido) return 1; /* y é "Nenhuma" e x é válida: x preferida */
    if (x->peso != y->peso)
        return x->peso < y->peso;
    uint32_t xmin = x->u < x->v ? x->u : x->v;
    uint32_t xmax = x->u < x->v ? x->v : x->u;
    uint32_t ymin = y->u < y->v ? y->u : y->v;
    uint32_t ymax = y->u < y->v ? y->v : y->u;
    if (xmin != ymin)
        return xmin < ymin;
    return xmax < ymax; /* se (min,max) iguais -> não preferida (retorna 0) */
}

/* ------------------------------------------------------------------------- *
 * Função da MPI_Op customizada: para cada elemento, mantém em 'inout' a aresta
 * preferida entre 'in' e 'inout'. É comutativa e associativa (é um MIN sob uma
 * ordem total), então pode ser registrada com commute=1.
 * ------------------------------------------------------------------------- */
static void op_aresta_min(void *in_, void *inout_, int *len, MPI_Datatype *dt)
{
    (void) dt;
    const ArestaMin *in = (const ArestaMin *) in_;
    ArestaMin *inout = (ArestaMin *) inout_;
    for (int i = 0; i < *len; i++)
        if (aresta_preferida(&in[i], &inout[i]))
            inout[i] = in[i];
}

/* ------------------------------------------------------------------------- *
 * Escreve (apenas no rank 0) o arquivo de saída: a MST como LISTA DE ARESTAS,
 * uma por linha, no formato do enunciado:
 *
 *     <u> → (<peso>) → <v>
 *
 * Uma aresta da MST por linha (sem caminhada Euler). O MESMO formato é usado
 * em sequencial.c, de modo que as saídas dos dois programas sejam idênticas.
 * Mantido em arquivo (não no stdout), pois podem ser ~V-1 arestas (milhões).
 * ------------------------------------------------------------------------- */
static void escrever_saida_rank0(const char *saida, const Aresta *mst,
                                 uint64_t mst_n)
{
    FILE *fo = fopen(saida, "w");
    if (fo == NULL) {
        fprintf(stderr, "Aviso: nao foi possivel abrir '%s' para escrita.\n", saida);
        return;
    }
    for (uint64_t i = 0; i < mst_n; i++)
        fprintf(fo, "%" PRIu32 " → (%.10g) → %" PRIu32 "\n",
                mst[i].u, mst[i].peso, mst[i].v);
    fclose(fo);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 2) {
        if (rank == 0)
            fprintf(stderr, "uso: mpirun -np <N> %s <arquivo_dados> [arquivo_saida]\n", argv[0]);
        MPI_Finalize();
        return 1;
    }
    const char *entrada = argv[1];
    const char *saida   = (argc >= 3) ? argv[2] : "mst_paralelo.txt";

    /* Log por rank. */
    char lognome[64];
    snprintf(lognome, sizeof lognome, "log_rank%d.txt", rank);
    FILE *logf = fopen(lognome, "w"); /* // VERIFICAR: se NULL, seguimos sem log */

    /* =====================================================================
     * D) TEMPO — início do TOTAL (barreira para alinhar todos os ranks).
     * ===================================================================== */
    MPI_Barrier(MPI_COMM_WORLD);
    double t_total_ini = MPI_Wtime();

    /* =====================================================================
     * A) SETUP + LEITURA via MPI-IO (cada rank lê sua fatia de bytes).
     * ===================================================================== */
    uint32_t V = 0;
    uint64_t E_total = 0, E_local = 0, byte_ini = 0, byte_fim = 0;
    int binario = tem_sufixo(entrada, ".bin");
    /* Leitor binário: por padrão usa a versão DISTRIBUÍDA (rank 0 lê e envia as
     * fatias pela rede), que NÃO exige o arquivo em todas as máquinas — adequada
     * a um cluster por SSH sem disco compartilhado. Se houver NFS/disco comum em
     * todos os nós, exporte GRAFO_MPIIO=1 para usar a leitura coletiva MPI-IO
     * (cada rank lê sua fatia em paralelo, mais rápida nesse caso). */
    int usar_mpiio = (getenv("GRAFO_MPIIO") != NULL);
    Aresta *loc = binario
        ? (usar_mpiio
            ? ler_grafo_binario_mpiio(entrada, &V, &E_total, rank, size,
                                      &E_local, &byte_ini, &byte_fim)
            : ler_grafo_binario_dist(entrada, &V, &E_total, rank, size,
                                     &E_local, &byte_ini, &byte_fim))
        : ler_grafo_mpiio(entrada, &V, &E_total, rank, size,
                          &E_local, &byte_ini, &byte_fim);
    if (loc == NULL || V == 0) {
        if (rank == 0)
            fprintf(stderr, "Erro: falha na leitura MPI-IO.\n");
        if (logf) fclose(logf);
        MPI_Finalize();
        return 1;
    }
    if (logf) {
        fprintf(logf, "==== log do rank %d de %d ====\n", rank, size);
        fprintf(logf, "arquivo de entrada.: %s\n", entrada);
        fprintf(logf, "V (cabecalho)......: %" PRIu32 "\n", V);
        fprintf(logf, "E total (cabecalho): %" PRIu64 "\n", E_total);
        fprintf(logf, "arestas locais.....: %" PRIu64 "\n", E_local);
        fprintf(logf, "bytes lidos........: [%" PRIu64 ", %" PRIu64 ")\n", byte_ini, byte_fim);
        fprintf(logf, "----- fases -----\n");
        fflush(logf);
    }

    /* =====================================================================
     * D) TEMPO — início do CÁLCULO (exclui a leitura/distribuição).
     * ===================================================================== */
    MPI_Barrier(MPI_COMM_WORLD);
    double t_calc_ini = MPI_Wtime();

    /* =====================================================================
     * B) NÚCLEO BORŮVKA (Estilo B).
     * ===================================================================== */

    /* DSU replicado: raiz(i) := i (igual em todos os ranks). */
    DSU dsu;
    if (!dsu_init(&dsu, V)) {
        if (rank == 0) fprintf(stderr, "Erro: dsu_init falhou (V=%" PRIu32 ").\n", V);
        free(loc); if (logf) fclose(logf); MPI_Finalize(); return 1;
    }

    /* best[c] = melhor aresta que sai do componente c (indexado pela raiz). */
    ArestaMin *best = (ArestaMin *) malloc((size_t) V * sizeof(ArestaMin));
    /* mst[] = arestas escolhidas (replicado; <= V-1). */
    Aresta    *mst  = (Aresta *)    malloc((size_t) V * sizeof(Aresta));
    if (best == NULL || mst == NULL) {
        if (rank == 0) fprintf(stderr, "Erro: sem memoria (best/mst).\n");
        free(best); free(mst); free(loc); dsu_free(&dsu);
        if (logf) fclose(logf);
        MPI_Finalize();
        return 1;
    }
    uint64_t mst_n = 0;
    double   soma  = 0.0; /* D2: soma em ponto flutuante */

    /* Tipo MPI para ArestaMin (uint8 valido + double peso + 2x uint32) via offsetof.
     * MPI_Type_create_resized fixa a extensão em sizeof(ArestaMin) para que o
     * Allreduce avance corretamente entre elementos do vetor, mesmo com padding. */
    MPI_Datatype TIPO_ARESTA_MIN;
    {
        int          blk[4]   = {1, 1, 1, 1};
        MPI_Aint     disp[4]  = {offsetof(ArestaMin, valido),
                                 offsetof(ArestaMin, peso),
                                 offsetof(ArestaMin, u),
                                 offsetof(ArestaMin, v)};
        MPI_Datatype tipos[4] = {MPI_UINT8_T, MPI_DOUBLE, MPI_UINT32_T, MPI_UINT32_T};
        MPI_Datatype tmp;
        MPI_Type_create_struct(4, blk, disp, tipos, &tmp);
        MPI_Type_create_resized(tmp, 0, (MPI_Aint) sizeof(ArestaMin), &TIPO_ARESTA_MIN);
        MPI_Type_commit(&TIPO_ARESTA_MIN);
        MPI_Type_free(&tmp);
    }
    /* MPI_Op customizada (comutativa). */
    MPI_Op OP_ARESTA_MIN;
    MPI_Op_create(op_aresta_min, 1 /* commute */, &OP_ARESTA_MIN);

    int      finalizado = 0;
    uint32_t fase = 0;
    int      d3_floresta = 0; /* sinaliza encerramento por grafo desconexo (D3) */

    while (!finalizado) {
        fase++;
        double tf0 = MPI_Wtime();

        /* (1) zera "melhor aresta por componente": valido = 0 ("Nenhuma"). */
        for (uint32_t c = 0; c < V; c++) {
            best[c].valido = 0;
            best[c].peso = 0;
            best[c].u = 0;
            best[c].v = 0;
        }

        /* (2) cada rank percorre SUAS arestas locais e, para (u,v) em
         *     componentes diferentes, atualiza o melhor de raiz(u) e raiz(v). */
        for (uint64_t i = 0; i < E_local; i++) {
            uint32_t ru = dsu_find(&dsu, loc[i].u);
            uint32_t rv = dsu_find(&dsu, loc[i].v);
            if (ru == rv)
                continue; /* mesma árvore: não conecta componentes distintos */

            ArestaMin cand = { 1 /* valido */, loc[i].peso, loc[i].u, loc[i].v };
            if (aresta_preferida(&cand, &best[ru])) best[ru] = cand;
            if (aresta_preferida(&cand, &best[rv])) best[rv] = cand;
        }

        /* (3) Allreduce com a MPI_Op customizada: combina os mínimos locais de
         *     todos os ranks. Após isto, TODOS têm 'best' global idêntico.
         *     MPI_IN_PLACE: usa o próprio 'best' como envio e recepção. */
        MPI_Allreduce(MPI_IN_PLACE, best, (int) V, TIPO_ARESTA_MIN,
                      OP_ARESTA_MIN, MPI_COMM_WORLD);

        /* (4) todos aplicam as MESMAS uniões na MESMA ordem (c crescente).
         *     dsu_union retorna 0 quando a aresta fecharia ciclo (ex.: a mesma
         *     aresta escolhida pelos dois componentes que ela liga, ou um par
         *     de componentes que já se uniu nesta fase) — evitando dupla
         *     contagem e preservando a FLORESTA. Determinístico -> idêntico em
         *     todos os ranks. */
        uint64_t fusoes = 0;
        for (uint32_t c = 0; c < V; c++) {
            if (!best[c].valido)
                continue; /* componente sem aresta de saída nesta fase */
            uint32_t u = best[c].u, v = best[c].v;
            if (dsu_union(&dsu, u, v)) {
                mst[mst_n].u = u;
                mst[mst_n].v = v;
                mst[mst_n].peso = best[c].peso;
                mst_n++;
                soma += best[c].peso; /* D2 */
                fusoes++;
            }
        }

        /* (5) parada — Estilo B: reduz coletivamente o flag "houve fusão?".
         *     fusoes é idêntico em todos os ranks (aplicam as mesmas uniões);
         *     ainda assim usamos MPI_Allreduce (exigência do Estilo B e
         *     proteção contra divergência acidental). */
        uint64_t fusoes_global = 0;
        MPI_Allreduce(&fusoes, &fusoes_global, 1, MPI_UINT64_T, MPI_SUM,
                      MPI_COMM_WORLD);

        double tf1 = MPI_Wtime();
        if (logf) {
            fprintf(logf, "fase %3" PRIu32 ": fusoes=%" PRIu64 " tempo=%.6f s\n",
                    fase, fusoes, tf1 - tf0);
            fflush(logf);
        }

        if (fusoes_global == 0) {
            /* Nenhuma aresta entre componentes diferentes:
             *   - se restou 1 componente -> condição LITERAL atingida;
             *   - se restou >1 componente -> grafo DESCONEXO (GUARDA D3):
             *     encerramos com a floresta geradora mínima. */
            finalizado = 1;
        }
    }

    /* Diagnóstico de componentes (idêntico em todos os ranks). */
    uint32_t ncomp = dsu_num_componentes(&dsu);
    if (ncomp > 1)
        d3_floresta = 1;

    /* =====================================================================
     * D) TEMPO — fim (barreira) e cálculo de TOTAL e SÓ-CÁLCULO.
     * ===================================================================== */
    MPI_Barrier(MPI_COMM_WORLD);
    double t_fim   = MPI_Wtime();
    double t_total = t_fim - t_total_ini; /* inclui leitura/distribuição */
    double t_calc  = t_fim - t_calc_ini;  /* exclui leitura/distribuição */

    if (logf) {
        fprintf(logf, "----- fim -----\n");
        fprintf(logf, "fases totais.......: %" PRIu32 "\n", fase);
        fprintf(logf, "arestas na floresta: %" PRIu64 "\n", mst_n);
        fprintf(logf, "componentes finais.: %" PRIu32 "\n", ncomp);
        fprintf(logf, "tempo TOTAL........: %.6f s\n", t_total);
        fprintf(logf, "tempo SO calculo...: %.6f s\n", t_calc);
        fclose(logf);
    }

    /* =====================================================================
     * C) SAÍDA — apenas o rank 0.
     * ===================================================================== */
    if (rank == 0) {
        /* 1) stdout: peso total (ponto flutuante). */
        printf("Peso total da MST: %.10g\n", soma);

        /* 2) arquivo: MST como lista de arestas (podem ser milhões -> não no stdout). */
        escrever_saida_rank0(saida, mst, mst_n);

        /* 3) stderr: tempos e diagnóstico. */
        fprintf(stderr, "----------------------------------------\n");
        fprintf(stderr, "Processos (np).....: %d\n", size);
        fprintf(stderr, "Vertices (V).......: %" PRIu32 "\n", V);
        fprintf(stderr, "Arestas  (E).......: %" PRIu64 "\n", E_total);
        fprintf(stderr, "Fases de Boruvka...: %" PRIu32 "\n", fase);
        fprintf(stderr, "Arestas na floresta: %" PRIu64 "\n", mst_n);
        fprintf(stderr, "Componentes finais.: %" PRIu32 "\n", ncomp);
        fprintf(stderr, "Arestas/caminho em.: %s\n", saida);
        fprintf(stderr, "Tempo TOTAL........: %.6f s (inclui I/O)\n", t_total);
        fprintf(stderr, "Tempo SO calculo...: %.6f s (exclui I/O)\n", t_calc);
        if (d3_floresta)
            fprintf(stderr,
                "[D3] Condicao literal '1 componente' NAO atendida: grafo DESCONEXO "
                "(%" PRIu32 " componentes). Retornada a FLORESTA geradora minima.\n", ncomp);
        fprintf(stderr, "----------------------------------------\n");
    }

    /* Liberação. */
    MPI_Op_free(&OP_ARESTA_MIN);
    MPI_Type_free(&TIPO_ARESTA_MIN);
    free(best);
    free(mst);
    free(loc);
    dsu_free(&dsu);

    MPI_Finalize();
    return 0;
}

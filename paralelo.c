#include <mpi.h>
#include "grafo.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <stddef.h>
#include <string.h>
#include <ctype.h>

typedef struct {
    uint8_t  valido;
    double   peso;
    uint32_t u;
    uint32_t v;
} ArestaMin;

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
    if (!x->valido) return 0;
    if (!y->valido) return 1;
    if (x->peso != y->peso)
        return x->peso < y->peso;
    uint32_t xmin = x->u < x->v ? x->u : x->v;
    uint32_t xmax = x->u < x->v ? x->v : x->u;
    uint32_t ymin = y->u < y->v ? y->u : y->v;
    uint32_t ymax = y->u < y->v ? y->v : y->u;
    if (xmin != ymin)
        return xmin < ymin;
    return xmax < ymax;
}

static void op_aresta_min(void *in_, void *inout_, int *len, MPI_Datatype *dt)
{
    (void) dt;
    const ArestaMin *in = (const ArestaMin *) in_;
    ArestaMin *inout = (ArestaMin *) inout_;
    for (int i = 0; i < *len; i++)
        if (aresta_preferida(&in[i], &inout[i]))
            inout[i] = in[i];
}

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

    char lognome[64];
    snprintf(lognome, sizeof lognome, "log_rank%d.txt", rank);
    FILE *logf = fopen(lognome, "w");

    MPI_Barrier(MPI_COMM_WORLD);
    double t_total_ini = MPI_Wtime();

    uint32_t V = 0;
    uint64_t E_total = 0, E_local = 0, byte_ini = 0, byte_fim = 0;
    int binario = tem_sufixo(entrada, ".bin");

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

    MPI_Barrier(MPI_COMM_WORLD);
    double t_calc_ini = MPI_Wtime();

    DSU dsu;
    if (!dsu_init(&dsu, V)) {
        if (rank == 0) fprintf(stderr, "Erro: dsu_init falhou (V=%" PRIu32 ").\n", V);
        free(loc); if (logf) fclose(logf); MPI_Finalize(); return 1;
    }

    ArestaMin *best = (ArestaMin *) malloc((size_t) V * sizeof(ArestaMin));

    Aresta    *mst  = (Aresta *)    malloc((size_t) V * sizeof(Aresta));

    uint32_t  *comp = (uint32_t *)  malloc((size_t) V * sizeof(uint32_t));
    if (best == NULL || mst == NULL || comp == NULL) {
        if (rank == 0) fprintf(stderr, "Erro: sem memoria (best/mst/comp).\n");
        free(best); free(mst); free(comp); free(loc); dsu_free(&dsu);
        if (logf) fclose(logf);
        MPI_Finalize();
        return 1;
    }
    uint64_t mst_n = 0;
    double   soma  = 0.0;

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

    MPI_Op OP_ARESTA_MIN;
    MPI_Op_create(op_aresta_min, 1 , &OP_ARESTA_MIN);

    int      finalizado = 0;
    uint32_t fase = 0;
    int      d3_floresta = 0;

    while (!finalizado) {
        fase++;
        double tf0 = MPI_Wtime();

        for (uint32_t c = 0; c < V; c++) {
            best[c].valido = 0;
            best[c].peso = 0;
            best[c].u = 0;
            best[c].v = 0;
        }

        for (uint32_t i = 0; i < V; i++) comp[i] = dsu_find(&dsu, i);

        for (uint64_t i = 0; i < E_local; i++) {
            uint32_t ru = comp[loc[i].u];
            uint32_t rv = comp[loc[i].v];
            if (ru == rv)
                continue;

            ArestaMin cand = { 1 , loc[i].peso, loc[i].u, loc[i].v };
            if (aresta_preferida(&cand, &best[ru])) best[ru] = cand;
            if (aresta_preferida(&cand, &best[rv])) best[rv] = cand;
        }

        MPI_Allreduce(MPI_IN_PLACE, best, (int) V, TIPO_ARESTA_MIN,
                      OP_ARESTA_MIN, MPI_COMM_WORLD);

        uint64_t fusoes = 0;
        for (uint32_t c = 0; c < V; c++) {
            if (!best[c].valido)
                continue;
            uint32_t u = best[c].u, v = best[c].v;
            if (dsu_union(&dsu, u, v)) {
                mst[mst_n].u = u;
                mst[mst_n].v = v;
                mst[mst_n].peso = best[c].peso;
                mst_n++;
                soma += best[c].peso;
                fusoes++;
            }
        }

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

            finalizado = 1;
        }
    }

    uint32_t ncomp = dsu_num_componentes(&dsu);
    if (ncomp > 1)
        d3_floresta = 1;

    MPI_Barrier(MPI_COMM_WORLD);
    double t_fim   = MPI_Wtime();
    double t_total = t_fim - t_total_ini;
    double t_calc  = t_fim - t_calc_ini;

    if (logf) {
        fprintf(logf, "----- fim -----\n");
        fprintf(logf, "fases totais.......: %" PRIu32 "\n", fase);
        fprintf(logf, "arestas na floresta: %" PRIu64 "\n", mst_n);
        fprintf(logf, "componentes finais.: %" PRIu32 "\n", ncomp);
        fprintf(logf, "tempo TOTAL........: %.6f s\n", t_total);
        fprintf(logf, "tempo SO calculo...: %.6f s\n", t_calc);
        fclose(logf);
    }

    if (rank == 0) {

        printf("Peso total da MST: %.10g\n", soma);

        escrever_saida_rank0(saida, mst, mst_n);

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

    MPI_Op_free(&OP_ARESTA_MIN);
    MPI_Type_free(&TIPO_ARESTA_MIN);
    free(best);
    free(mst);
    free(comp);
    free(loc);
    dsu_free(&dsu);

    MPI_Finalize();
    return 0;
}

/* AGM por Boruvka SEQUENCIAL (Parte I). Dois modos:
 *  - em-RAM: grafo inteiro em vetor (texto ou .bin pequeno);
 *  - streaming: rele o .bin em blocos a cada fase (memoria O(bloco)).
 * Uso: ./sequencial [--stream] <arquivo_dados> [arquivo_saida]. */
#define _POSIX_C_SOURCE 200112L
#define _FILE_OFFSET_BITS 64

#include "grafo.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <inttypes.h>
#include <sys/types.h>

/* "Nenhuma" = sentinela de INDICE de aresta (nao de peso). Seguro pois E < UINT32_MAX. */
#define ARESTA_NENHUMA UINT32_MAX

static double seg(struct timespec a, struct timespec b)
{
    return (double) (b.tv_sec - a.tv_sec) + (double) (b.tv_nsec - a.tv_nsec) / 1e9;
}

/* True se 'nome' termina em 'suf' (case-insensitive): escolhe .bin vs texto. */
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

/* Desempate (ordem total): 'cand' preferida sobre 'atual' sse menor peso; em
 * empate, menor par (min,max) normalizado. Mesma regra do paralelo -> resultado
 * identico. 'atual' pode ser ARESTA_NENHUMA. Modo em-RAM (indices de aresta). */
static int for_preferido_sobre(const Aresta *arestas, uint32_t cand, uint32_t atual)
{
    if (atual == ARESTA_NENHUMA)
        return 1;

    double pc = arestas[cand].peso;
    double pa = arestas[atual].peso;
    if (pc < pa) return 1;
    if (pc > pa) return 0;

    uint32_t cu = arestas[cand].u, cv = arestas[cand].v;
    uint32_t au = arestas[atual].u, av = arestas[atual].v;
    uint32_t cmin = cu < cv ? cu : cv, cmax = cu < cv ? cv : cu;
    uint32_t amin = au < av ? au : av, amax = au < av ? av : au;

    if (cmin != amin) return cmin < amin;
    return cmax < amax;
}

/* Melhor aresta por componente no modo streaming (guarda a aresta explicita,
 * nao um indice, pois no streaming nao ha vetor global de arestas). */
typedef struct {
    uint8_t  valido;
    double   peso;
    uint32_t u;
    uint32_t v;
} MenorE;

/* Mesmo desempate de for_preferido_sobre, versao struct. */
static int menor_preferida(const MenorE *cand, const MenorE *atual)
{
    if (!cand->valido)  return 0;
    if (!atual->valido) return 1;
    if (cand->peso < atual->peso) return 1;
    if (cand->peso > atual->peso) return 0;
    uint32_t cmin = cand->u < cand->v ? cand->u : cand->v;
    uint32_t cmax = cand->u < cand->v ? cand->v : cand->u;
    uint32_t amin = atual->u < atual->v ? atual->u : atual->v;
    uint32_t amax = atual->u < atual->v ? atual->v : atual->u;
    if (cmin != amin) return cmin < amin;
    return cmax < amax;
}

/* Nº de blocos em que cada fase rele o .bin (12 -> ~1 GB/bloco no graph.bin). */
#define NUM_BLOCOS_STREAM 12

/* Considera a aresta (u,v) p/ os componentes ru e rv (ja sabidos distintos). */
static inline void atualiza_best(MenorE *best, uint32_t ru, uint32_t rv,
                                 uint32_t u, uint32_t v, double peso)
{
    MenorE cand = { 1, peso, u, v };
    if (menor_preferida(&cand, &best[ru])) best[ru] = cand;
    if (menor_preferida(&cand, &best[rv])) best[rv] = cand;
}

/* Nucleo Boruvka em STREAMING: a cada fase rele 'entrada' em blocos, acha a
 * melhor aresta de saida por componente, aplica as unioes. */
static int boruvka_stream(const char *entrada, uint32_t V, uint64_t E, DSU *dsu,
                          Aresta *mst, uint64_t *p_mst_n, double *p_soma,
                          uint32_t *p_fases, int *p_d3)
{
    uint64_t BUF = (E + NUM_BLOCOS_STREAM - 1) / NUM_BLOCOS_STREAM;
    if (BUF == 0) BUF = 1;

    fprintf(stderr,
            "Streaming..........: %d blocos/fase, ~%.0f MB por bloco "
            "(BUF=%" PRIu64 " arestas)\n",
            NUM_BLOCOS_STREAM,
            (double) BUF * sizeof(Aresta) / (1024.0 * 1024.0), BUF);

    MenorE   *best = (MenorE *)   malloc((size_t) V * sizeof(MenorE));
    uint32_t *comp = (uint32_t *) malloc((size_t) V * sizeof(uint32_t));
    Aresta   *buf  = (Aresta *)   malloc((size_t) BUF * sizeof(Aresta));
    if (best == NULL || comp == NULL || buf == NULL) {
        fprintf(stderr, "boruvka_stream: falha de alocacao (best/comp/buf).\n");
        free(best); free(comp); free(buf);
        return 0;
    }

    uint64_t cap = 0;
    {
        const char *env = getenv("SEQ_RAM_EDGES");
        if (env != NULL) cap = strtoull(env, NULL, 10);
    }
    Aresta  *ram    = NULL;
    uint64_t ram_n  = 0;
    int      in_ram = 0;
    if (cap > 0) {
        ram = (Aresta *) malloc((size_t) cap * sizeof(Aresta));
        if (ram == NULL) {
            fprintf(stderr, "  [captura] sem memoria para %" PRIu64
                    " arestas; seguindo em SO-LEITURA.\n", cap);
            cap = 0;
        } else {
            fprintf(stderr, "  [captura] habilitada: ate %" PRIu64
                    " arestas (~%.0f MB) podem migrar para RAM.\n",
                    cap, (double) cap * sizeof(Aresta) / (1024.0 * 1024.0));
        }
    }

    double   soma  = 0.0;
    uint64_t mst_n = 0;
    uint32_t fases = 0;
    int      d3    = 0;
    int concluido  = 0;

    while (!concluido) {
        fases++;

        /* comp[i] := raiz(i) 1x por fase (vetor plano; nenhuma uniao ocorre na
         * varredura, entao comp[] fica estavel). Depois zera o best. */
        for (uint32_t i = 0; i < V; i++) comp[i] = dsu_find(dsu, i);
        for (uint32_t c = 0; c < V; c++) best[c].valido = 0;

        if (in_ram) {
            /* arestas vivas ja em RAM: varre so elas e COMPACTA (remove as que
             * cairam na mesma componente) -> ram_n encolhe fase a fase. */
            uint64_t w = 0;
            for (uint64_t i = 0; i < ram_n; i++) {
                uint32_t cu = comp[ram[i].u], cv = comp[ram[i].v];
                if (cu == cv) continue;
                Aresta a = ram[i];
                ram[w++] = a; /* sobrevivente (w <= i) */
                atualiza_best(best, cu, cv, a.u, a.v, a.peso);
            }
            ram_n = w;
            fprintf(stderr, "  fase %2" PRIu32 ": %" PRIu64 " arestas em RAM\n",
                    fases, ram_n);
        } else {
            /* rele o .bin em blocos; opcionalmente captura as arestas vivas */
            FILE *fr = grafo_bin_abrir(entrada);
            if (fr == NULL) { free(best); free(comp); free(buf); free(ram); return 0; }

            int      capturando = (cap > 0);
            uint64_t cap_n      = 0;
            int      estourou   = 0;
            uint32_t blocos     = 0;
            uint64_t n;
            while ((n = grafo_bin_ler_bloco(fr, buf, BUF)) > 0) {
                for (uint64_t i = 0; i < n; i++) {
                    uint32_t cu = comp[buf[i].u], cv = comp[buf[i].v];
                    if (cu == cv) continue;
                    atualiza_best(best, cu, cv, buf[i].u, buf[i].v, buf[i].peso);
                    if (capturando) {
                        if (cap_n < cap) ram[cap_n++] = buf[i];
                        else { capturando = 0; estourou = 1; } /* nao coube: aborta captura */
                    }
                }
                blocos++;
            }
            fclose(fr);

            if (cap > 0 && !estourou) {
                /* todas as vivas couberam: fim das releituras de disco */
                ram_n  = cap_n;
                in_ram = 1;
                free(buf); buf = NULL;
                fprintf(stderr, "  fase %2" PRIu32 ": %" PRIu32 " blocos lidos; "
                        "%" PRIu64 " arestas vivas migradas para RAM (fim das releituras).\n",
                        fases, blocos, ram_n);
            } else {
                fprintf(stderr, "  fase %2" PRIu32 ": %" PRIu32 " blocos lidos%s\n",
                        fases, blocos, (cap > 0) ? " (ainda nao coube em RAM)" : "");
            }
        }

        /* parada: nenhum componente tem aresta de saida */
        int algum = 0;
        for (uint32_t c = 0; c < V; c++)
            if (best[c].valido) { algum = 1; break; }

        if (!algum) {
            concluido = 1;
        } else {
            /* aplica as unioes (dsu_union filtra ciclos/dupla escolha) */
            uint64_t fusoes = 0;
            for (uint32_t c = 0; c < V; c++) {
                if (!best[c].valido) continue;
                /* [SOMA PARCIAL] passa best[c].peso: o DSU acumula o peso do
                 * componente resultante (peso dos dois lados + esta aresta). */
                if (dsu_union(dsu, best[c].u, best[c].v, best[c].peso)) {
                    mst[mst_n].u = best[c].u;
                    mst[mst_n].v = best[c].v;
                    mst[mst_n].peso = best[c].peso;
                    mst_n++;
                    soma += best[c].peso;
                    fusoes++;
                }
            }
            if (fusoes == 0) { d3 = 1; concluido = 1; } /* guarda D3: sem progresso */
        }
    }

    free(best);
    free(comp);
    free(buf);
    free(ram);

    *p_mst_n = mst_n;
    *p_soma  = soma;
    *p_fases = fases;
    *p_d3    = d3;
    return 1;
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "uso: %s [--stream] <arquivo_dados> [arquivo_saida]\n", argv[0]);
        fprintf(stderr, "  --stream: le o .bin em BLOCOS (memoria O(bloco)); para\n"
                        "            arquivos grandes que nao cabem na RAM. Auto-ligado\n"
                        "            para .bin > 2 GB.\n");
        return 1;
    }

    /* flag opcional --stream/-s antes do arquivo */
    int stream_flag = 0;
    int ai = 1;
    if (argv[1][0] == '-') {
        if (strcmp(argv[1], "--stream") == 0 || strcmp(argv[1], "-s") == 0)
            stream_flag = 1;
        else
            fprintf(stderr, "Aviso: opcao '%s' desconhecida; ignorada.\n", argv[1]);
        ai = 2;
    }
    if (argc <= ai) {
        fprintf(stderr, "uso: %s [--stream] <arquivo_dados> [arquivo_saida]\n", argv[0]);
        return 1;
    }
    const char *entrada = argv[ai];
    const char *saida   = (argc >= ai + 2) ? argv[ai + 1] : "mst_sequencial.txt";

    struct timespec t0, t1;

    /* ===== 1) leitura ===== */
    uint32_t V = 0;
    uint64_t E = 0;
    int binario = tem_sufixo(entrada, ".bin");

    /* decide streaming: por flag, ou auto se o .bin passa de 2 GB */
    int stream = 0;
    if (binario) {
        if (stream_flag) {
            stream = 1;
        } else {
            FILE *ft = grafo_bin_abrir(entrada);
            if (ft) {
                if (fseeko(ft, 0, SEEK_END) == 0) {
                    off_t sz = ftello(ft);
                    if (sz > (off_t) 2 * 1024 * 1024 * 1024) stream = 1;
                }
                fclose(ft);
            }
        }
    } else if (stream_flag) {
        fprintf(stderr, "Aviso: --stream so vale para arquivos .bin; ignorado.\n");
    }

    clock_gettime(CLOCK_MONOTONIC, &t0);
    Aresta *arestas = NULL;
    if (stream) {
        /* streaming: nao carrega o grafo, so descobre V e E */
        if (!grafo_bin_info(entrada, &V, &E)) V = 0;
    } else {
        arestas = binario ? ler_grafo_binario(entrada, &V, &E)
                          : ler_grafo_texto(entrada, &V, &E);
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double tempo_io = seg(t0, t1);
    fprintf(stderr, "Leitor.............: %s\n",
            stream ? "binario STREAMING (.bin, blocos)"
                   : (binario ? "binario (.bin, em RAM)" : "texto"));

    if (V == 0 || (!stream && arestas == NULL)) {
        fprintf(stderr, "Erro: falha ao ler o grafo.\n");
        free(arestas);
        return 1;
    }

    /* ===== 2) calculo (Boruvka) ===== */
    clock_gettime(CLOCK_MONOTONIC, &t0);

    DSU dsu;
    if (!dsu_init(&dsu, V)) {
        fprintf(stderr, "Erro: falha ao inicializar DSU (V=%" PRIu32 ").\n", V);
        free(arestas);
        return 1;
    }

    /* mst[] = arestas da floresta resultante (<= V-1). */
    Aresta *mst = (Aresta *) malloc((size_t) V * sizeof(Aresta));
    if (mst == NULL) {
        fprintf(stderr, "Erro: falha de alocacao (mst).\n");
        dsu_free(&dsu); free(arestas);
        return 1;
    }

    double   soma   = 0.0;
    uint64_t mst_n  = 0;
    uint32_t fases  = 0;
    int desconexo_sem_fusao = 0;

    if (stream) {
        if (!boruvka_stream(entrada, V, E, &dsu, mst, &mst_n, &soma, &fases,
                            &desconexo_sem_fusao)) {
            fprintf(stderr, "Erro: falha no calculo em streaming.\n");
            free(mst); dsu_free(&dsu);
            return 1;
        }
    } else {
        /* ----- modo em-RAM: grafo inteiro em arestas[] ----- */
        /* comp[i]=raiz(i) por fase; menor[c]=indice da melhor aresta do comp. c. */
        uint32_t *comp  = (uint32_t *) malloc((size_t) V * sizeof(uint32_t));
        uint32_t *menor = (uint32_t *) malloc((size_t) V * sizeof(uint32_t));
        if (comp == NULL || menor == NULL) {
            fprintf(stderr, "Erro: falha de alocacao no calculo.\n");
            free(comp); free(menor); free(mst);
            dsu_free(&dsu); free(arestas);
            return 1;
        }

        uint64_t e_viv = E; /* arestas ainda vivas; compacta a cada fase (E fica p/ log) */
        int concluido = 0;
        while (!concluido) {
            fases++;

            for (uint32_t i = 0; i < V; i++)
                comp[i] = dsu_find(&dsu, i);

            for (uint32_t i = 0; i < V; i++)
                menor[i] = ARESTA_NENHUMA;
            uint64_t w = 0;
            for (uint64_t i = 0; i < e_viv; i++) {
                uint32_t cu = comp[arestas[i].u];
                uint32_t cv = comp[arestas[i].v];
                if (cu == cv)
                    continue;

                arestas[w] = arestas[i]; /* compacta in-place (w <= i) */
                if (for_preferido_sobre(arestas, (uint32_t) w, menor[cu]))
                    menor[cu] = (uint32_t) w;
                if (for_preferido_sobre(arestas, (uint32_t) w, menor[cv]))
                    menor[cv] = (uint32_t) w;
                w++;
            }
            e_viv = w; /* proximas fases so varrem as vivas */

            int algum = 0;
            for (uint32_t c = 0; c < V; c++) {
                if (menor[c] != ARESTA_NENHUMA) { algum = 1; break; }
            }

            if (!algum) {
                /* parada literal (Parte I): nenhuma arvore pode ser mesclada */
                concluido = 1;
            } else {
                /* adiciona a melhor aresta de cada componente; dsu_union evita ciclo */
                uint64_t fusoes = 0;
                for (uint32_t c = 0; c < V; c++) {
                    if (menor[c] == ARESTA_NENHUMA)
                        continue;
                    Aresta e = arestas[menor[c]];
                    /* [SOMA PARCIAL] passa e.peso: o DSU acumula o peso do
                     * componente resultante (peso dos dois lados + esta aresta). */
                    if (dsu_union(&dsu, e.u, e.v, e.peso)) {
                        mst[mst_n++] = e;
                        soma += e.peso;
                        fusoes++;
                    }
                }

                if (fusoes == 0) { /* guarda D3: sem progresso -> encerra */
                    desconexo_sem_fusao = 1;
                    concluido = 1;
                }
            }
        }

        free(comp);
        free(menor);
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double tempo_calc = seg(t0, t1);

    uint32_t n_comp = dsu_num_componentes(&dsu);

   
    double soma_parcial_dsu = 0.0;
    for (uint32_t i = 0; i < V; i++)
        if (dsu_find(&dsu, i) == i)
            soma_parcial_dsu += dsu.peso[i];

    /* ===== 3) saida ===== */
    printf("Peso total da MST: %.10g\n", soma);

    /* arquivo: MST como lista de arestas (podem ser milhoes -> nao no stdout) */
    FILE *fo = fopen(saida, "w");
    if (fo == NULL) {
        fprintf(stderr, "Aviso: nao foi possivel abrir '%s' para escrita.\n", saida);
    } else {
        for (uint64_t i = 0; i < mst_n; i++)
            fprintf(fo, "%" PRIu32 " → (%.10g) → %" PRIu32 "\n",
                    mst[i].u, mst[i].peso, mst[i].v);
        fclose(fo);
    }

    fprintf(stderr, "----------------------------------------\n");
    fprintf(stderr, "Vertices (V).......: %" PRIu32 "\n", V);
    fprintf(stderr, "Arestas  (E).......: %" PRIu64 "\n", E);
    fprintf(stderr, "Fases de Boruvka...: %" PRIu32 "\n", fases);
    fprintf(stderr, "Arestas na floresta: %" PRIu64 "\n", mst_n);
    fprintf(stderr, "Componentes finais.: %" PRIu32 "\n", n_comp);
    /* [SOMA PARCIAL] peso total lido do DSU + checagem de consistencia. */
    fprintf(stderr, "Soma parcial (DSU).: %.10g%s\n", soma_parcial_dsu,
            (soma_parcial_dsu == soma) ? " (confere com o soma)"
                                       : " (DIVERGE do soma!)");
    fprintf(stderr, "Arestas gravadas em: %s\n", saida);
    fprintf(stderr, "Tempo de I/O.......: %.6f s\n", tempo_io);
    fprintf(stderr, "Tempo de calculo...: %.6f s\n", tempo_calc);
    if (n_comp > 1) {
        /* D3: grafo desconexo -> retorna floresta, nao 1 arvore */
        fprintf(stderr,
            "[D3] Condicao literal '1 componente' NAO atendida: grafo DESCONEXO "
            "(%" PRIu32 " componentes). Retornada a FLORESTA geradora minima.\n", n_comp);
        if (desconexo_sem_fusao)
            fprintf(stderr,
                "[D3] (encerrado pela guarda: fase sem nenhuma fusao)\n");
    }
    fprintf(stderr, "----------------------------------------\n");

    free(mst);
    dsu_free(&dsu);
    free(arestas);
    return 0;
}

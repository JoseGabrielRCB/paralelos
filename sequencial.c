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
#define _POSIX_C_SOURCE 200112L   /* clock_gettime + fseeko/ftello */
#define _FILE_OFFSET_BITS 64      /* off_t de 64 bits (arquivos > 2 GB) */

#include "grafo.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <inttypes.h>
#include <sys/types.h>

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

/* Verdadeiro se 'nome' termina em 'suf' (case-insensitive). Usado para escolher
 * o leitor: arquivos ".bin" => binário (graph.bin); senão => texto. */
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

    double pc = arestas[cand].peso;
    double pa = arestas[atual].peso;
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

/* ------------------------------------------------------------------------- *
 * MODO STREAMING SÓ-LEITURA (forma LITERAL da Parte I, para arquivos binários
 * que não cabem na RAM, ex.: graph.bin de 12,8 GB).
 *
 * O pseudocódigo literal varre o E ORIGINAL a cada fase ("para cada aresta uv
 * em E ..."). Fazemos exatamente isso: a cada fase RE-LEMOS o arquivo de entrada
 * em BLOCOS, calculando a menor aresta por componente. A "menor aresta por
 * componente" é guardada como struct EXPLÍCITA (peso+u+v+valido), igual ao
 * paralelo.c — o desempate é o MESMO (peso; em empate, menor par (min,max)),
 * então o resultado é idêntico ao modo em-RAM e ao paralelo.
 *
 * IMPORTANTE: este modo NÃO grava nada em disco (sem arquivos de rascunho), há
 * apenas LEITURA. Assim o uso de memória é O(bloco) + O(V) e a page cache é
 * totalmente recuperável, sem acúmulo de dirty pages — crucial em VMs pequenas
 * (ex.: WSL2 com 3,5 GB), onde escritas pesadas inundariam o cache e travariam
 * a VM. O custo é reler o arquivo a cada fase (~E*fases bytes lidos), aceitável
 * por ser I/O sequencial puro.
 * ------------------------------------------------------------------------- */
typedef struct {
    uint8_t  valido; /* 0 = "Nenhuma" */
    double   peso;
    uint32_t u;
    uint32_t v;
} MenorE;

/* 'cand' é preferida sobre 'atual'? (mesma ordem total de for_preferido_sobre) */
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

/* Em quantos blocos o arquivo ORIGINAL é subdividido por fase (cada fase relê
 * o arquivo de entrada em NUM_BLOCOS_STREAM pedaços). Para o graph.bin
 * (800M arestas = 12,8 GB), 12 blocos => ~1 GB por bloco. */
#define NUM_BLOCOS_STREAM 12

/* Atualiza best[] com uma aresta candidata (u,v,peso) já sabidamente entre
 * componentes diferentes ru!=rv. */
static inline void atualiza_best(MenorE *best, uint32_t ru, uint32_t rv,
                                 uint32_t u, uint32_t v, double peso)
{
    MenorE cand = { 1, peso, u, v };
    if (menor_preferida(&cand, &best[ru])) best[ru] = cand;
    if (menor_preferida(&cand, &best[rv])) best[rv] = cand;
}

/* Núcleo Borůvka em streaming SÓ-LEITURA (forma literal da Parte I). A cada fase
 * relê o arquivo 'entrada' em blocos. Preenche mst[]/soma/fases; usa o 'dsu' já
 * inicializado. Retorna 1 em sucesso, 0 em falha de alocação/I-O. */
static int boruvka_stream(const char *entrada, uint32_t V, uint64_t E, DSU *dsu,
                          Aresta *mst, uint64_t *p_mst_n, double *p_soma,
                          uint32_t *p_fases, int *p_d3)
{
    /* Bloco de leitura = teto(E / 12) arestas (subdivide a fonte em 12 pedaços). */
    uint64_t BUF = (E + NUM_BLOCOS_STREAM - 1) / NUM_BLOCOS_STREAM;
    if (BUF == 0) BUF = 1;

    fprintf(stderr,
            "Streaming..........: %d blocos/fase, ~%.0f MB por bloco "
            "(BUF=%" PRIu64 " arestas)\n",
            NUM_BLOCOS_STREAM,
            (double) BUF * sizeof(Aresta) / (1024.0 * 1024.0), BUF);

    MenorE   *best = (MenorE *)   malloc((size_t) V * sizeof(MenorE));
    uint32_t *comp = (uint32_t *) malloc((size_t) V * sizeof(uint32_t)); /* raiz por fase */
    Aresta   *buf  = (Aresta *)   malloc((size_t) BUF * sizeof(Aresta)); /* leitura */
    if (best == NULL || comp == NULL || buf == NULL) {
        fprintf(stderr, "boruvka_stream: falha de alocacao (best/comp/buf).\n");
        free(best); free(comp); free(buf);
        return 0;
    }

    /* CAPTURA-EM-RAM (opcional, opt-in por SEQ_RAM_EDGES=<n_arestas>).
     * O modo SÓ-LEITURA relê o arquivo a cada fase => I/O ~ E * fases (dominante
     * para o graph.bin). Como o DSU só FUNDE componentes (nunca separa), uma
     * aresta cujos extremos já caíram na mesma componente NUNCA mais será útil.
     * Logo, assim que as arestas que AINDA cruzam componentes couberem na RAM,
     * guardamos esse subconjunto e paramos de reler o disco — as fases seguintes
     * varrem só ele. O resultado é IDÊNTICO (descartar aresta interna é sempre
     * seguro). Desligado por padrão para preservar a memória O(bloco) em VMs
     * pequenas (WSL2); ligue com, p.ex., SEQ_RAM_EDGES=80000000 (~1,3 GB). */
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

        /* comp[i] := raiz(i) UMA vez por fase: troca 2 dsu_find por aresta (busca
         * de ponteiros / cache-miss) por 2 leituras de vetor plano. Como nenhuma
         * união ocorre durante a varredura, comp[] é estável e o resultado é
         * idêntico ao de chamar dsu_find por aresta. */
        for (uint32_t i = 0; i < V; i++) comp[i] = dsu_find(dsu, i);
        for (uint32_t c = 0; c < V; c++) best[c].valido = 0;

        if (in_ram) {
            /* Arestas vivas já estão em RAM: varre só elas, sem tocar no disco. */
            for (uint64_t i = 0; i < ram_n; i++) {
                uint32_t cu = comp[ram[i].u], cv = comp[ram[i].v];
                if (cu == cv) continue;
                atualiza_best(best, cu, cv, ram[i].u, ram[i].v, ram[i].peso);
            }
            fprintf(stderr, "  fase %2" PRIu32 ": %" PRIu64 " arestas em RAM\n",
                    fases, ram_n);
        } else {
            /* (a)-(c) RE-LÊ o E original em blocos (apenas leitura, sem rascunho);
             * para cada aresta uv entre componentes diferentes, atualiza o menor
             * do componente de u e o de v. Opcionalmente captura essas arestas. */
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
                    if (cu == cv) continue; /* mesma componente: descartada */
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
                /* Todas as arestas que cruzam componentes couberam em 'ram':
                 * a partir daqui NÃO relemos mais o disco. */
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

        /* (d) todos "Nenhuma" -> concluído. */
        int algum = 0;
        for (uint32_t c = 0; c < V; c++)
            if (best[c].valido) { algum = 1; break; }

        if (!algum) {
            concluido = 1;
        } else {
            /* (e) aplica as uniões (mesma ordem c crescente — resultado idêntico). */
            uint64_t fusoes = 0;
            for (uint32_t c = 0; c < V; c++) {
                if (!best[c].valido) continue;
                if (dsu_union(dsu, best[c].u, best[c].v)) {
                    mst[mst_n].u = best[c].u;
                    mst[mst_n].v = best[c].v;
                    mst[mst_n].peso = best[c].peso;
                    mst_n++;
                    soma += best[c].peso;
                    fusoes++;
                }
            }
            if (fusoes == 0) { d3 = 1; concluido = 1; }
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

    /* Flag opcional --stream (ou -s) antes do arquivo. */
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

    /* ===================== 1) LEITURA (I/O) ===================== */
    uint32_t V = 0;
    uint64_t E = 0;
    int binario = tem_sufixo(entrada, ".bin");

    /* Decide STREAMING: por flag, ou auto se o .bin passa de 2 GB (nao caberia
     * na RAM como vetor unico). Streaming exige binario. */
    int stream = 0;
    if (binario) {
        if (stream_flag) {
            stream = 1;
        } else {
            FILE *ft = grafo_bin_abrir(entrada);
            if (ft) {
                if (fseeko(ft, 0, SEEK_END) == 0) {
                    off_t sz = ftello(ft);
                    if (sz > (off_t) 2 * 1024 * 1024 * 1024) stream = 1; /* > 2 GB */
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
        /* Streaming: nao carrega o grafo; so descobre V e E num passe. */
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

    /* ===================== 2) CÁLCULO (Borůvka) ===================== */
    clock_gettime(CLOCK_MONOTONIC, &t0);

    DSU dsu;
    if (!dsu_init(&dsu, V)) {
        fprintf(stderr, "Erro: falha ao inicializar DSU (V=%" PRIu32 ").\n", V);
        free(arestas);
        return 1;
    }

    /* mst[] = arestas adicionadas a E' (a floresta resultante); <= V-1 arestas. */
    Aresta *mst = (Aresta *) malloc((size_t) V * sizeof(Aresta));
    if (mst == NULL) {
        fprintf(stderr, "Erro: falha de alocacao (mst).\n");
        dsu_free(&dsu); free(arestas);
        return 1;
    }

    double   soma   = 0.0; /* D2: soma dos pesos em ponto flutuante */
    uint64_t mst_n  = 0;  /* número de arestas em E' */
    uint32_t fases  = 0;
    int desconexo_sem_fusao = 0; /* sinaliza a guarda D3 (fase sem fusão) */

    if (stream) {
        /* ----- MODO STREAMING: re-lê o arquivo em blocos a cada fase. ----- */
        if (!boruvka_stream(entrada, V, E, &dsu, mst, &mst_n, &soma, &fases,
                            &desconexo_sem_fusao)) {
            fprintf(stderr, "Erro: falha no calculo em streaming.\n");
            free(mst); dsu_free(&dsu);
            return 1;
        }
    } else {
        /* ----- MODO EM-RAM: grafo inteiro em 'arestas[]'. ----- */
        /* comp[i]  = componente (raiz) do vértice i, recalculado a cada fase.
         * menor[c] = índice da aresta de menor peso do componente c ("Nenhuma"
         *            se ARESTA_NENHUMA). */
        uint32_t *comp  = (uint32_t *) malloc((size_t) V * sizeof(uint32_t));
        uint32_t *menor = (uint32_t *) malloc((size_t) V * sizeof(uint32_t));
        if (comp == NULL || menor == NULL) {
            fprintf(stderr, "Erro: falha de alocacao no calculo.\n");
            free(comp); free(menor); free(mst);
            dsu_free(&dsu); free(arestas);
            return 1;
        }

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
                        soma += e.peso; /* D2 */
                        fusoes++;
                    }
                }

                /* GUARDA D3: se havia arestas candidatas mas nenhuma fusão
                 * ocorreu, não há progresso possível — encerra retornando a
                 * floresta. (Com DSU válido isto não deveria ocorrer.) */
                if (fusoes == 0) {
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

    /* ===================== 3) SAÍDA ===================== */

    /* 3.1 stdout: peso total (ponto flutuante). */
    printf("Peso total da MST: %.10g\n", soma);

    /* 3.2 arquivo: MST como LISTA DE ARESTAS, uma por linha, no formato do
     *     enunciado "<u> → (<peso>) → <v>" (mesmo formato do paralelo.c, para
     *     saídas idênticas). Podem ser milhões -> NÃO vão ao stdout. */
    FILE *fo = fopen(saida, "w");
    if (fo == NULL) {
        fprintf(stderr, "Aviso: nao foi possivel abrir '%s' para escrita.\n", saida);
    } else {
        for (uint64_t i = 0; i < mst_n; i++)
            fprintf(fo, "%" PRIu32 " → (%.10g) → %" PRIu32 "\n",
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
    /* comp/menor já foram liberados no modo em-RAM; arestas é NULL em streaming. */
    free(mst);
    dsu_free(&dsu);
    free(arestas);
    return 0;
}

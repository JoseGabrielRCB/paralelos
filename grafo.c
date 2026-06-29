/*
 * grafo.c — Implementação do DSU (Union-Find) e do I/O (texto e binário).
 */
/* Para fseeko/ftello e off_t de 64 bits (arquivos > 2 GB, ex.: graph.bin). */
#define _FILE_OFFSET_BITS 64
#define _POSIX_C_SOURCE 200112L

#include "grafo.h"
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <inttypes.h>
#include <sys/types.h>

/* ------------------------------------------------------------------------- *
 * DSU
 * ------------------------------------------------------------------------- */

int dsu_init(DSU *dsu, uint32_t n)
{
    if (dsu == NULL || n == 0)
        return 0;

    dsu->pai  = (uint32_t *) malloc((size_t) n * sizeof(uint32_t));
    dsu->rank = (uint32_t *) calloc((size_t) n, sizeof(uint32_t)); /* rank=0 */
    if (dsu->pai == NULL || dsu->rank == NULL) {
        /* Falha de alocação: libera o que tiver sido alocado e zera. */
        free(dsu->pai);
        free(dsu->rank);
        dsu->pai = NULL;
        dsu->rank = NULL;
        dsu->n = 0;
        return 0;
    }

    /* Cada vértice começa como seu próprio conjunto (singleton). */
    for (uint32_t i = 0; i < n; i++)
        dsu->pai[i] = i;

    dsu->n = n;
    return 1;
}

void dsu_free(DSU *dsu)
{
    if (dsu == NULL)
        return;
    free(dsu->pai);
    free(dsu->rank);
    dsu->pai = NULL;
    dsu->rank = NULL;
    dsu->n = 0;
}

uint32_t dsu_find(DSU *dsu, uint32_t x)
{
    /* COMPRESSÃO DE CAMINHO: faz cada nó visitado apontar direto p/ a raiz.
     * Implementação iterativa (segura para V grande — evita estouro de pilha). */
    uint32_t raiz = x;
    while (dsu->pai[raiz] != raiz)
        raiz = dsu->pai[raiz];

    /* Segunda passada: aponta todos do caminho para a raiz encontrada. */
    while (dsu->pai[x] != raiz) {
        uint32_t prox = dsu->pai[x];
        dsu->pai[x] = raiz;
        x = prox;
    }
    return raiz;
}

int dsu_union(DSU *dsu, uint32_t a, uint32_t b)
{
    uint32_t ra = dsu_find(dsu, a);
    uint32_t rb = dsu_find(dsu, b);

    if (ra == rb)
        return 0; /* já no mesmo conjunto: fusão redundante */

    /* UNION BY RANK: pendura a árvore de menor rank sob a de maior rank. */
    if (dsu->rank[ra] < dsu->rank[rb]) {
        dsu->pai[ra] = rb;
    } else if (dsu->rank[ra] > dsu->rank[rb]) {
        dsu->pai[rb] = ra;
    } else {
        /* Ranks iguais: escolhe uma raiz e incrementa seu rank. */
        dsu->pai[rb] = ra;
        dsu->rank[ra]++;
    }
    return 1; /* houve fusão */
}

uint32_t dsu_num_componentes(DSU *dsu)
{
    /* Um componente por vértice que é raiz de si mesmo. */
    uint32_t comp = 0;
    for (uint32_t i = 0; i < dsu->n; i++)
        if (dsu_find(dsu, i) == i)
            comp++;
    return comp;
}

/* ------------------------------------------------------------------------- *
 * I/O — STUBS (Fase 0). Implementação real nas fases seguintes.
 * ------------------------------------------------------------------------- */

/* Retorna a próxima linha NÃO vazia (que tenha algum caractere não-branco) de f,
 * gravada em buf. Retorna NULL no fim de arquivo. Usada só para o cabeçalho. */
static char *proxima_linha_nao_vazia(FILE *f, char *buf, int n)
{
    while (fgets(buf, n, f) != NULL) {
        for (const char *p = buf; *p; p++)
            if (!isspace((unsigned char) *p))
                return buf; /* achou conteúdo */
        /* linha só com espaços/quebra: ignora e tenta a próxima */
    }
    return NULL; /* EOF */
}

Aresta *ler_grafo_texto(const char *caminho, uint32_t *V, uint64_t *E)
{
    FILE *f = fopen(caminho, "r");
    if (f == NULL) {
        fprintf(stderr, "ler_grafo_texto: nao foi possivel abrir '%s'.\n", caminho);
        return NULL;
    }

    char linha[256];

    /* --- Cabeçalho: linha 1 = V, linha 2 = E (confiar no E do cabeçalho). --- */
    if (proxima_linha_nao_vazia(f, linha, sizeof linha) == NULL ||
        sscanf(linha, "%" SCNu32, V) != 1) {
        fprintf(stderr, "ler_grafo_texto: cabecalho invalido (V).\n");
        fclose(f);
        return NULL;
    }
    uint64_t E_cab = 0;
    if (proxima_linha_nao_vazia(f, linha, sizeof linha) == NULL ||
        sscanf(linha, "%" SCNu64, &E_cab) != 1) {
        fprintf(stderr, "ler_grafo_texto: cabecalho invalido (E).\n");
        fclose(f);
        return NULL;
    }

    /* Aloca exatamente E arestas (confiando no cabeçalho). */
    Aresta *arestas = (Aresta *) malloc((size_t) E_cab * sizeof(Aresta));
    if (arestas == NULL) {
        fprintf(stderr, "ler_grafo_texto: falha ao alocar %" PRIu64 " arestas.\n", E_cab);
        fclose(f);
        return NULL;
    }

    /* --- Arestas: "u v peso". Linhas em branco/mal-formadas sao IGNORADAS
     *     (o arquivo tem 1 linha vazia no final). Para de ler ao atingir E. --- */
    uint64_t lidas = 0;
    uint64_t ignoradas = 0;
    while (lidas < E_cab && fgets(linha, sizeof linha, f) != NULL) {
        uint32_t u, v;
        double peso; /* D2: peso em ponto flutuante */
        if (sscanf(linha, "%" SCNu32 " %" SCNu32 " %lf", &u, &v, &peso) != 3) {
            ignoradas++; /* linha em branco ou mal-formada */
            continue;
        }
        /* Índices base 0: validos em 0..V-1. Fora disso = mal-formada. */
        if (u >= *V || v >= *V) {
            ignoradas++;
            continue;
        }
        arestas[lidas].u = u;
        arestas[lidas].v = v;
        arestas[lidas].peso = peso;
        lidas++;
    }
    fclose(f);

    if (ignoradas > 0)
        fprintf(stderr, "ler_grafo_texto: %" PRIu64 " linha(s) ignorada(s) (vazias/mal-formadas).\n",
                ignoradas);
    if (lidas < E_cab)
        fprintf(stderr, "ler_grafo_texto: AVISO — lidas %" PRIu64 " de %" PRIu64 " arestas do cabecalho.\n",
                lidas, E_cab);

    *E = lidas; /* número de arestas efetivamente carregadas */
    return arestas;
}

/* ------------------------------------------------------------------------- *
 * ler_grafo_binario — leitor BINÁRIO sequencial do formato de `graph.bin`.
 *
 * Cada registro tem 16 bytes (u uint32, v uint32, peso double), exatamente o
 * layout de struct Aresta, então a leitura é uma cópia direta para o vetor.
 * Não há cabeçalho: E = tamanho/16 e V = (maior índice de vértice) + 1.
 * ------------------------------------------------------------------------- */
Aresta *ler_grafo_binario(const char *caminho, uint32_t *V, uint64_t *E)
{
    FILE *f = fopen(caminho, "rb");
    if (f == NULL) {
        fprintf(stderr, "ler_grafo_binario: nao foi possivel abrir '%s'.\n", caminho);
        return NULL;
    }

    /* Tamanho do arquivo (fseek/ftello para suportar > 2 GB). */
    if (fseeko(f, 0, SEEK_END) != 0) {
        fprintf(stderr, "ler_grafo_binario: fseeko falhou.\n");
        fclose(f);
        return NULL;
    }
    off_t fsize = ftello(f);
    rewind(f);
    if (fsize < 0 || (uint64_t) fsize % sizeof(Aresta) != 0) {
        fprintf(stderr, "ler_grafo_binario: tamanho (%lld) nao e multiplo de %zu.\n",
                (long long) fsize, sizeof(Aresta));
        fclose(f);
        return NULL;
    }
    uint64_t n = (uint64_t) fsize / sizeof(Aresta);

    Aresta *arestas = (Aresta *) malloc((size_t) n * sizeof(Aresta));
    if (arestas == NULL) {
        fprintf(stderr, "ler_grafo_binario: falha ao alocar %" PRIu64 " arestas (%.1f GB).\n",
                n, (double) n * sizeof(Aresta) / 1e9);
        fclose(f);
        return NULL;
    }

    /* Leitura em blocos (fread de >2 GB de uma vez pode falhar em algumas libc). */
    const uint64_t BLOCO = 16ULL * 1024 * 1024; /* 16M arestas por bloco */
    uint64_t lidas = 0;
    while (lidas < n) {
        uint64_t querer = (n - lidas < BLOCO) ? (n - lidas) : BLOCO;
        size_t got = fread(arestas + lidas, sizeof(Aresta), (size_t) querer, f);
        if (got == 0) break;
        lidas += got;
    }
    fclose(f);

    if (lidas != n) {
        fprintf(stderr, "ler_grafo_binario: AVISO — lidas %" PRIu64 " de %" PRIu64 " arestas.\n",
                lidas, n);
        n = lidas;
    }

    /* V = maior índice de vértice + 1 (não há cabeçalho). */
    uint32_t vmax = 0;
    for (uint64_t i = 0; i < n; i++) {
        if (arestas[i].u > vmax) vmax = arestas[i].u;
        if (arestas[i].v > vmax) vmax = arestas[i].v;
    }

    *V = (n > 0) ? (vmax + 1) : 0;
    *E = n;
    return arestas;
}

/* ------------------------------------------------------------------------- *
 * Leitura BINÁRIA em STREAMING — ver descrição em grafo.h.
 * ------------------------------------------------------------------------- */

FILE *grafo_bin_abrir(const char *caminho)
{
    FILE *f = fopen(caminho, "rb");
    if (f == NULL)
        fprintf(stderr, "grafo_bin_abrir: nao foi possivel abrir '%s'.\n", caminho);
    return f;
}

uint64_t grafo_bin_ler_bloco(FILE *f, Aresta *buf, uint64_t cap)
{
    if (f == NULL || buf == NULL || cap == 0)
        return 0;
    /* Cada Aresta tem 16 bytes idênticos ao registro em disco: fread direto. */
    return (uint64_t) fread(buf, sizeof(Aresta), (size_t) cap, f);
}

int grafo_bin_info(const char *caminho, uint32_t *V, uint64_t *E)
{
    FILE *f = grafo_bin_abrir(caminho);
    if (f == NULL)
        return 0;

    /* Buffer temporário de 1M arestas (16 MB) para o passe de varredura. */
    const uint64_t CAP = 1024ULL * 1024;
    Aresta *buf = (Aresta *) malloc((size_t) CAP * sizeof(Aresta));
    if (buf == NULL) {
        fprintf(stderr, "grafo_bin_info: falha ao alocar buffer.\n");
        fclose(f);
        return 0;
    }

    uint64_t total = 0;
    uint32_t vmax = 0;
    uint64_t n;
    while ((n = grafo_bin_ler_bloco(f, buf, CAP)) > 0) {
        for (uint64_t i = 0; i < n; i++) {
            if (buf[i].u > vmax) vmax = buf[i].u;
            if (buf[i].v > vmax) vmax = buf[i].v;
        }
        total += n;
    }

    free(buf);
    fclose(f);

    *E = total;
    *V = (total > 0) ? (vmax + 1) : 0;
    return 1;
}

#ifdef HAVE_MPI
/* ===========================================================================
 * Implementação MPI-IO (compilada apenas com mpicc -DHAVE_MPI).
 *
 * Desvio "Estilo B / MPI-IO" autorizado: a leitura/distribuição é feita por
 * MPI-IO; cada rank lê só a SUA fatia de bytes do arquivo (escalável, sem o
 * rank 0 ler tudo e espalhar). D1: lista de arestas (não a matriz W(u,v)).
 * =========================================================================== */
#include <mpi.h>

Aresta *ler_grafo_mpiio(const char *caminho, uint32_t *V, uint64_t *E_total,
                        int rank, int nprocs, uint64_t *E_local,
                        uint64_t *byte_ini, uint64_t *byte_fim)
{
    /* Comprimento máximo de uma linha "u v peso" (3 uint32 = até ~33 chars);
     * usamos folga. Serve de OVERLAP para completar a linha que cruza a borda. */
    const uint64_t MAX_LINHA = 128;

    /* ----------------------------------------------------------------------
     * 1) CABEÇALHO: só o rank 0 lê V, E e o offset onde começam as arestas.
     *    Em seguida difunde os três valores a todos (MPI_Bcast).
     * ---------------------------------------------------------------------- */
    uint32_t Vh = 0;
    uint64_t Eh = 0;
    uint64_t data_start = 0; /* byte de início da 1ª aresta (começo da linha 3) */

    if (rank == 0) {
        FILE *f = fopen(caminho, "r");
        if (f != NULL) {
            char linha[256];
            /* lê as 2 primeiras linhas NÃO vazias: V e E */
            if (proxima_linha_nao_vazia(f, linha, sizeof linha))
                sscanf(linha, "%" SCNu32, &Vh);
            if (proxima_linha_nao_vazia(f, linha, sizeof linha))
                sscanf(linha, "%" SCNu64, &Eh);
            /* Após ler as 2 linhas, a posição do arquivo é exatamente o início
             * da linha 3 (1ª aresta). ftell dá esse offset em bytes.
             * // VERIFICAR: assume quebras de linha '\n' (sem CRLF). O arquivo
             *    de dados verificado usa LF; com CRLF o offset ainda seria
             *    correto, pois ftell conta os bytes realmente consumidos. */
            long pos = ftell(f);
            data_start = (pos >= 0) ? (uint64_t) pos : 0;
            fclose(f);
        } else {
            fprintf(stderr, "ler_grafo_mpiio: rank 0 nao abriu '%s'.\n", caminho);
        }
    }

    /* Difusão coletiva do cabeçalho (todos os ranks precisam de V e do offset). */
    MPI_Bcast(&Vh,         1, MPI_UINT32_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Eh,         1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(&data_start, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    if (Vh == 0)
        return NULL; /* cabeçalho inválido em todos os ranks */

    *V = Vh;
    *E_total = Eh;

    /* ----------------------------------------------------------------------
     * 2) ABERTURA COLETIVA via MPI-IO e tamanho total do arquivo.
     * ---------------------------------------------------------------------- */
    MPI_File fh;
    if (MPI_File_open(MPI_COMM_WORLD, caminho, MPI_MODE_RDONLY,
                      MPI_INFO_NULL, &fh) != MPI_SUCCESS) {
        if (rank == 0)
            fprintf(stderr, "ler_grafo_mpiio: MPI_File_open falhou.\n");
        return NULL;
    }
    MPI_Offset fsize_off = 0;
    MPI_File_get_size(fh, &fsize_off);
    uint64_t fsize = (uint64_t) fsize_off;

    /* ----------------------------------------------------------------------
     * 3) PARTIÇÃO em 'nprocs' blocos contíguos da região de arestas
     *    [data_start, fsize). Fronteiras nominais (tiling exato):
     *        b(i) = data_start + (L * i) / nprocs
     *    O rank i é dono das LINHAS cujo INÍCIO (offset do 1º byte) cai em
     *    [b(i), b(i+1)). Assim cada linha pertence a exatamente um rank.
     * ---------------------------------------------------------------------- */
    uint64_t L   = (fsize > data_start) ? (fsize - data_start) : 0;
    uint64_t bi  = data_start + (L * (uint64_t) rank)       / (uint64_t) nprocs;
    uint64_t bi1 = data_start + (L * (uint64_t) (rank + 1)) / (uint64_t) nprocs;

    /* ALINHAMENTO DE BORDAS:
     *  - Borda inicial: a linha que cruza b(i) pertence ao rank i-1. Para
     *    decidir corretamente (inclusive quando b(i) cai EXATAMENTE no início
     *    de uma linha), lemos 1 byte de "prefixo" antes de b(i): se esse byte
     *    for '\n', então b(i) já é início de linha (não há nada a descartar);
     *    senão, há uma linha parcial a descartar até o primeiro '\n'.
     *  - Borda final: lemos MAX_LINHA bytes ALÉM de b(i+1) para conseguir
     *    completar a última linha que começa antes de b(i+1) e cruza a borda.
     * O rank 0 começa em data_start, que já é início de linha (linha 3). */
    int tem_prefixo = (rank != 0 && bi > data_start);
    uint64_t read_lo = tem_prefixo ? (bi - 1) : bi;
    uint64_t read_hi = bi1 + MAX_LINHA;
    if (read_hi > fsize) read_hi = fsize;
    uint64_t nbytes = (read_hi > read_lo) ? (read_hi - read_lo) : 0;

    *byte_ini = read_lo;
    *byte_fim = read_hi;

    /* Buffer da fatia (+1 para terminador nulo).
     * LIMITAÇÃO (nbytes em int): o parâmetro 'count' de MPI_File_read_at_all é
     * um int, logo a fatia de um rank deve caber em ~2 GB. Para o dataset deste
     * trabalho (~70 MB) e qualquer nº razoável de ranks isso sempre vale. Para
     * fatias > 2 GB seria preciso usar MPI_Count (MPI_File_read_at_all_c, MPI-4)
     * ou ler em vários blocos de até INT_MAX bytes — não exigido nesta escala. */
    char *buf = (char *) malloc((size_t) nbytes + 1);
    if (buf == NULL) {
        /* LIMITAÇÃO (alocação antes de coletiva): se UM rank falhasse aqui e
         * retornasse antes do MPI_File_read_at_all (coletivo), os demais ficariam
         * bloqueados. Nesta escala (fatia ~ dezenas de MB) assume-se que o malloc
         * sempre tem sucesso; um tratamento robusto exigiria acordar a falha
         * entre todos os ranks (ex.: MPI_Allreduce do código de erro) antes da
         * coletiva. */
        MPI_File_close(&fh);
        return NULL;
    }

    /* ----------------------------------------------------------------------
     * 4) LEITURA COLETIVA da fatia de bytes. Todos os ranks participam
     *    (read_at_all é coletiva); offsets/contagens podem diferir entre ranks.
     * ---------------------------------------------------------------------- */
    MPI_Status st;
    MPI_File_read_at_all(fh, (MPI_Offset) read_lo, buf, (int) nbytes,
                         MPI_CHAR, &st);
    buf[nbytes] = '\0';
    MPI_File_close(&fh);

    /* ----------------------------------------------------------------------
     * 5) PARSING ALINHADO: percorre as linhas cujo INÍCIO cai em [b(i), b(i+1)).
     * ---------------------------------------------------------------------- */
    /* Define o índice (no buffer) da primeira linha que pertence a este rank. */
    uint64_t start_idx;
    if (rank == 0) {
        start_idx = 0; /* data_start já é início de linha */
    } else if (tem_prefixo) {
        /* buf[0] é o byte em (b(i)-1). */
        if (buf[0] == '\n') {
            start_idx = 1; /* b(i) é exatamente um início de linha */
        } else {
            /* há linha parcial (pertence ao rank i-1): pula até após o 1º '\n' */
            uint64_t k = 1;
            while (k < nbytes && buf[k] != '\n') k++;
            start_idx = (k < nbytes) ? (k + 1) : nbytes;
        }
    } else {
        /* caso degenerado (bi == data_start em rank>0): trata como rank 0. */
        start_idx = 0;
    }

    /* Vetor local de arestas com capacidade crescente. */
    uint64_t cap = (nbytes / 16) + 16; /* estimativa: ~16 bytes por linha */
    Aresta *loc = (Aresta *) malloc((size_t) cap * sizeof(Aresta));
    if (loc == NULL) {
        free(buf);
        *E_local = 0;
        return NULL;
    }
    uint64_t n = 0;
    uint64_t ignoradas = 0;

    uint64_t idx = start_idx;
    while (idx < nbytes) {
        uint64_t abs_inicio = read_lo + idx; /* offset absoluto do início da linha */
        if (abs_inicio >= bi1)
            break; /* esta linha já pertence ao próximo rank */

        /* Delimita a linha [idx, j): vai até '\n' ou fim do buffer. A linha que
         * cruza b(i+1) cabe inteira graças ao OVERLAP (MAX_LINHA). */
        uint64_t j = idx;
        while (j < nbytes && buf[j] != '\n') j++;
        if (j < nbytes) buf[j] = '\0'; /* termina a string da linha */

        uint32_t u, v, w;
        if (sscanf(buf + idx, "%" SCNu32 " %" SCNu32 " %" SCNu32, &u, &v, &w) == 3
            && u < Vh && v < Vh) {
            if (n == cap) {
                cap *= 2;
                Aresta *novo = (Aresta *) realloc(loc, (size_t) cap * sizeof(Aresta));
                if (novo == NULL) { free(loc); free(buf); *E_local = 0; return NULL; }
                loc = novo;
            }
            loc[n].u = u;
            loc[n].v = v;
            loc[n].peso = w;
            n++;
        } else {
            ignoradas++; /* linha em branco/mal-formada (ex.: a linha vazia final) */
        }

        idx = j + 1; /* avança para a próxima linha (após o '\n') */
    }

    free(buf);

    fprintf(stderr, "[rank %d] arestas locais=%" PRIu64
            " ignoradas=%" PRIu64 " bytes=[%" PRIu64 ",%" PRIu64 ")\n",
            rank, n, ignoradas, read_lo, read_hi);

    *E_local = n;
    return loc;
}

/* ===========================================================================
 * ler_grafo_binario_mpiio — leitor BINÁRIO via MPI-IO (Estilo B).
 *
 * Como cada registro tem 16 bytes FIXOS, a partição é aritmética exata e NÃO
 * exige alinhamento de bordas de linha (ao contrário do leitor texto). O rank i
 * lê os registros [i*E/np, (i+1)*E/np) com MPI_File_read_at_all e reinterpreta
 * os bytes diretamente como vetor de Aresta (mesmo layout de 16 bytes).
 *
 * Sem cabeçalho: E_total = tamanho/16; V = MPI_Allreduce(MAX) do maior índice
 * de vértice local + 1.
 * =========================================================================== */
Aresta *ler_grafo_binario_mpiio(const char *caminho, uint32_t *V, uint64_t *E_total,
                                int rank, int nprocs, uint64_t *E_local,
                                uint64_t *rec_ini, uint64_t *rec_fim)
{
    const uint64_t REC = sizeof(Aresta); /* 16 bytes por aresta */

    /* Tipo MPI de um registro: 16 bytes contíguos. Usar um tipo de 16 bytes
     * (em vez de MPI_CHAR com count em bytes) mantém o 'count' pequeno e evita
     * estourar o int de count quando a fatia passa de 2 GB. */
    MPI_Datatype TIPO_REC;
    MPI_Type_contiguous((int) REC, MPI_BYTE, &TIPO_REC);
    MPI_Type_commit(&TIPO_REC);

    /* Abertura coletiva e tamanho do arquivo. */
    MPI_File fh;
    if (MPI_File_open(MPI_COMM_WORLD, caminho, MPI_MODE_RDONLY,
                      MPI_INFO_NULL, &fh) != MPI_SUCCESS) {
        if (rank == 0)
            fprintf(stderr, "ler_grafo_binario_mpiio: MPI_File_open falhou ('%s').\n", caminho);
        MPI_Type_free(&TIPO_REC);
        return NULL;
    }
    MPI_Offset fsize_off = 0;
    MPI_File_get_size(fh, &fsize_off);
    uint64_t fsize = (uint64_t) fsize_off;

    if (fsize % REC != 0 && rank == 0)
        fprintf(stderr, "ler_grafo_binario_mpiio: AVISO — tamanho %" PRIu64
                " nao e multiplo de %" PRIu64 " (registros truncados ignorados).\n",
                fsize, REC);
    uint64_t E = fsize / REC;
    *E_total = E;

    /* Partição em registros (divisão exata por offset). */
    uint64_t r_ini = (E * (uint64_t) rank)       / (uint64_t) nprocs;
    uint64_t r_fim = (E * (uint64_t) (rank + 1)) / (uint64_t) nprocs;
    uint64_t nrec  = (r_fim > r_ini) ? (r_fim - r_ini) : 0;
    *rec_ini = r_ini;
    *rec_fim = r_fim;

    /* Buffer local = vetor de Aresta (layout idêntico ao disco). */
    Aresta *loc = (Aresta *) malloc((size_t) nrec * sizeof(Aresta));
    if (loc == NULL && nrec > 0) {
        fprintf(stderr, "[rank %d] ler_grafo_binario_mpiio: malloc de %" PRIu64
                " arestas falhou.\n", rank, nrec);
        MPI_File_close(&fh);
        MPI_Type_free(&TIPO_REC);
        *E_local = 0;
        return NULL;
    }

    /* Leitura coletiva da fatia: offset em BYTES, count em REGISTROS.
     * 'nrec' <= E (<= ~800M) cabe folgado em int. */
    MPI_Status stt;
    MPI_File_read_at_all(fh, (MPI_Offset) (r_ini * REC), loc, (int) nrec,
                         TIPO_REC, &stt);
    MPI_File_close(&fh);
    MPI_Type_free(&TIPO_REC);

    /* V global = maior índice de vértice + 1 (Allreduce MAX). */
    uint32_t vmax_local = 0;
    for (uint64_t i = 0; i < nrec; i++) {
        if (loc[i].u > vmax_local) vmax_local = loc[i].u;
        if (loc[i].v > vmax_local) vmax_local = loc[i].v;
    }
    uint32_t vmax_global = 0;
    MPI_Allreduce(&vmax_local, &vmax_global, 1, MPI_UINT32_T, MPI_MAX, MPI_COMM_WORLD);
    *V = (E > 0) ? (vmax_global + 1) : 0;

    fprintf(stderr, "[rank %d] arestas locais=%" PRIu64
            " registros=[%" PRIu64 ",%" PRIu64 ") V=%" PRIu32 "\n",
            rank, nrec, r_ini, r_fim, *V);

    *E_local = nrec;
    return loc;
}

#else /* !HAVE_MPI — stub usado nos alvos compilados com gcc (teste/sequencial) */

Aresta *ler_grafo_mpiio(const char *caminho, uint32_t *V, uint64_t *E_total,
                        int rank, int nprocs, uint64_t *E_local,
                        uint64_t *byte_ini, uint64_t *byte_fim)
{
    (void) caminho; (void) V; (void) E_total; (void) rank; (void) nprocs;
    (void) E_local; (void) byte_ini; (void) byte_fim;
    fprintf(stderr, "ler_grafo_mpiio: compilado SEM -DHAVE_MPI (use mpicc).\n");
    return NULL;
}

Aresta *ler_grafo_binario_mpiio(const char *caminho, uint32_t *V, uint64_t *E_total,
                                int rank, int nprocs, uint64_t *E_local,
                                uint64_t *rec_ini, uint64_t *rec_fim)
{
    (void) caminho; (void) V; (void) E_total; (void) rank; (void) nprocs;
    (void) E_local; (void) rec_ini; (void) rec_fim;
    fprintf(stderr, "ler_grafo_binario_mpiio: compilado SEM -DHAVE_MPI (use mpicc).\n");
    return NULL;
}

#endif /* HAVE_MPI */

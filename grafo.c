/* DSU (Union-Find) + leitores de I/O (texto, binario, streaming e MPI). */
#define _FILE_OFFSET_BITS 64
#define _POSIX_C_SOURCE 200112L

#include "grafo.h"
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <inttypes.h>
#include <sys/types.h>

/* ===================== DSU ===================== */

int dsu_init(DSU *dsu, uint32_t n)
{
    if (dsu == NULL || n == 0)
        return 0;

    dsu->pai  = (uint32_t *) malloc((size_t) n * sizeof(uint32_t));
    dsu->rank = (uint32_t *) calloc((size_t) n, sizeof(uint32_t));
    /* [SOMA PARCIAL] calloc -> todo componente comeca com peso 0.0 (singleton). */
    dsu->peso = (double *)   calloc((size_t) n, sizeof(double));
    if (dsu->pai == NULL || dsu->rank == NULL || dsu->peso == NULL) {
        free(dsu->pai);
        free(dsu->rank);
        free(dsu->peso);
        dsu->pai = NULL;
        dsu->rank = NULL;
        dsu->peso = NULL;
        dsu->n = 0;
        return 0;
    }

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
    free(dsu->peso);            /* [SOMA PARCIAL] libera o vetor de pesos */
    dsu->pai = NULL;
    dsu->rank = NULL;
    dsu->peso = NULL;           /* [SOMA PARCIAL] */
    dsu->n = 0;
}

uint32_t dsu_find(DSU *dsu, uint32_t x)
{
    /* Iterativo (V grande -> sem recursao). 1a passada acha a raiz. */
    uint32_t raiz = x;
    while (dsu->pai[raiz] != raiz)
        raiz = dsu->pai[raiz];

    /* 2a passada: compressao de caminho (aponta tudo direto p/ a raiz). */
    while (dsu->pai[x] != raiz) {
        uint32_t prox = dsu->pai[x];
        dsu->pai[x] = raiz;
        x = prox;
    }
    return raiz;
}

int dsu_union(DSU *dsu, uint32_t a, uint32_t b, double peso_aresta)
{
    uint32_t ra = dsu_find(dsu, a);
    uint32_t rb = dsu_find(dsu, b);

    if (ra == rb)
        return 0;

    /* [SOMA PARCIAL] soma parcial do componente resultante: os dois lados mais a
     * aresta que os une. Calculada antes de decidir a raiz porque independe de
     * quem fica como raiz; depois e gravada na raiz vencedora. */
    double soma_comp = dsu->peso[ra] + dsu->peso[rb] + peso_aresta;

    /* Union by rank: pendura a arvore mais baixa sob a mais alta. */
    if (dsu->rank[ra] < dsu->rank[rb]) {
        dsu->pai[ra] = rb;
        dsu->peso[rb] = soma_comp;   /* [SOMA PARCIAL] rb virou a raiz */
    } else if (dsu->rank[ra] > dsu->rank[rb]) {
        dsu->pai[rb] = ra;
        dsu->peso[ra] = soma_comp;   /* [SOMA PARCIAL] ra virou a raiz */
    } else {
        dsu->pai[rb] = ra;
        dsu->rank[ra]++;
        dsu->peso[ra] = soma_comp;   /* [SOMA PARCIAL] ra virou a raiz */
    }
    return 1;
}

/* [SOMA PARCIAL] Peso acumulado do componente de 'x' (le sempre pela raiz). */
double dsu_peso_componente(DSU *dsu, uint32_t x)
{
    return dsu->peso[dsu_find(dsu, x)];
}

uint32_t dsu_num_componentes(DSU *dsu)
{
    uint32_t comp = 0;
    for (uint32_t i = 0; i < dsu->n; i++)
        if (dsu_find(dsu, i) == i)
            comp++;
    return comp;
}

/* ===================== I/O sequencial ===================== */

/* Pula linhas so com espacos; usada apenas para o cabecalho. */
static char *proxima_linha_nao_vazia(FILE *f, char *buf, int n)
{
    while (fgets(buf, n, f) != NULL) {
        for (const char *p = buf; *p; p++)
            if (!isspace((unsigned char) *p))
                return buf;
    }
    return NULL;
}

Aresta *ler_grafo_texto(const char *caminho, uint32_t *V, uint64_t *E)
{
    FILE *f = fopen(caminho, "r");
    if (f == NULL) {
        fprintf(stderr, "ler_grafo_texto: nao foi possivel abrir '%s'.\n", caminho);
        return NULL;
    }

    char linha[256];

    /* Cabecalho: V e E (confia no E do cabecalho p/ dimensionar o vetor). */
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

    Aresta *arestas = (Aresta *) malloc((size_t) E_cab * sizeof(Aresta));
    if (arestas == NULL) {
        fprintf(stderr, "ler_grafo_texto: falha ao alocar %" PRIu64 " arestas.\n", E_cab);
        fclose(f);
        return NULL;
    }

    /* Le "u v peso"; ignora linhas vazias/mal-formadas ou com indice fora de 0..V-1. */
    uint64_t lidas = 0;
    uint64_t ignoradas = 0;
    while (lidas < E_cab && fgets(linha, sizeof linha, f) != NULL) {
        uint32_t u, v;
        double peso;
        if (sscanf(linha, "%" SCNu32 " %" SCNu32 " %lf", &u, &v, &peso) != 3) {
            ignoradas++;
            continue;
        }
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

    *E = lidas;
    return arestas;
}

Aresta *ler_grafo_binario(const char *caminho, uint32_t *V, uint64_t *E)
{
    FILE *f = fopen(caminho, "rb");
    if (f == NULL) {
        fprintf(stderr, "ler_grafo_binario: nao foi possivel abrir '%s'.\n", caminho);
        return NULL;
    }

    /* E = tamanho/16 (fseeko/ftello p/ suportar > 2 GB). */
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

    /* fread em blocos (um fread > 2 GB pode falhar em algumas libc). */
    const uint64_t BLOCO = 16ULL * 1024 * 1024;
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

    /* V = maior indice de vertice + 1 (nao ha cabecalho). */
    uint32_t vmax = 0;
    for (uint64_t i = 0; i < n; i++) {
        if (arestas[i].u > vmax) vmax = arestas[i].u;
        if (arestas[i].v > vmax) vmax = arestas[i].v;
    }

    *V = (n > 0) ? (vmax + 1) : 0;
    *E = n;
    return arestas;
}

/* ===================== I/O binario em streaming ===================== */

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
    return (uint64_t) fread(buf, sizeof(Aresta), (size_t) cap, f);
}

int grafo_bin_info(const char *caminho, uint32_t *V, uint64_t *E)
{
    FILE *f = grafo_bin_abrir(caminho);
    if (f == NULL)
        return 0;

    const uint64_t CAP = 1024ULL * 1024;
    Aresta *buf = (Aresta *) malloc((size_t) CAP * sizeof(Aresta));
    if (buf == NULL) {
        fprintf(stderr, "grafo_bin_info: falha ao alocar buffer.\n");
        fclose(f);
        return 0;
    }

    /* Varre tudo em blocos so contando E e rastreando o maior indice. */
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
/* ===== Implementacao MPI (so com mpicc -DHAVE_MPI). ===== */
#include <mpi.h>

/* Leitor TEXTO via MPI-IO. O trabalhoso aqui e o alinhamento de bordas:
 * cada rank recebe uma fatia de bytes que pode cortar linhas ao meio. */
Aresta *ler_grafo_mpiio(const char *caminho, uint32_t *V, uint64_t *E_total,
                        int rank, int nprocs, uint64_t *E_local,
                        uint64_t *byte_ini, uint64_t *byte_fim)
{
    /* Folga lida alem da borda p/ completar a linha que cruza o limite. */
    const uint64_t MAX_LINHA = 128;

    /* (1) so o rank 0 le V, E e o offset da 1a aresta; difunde a todos. */
    uint32_t Vh = 0;
    uint64_t Eh = 0;
    uint64_t data_start = 0;

    if (rank == 0) {
        FILE *f = fopen(caminho, "r");
        if (f != NULL) {
            char linha[256];
            if (proxima_linha_nao_vazia(f, linha, sizeof linha))
                sscanf(linha, "%" SCNu32, &Vh);
            if (proxima_linha_nao_vazia(f, linha, sizeof linha))
                sscanf(linha, "%" SCNu64, &Eh);
            /* apos ler as 2 linhas, ftell = inicio da 1a aresta */
            long pos = ftell(f);
            data_start = (pos >= 0) ? (uint64_t) pos : 0;
            fclose(f);
        } else {
            fprintf(stderr, "ler_grafo_mpiio: rank 0 nao abriu '%s'.\n", caminho);
        }
    }

    MPI_Bcast(&Vh,         1, MPI_UINT32_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Eh,         1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(&data_start, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    if (Vh == 0)
        return NULL;

    *V = Vh;
    *E_total = Eh;

    /* (2) abertura coletiva e tamanho do arquivo. */
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

    /* (3) fronteiras nominais [bi, bi1) da regiao de arestas deste rank. */
    uint64_t L   = (fsize > data_start) ? (fsize - data_start) : 0;
    uint64_t bi  = data_start + (L * (uint64_t) rank)       / (uint64_t) nprocs;
    uint64_t bi1 = data_start + (L * (uint64_t) (rank + 1)) / (uint64_t) nprocs;

    /* Le 1 byte antes (prefixo) p/ saber se bi cai no inicio de uma linha,
     * e MAX_LINHA depois p/ fechar a ultima linha que cruza bi1. */
    int tem_prefixo = (rank != 0 && bi > data_start);
    uint64_t read_lo = tem_prefixo ? (bi - 1) : bi;
    uint64_t read_hi = bi1 + MAX_LINHA;
    if (read_hi > fsize) read_hi = fsize;
    uint64_t nbytes = (read_hi > read_lo) ? (read_hi - read_lo) : 0;

    *byte_ini = read_lo;
    *byte_fim = read_hi;

    /* nbytes cabe em int (count de MPI-IO); ok nesta escala (dataset ~70 MB). */
    char *buf = (char *) malloc((size_t) nbytes + 1);
    if (buf == NULL) {
        MPI_File_close(&fh);
        return NULL;
    }

    /* (4) leitura coletiva da fatia de bytes (offset/count diferem por rank). */
    MPI_Status st;
    MPI_File_read_at_all(fh, (MPI_Offset) read_lo, buf, (int) nbytes,
                         MPI_CHAR, &st);
    buf[nbytes] = '\0';
    MPI_File_close(&fh);

    /* (5) acha o inicio da 1a linha que pertence a este rank. Se o prefixo nao
     * for '\n', a linha parcial e do rank anterior -> pula ate o proximo '\n'. */
    uint64_t start_idx;
    if (rank == 0) {
        start_idx = 0;
    } else if (tem_prefixo) {
        if (buf[0] == '\n') {
            start_idx = 1;
        } else {
            uint64_t k = 1;
            while (k < nbytes && buf[k] != '\n') k++;
            start_idx = (k < nbytes) ? (k + 1) : nbytes;
        }
    } else {
        start_idx = 0;
    }

    uint64_t cap = (nbytes / 16) + 16;
    Aresta *loc = (Aresta *) malloc((size_t) cap * sizeof(Aresta));
    if (loc == NULL) {
        free(buf);
        *E_local = 0;
        return NULL;
    }
    uint64_t n = 0;
    uint64_t ignoradas = 0;

    /* Percorre so as linhas cujo inicio cai em [bi, bi1) (cada linha = 1 rank). */
    uint64_t idx = start_idx;
    while (idx < nbytes) {
        uint64_t abs_inicio = read_lo + idx;
        if (abs_inicio >= bi1)
            break;

        uint64_t j = idx;
        while (j < nbytes && buf[j] != '\n') j++;
        if (j < nbytes) buf[j] = '\0';

        /* NOTA: peso lido como uint32 aqui (difere do double do texto sequencial). */
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
            ignoradas++;
        }

        idx = j + 1;
    }

    free(buf);

    fprintf(stderr, "[rank %d] arestas locais=%" PRIu64
            " ignoradas=%" PRIu64 " bytes=[%" PRIu64 ",%" PRIu64 ")\n",
            rank, n, ignoradas, read_lo, read_hi);

    *E_local = n;
    return loc;
}

/* Leitor BINARIO via MPI-IO COLETIVO. Simples: registro de 16 bytes fixos ->
 * particao aritmetica exata, cada rank le [i*E/np, (i+1)*E/np). */
Aresta *ler_grafo_binario_mpiio(const char *caminho, uint32_t *V, uint64_t *E_total,
                                int rank, int nprocs, uint64_t *E_local,
                                uint64_t *rec_ini, uint64_t *rec_fim)
{
    const uint64_t REC = sizeof(Aresta);

    /* Tipo de 16 bytes: mantem 'count' em REGISTROS (cabe em int mesmo > 2 GB). */
    MPI_Datatype TIPO_REC;
    MPI_Type_contiguous((int) REC, MPI_BYTE, &TIPO_REC);
    MPI_Type_commit(&TIPO_REC);

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

    uint64_t r_ini = (E * (uint64_t) rank)       / (uint64_t) nprocs;
    uint64_t r_fim = (E * (uint64_t) (rank + 1)) / (uint64_t) nprocs;
    uint64_t nrec  = (r_fim > r_ini) ? (r_fim - r_ini) : 0;
    *rec_ini = r_ini;
    *rec_fim = r_fim;

    /* Buffer local = vetor de Aresta (mesmo layout do disco: le direto). */
    Aresta *loc = (Aresta *) malloc((size_t) nrec * sizeof(Aresta));
    if (loc == NULL && nrec > 0) {
        fprintf(stderr, "[rank %d] ler_grafo_binario_mpiio: malloc de %" PRIu64
                " arestas falhou.\n", rank, nrec);
        MPI_File_close(&fh);
        MPI_Type_free(&TIPO_REC);
        *E_local = 0;
        return NULL;
    }

    MPI_Status stt;
    MPI_File_read_at_all(fh, (MPI_Offset) (r_ini * REC), loc, (int) nrec,
                         TIPO_REC, &stt);
    MPI_File_close(&fh);
    MPI_Type_free(&TIPO_REC);

    /* V global = maior indice + 1 (reduz o maximo local de cada rank). */
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

/* Leitor BINARIO SEM disco compartilhado (padrao). Mesma particao do coletivo,
 * mas so o rank 0 abre o arquivo: le em blocos e ENVIA a fatia de cada rank. */
Aresta *ler_grafo_binario_dist(const char *caminho, uint32_t *V, uint64_t *E_total,
                               int rank, int nprocs, uint64_t *E_local,
                               uint64_t *rec_ini, uint64_t *rec_fim)
{
    const uint64_t REC   = sizeof(Aresta);
    const uint64_t CHUNK = 4ull * 1024 * 1024;
    const int      TAG   = 77;

    MPI_Datatype TIPO_REC;
    MPI_Type_contiguous((int) REC, MPI_BYTE, &TIPO_REC);
    MPI_Type_commit(&TIPO_REC);

    /* (1) so o rank 0 descobre E pelo tamanho; difunde. */
    uint64_t E = 0;
    FILE *f = NULL;
    if (rank == 0) {
        f = fopen(caminho, "rb");
        if (f == NULL) {
            fprintf(stderr, "ler_grafo_binario_dist: rank 0 nao abriu '%s'.\n", caminho);
        } else if (fseeko(f, 0, SEEK_END) == 0) {
            off_t sz = ftello(f);
            if (sz > 0) {
                if ((uint64_t) sz % REC != 0)
                    fprintf(stderr, "ler_grafo_binario_dist: AVISO — tamanho %lld nao e"
                            " multiplo de %" PRIu64 " (resto ignorado).\n",
                            (long long) sz, REC);
                E = (uint64_t) sz / REC;
            }
        }
    }
    MPI_Bcast(&E, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    if (E == 0) {
        if (f != NULL) fclose(f);
        MPI_Type_free(&TIPO_REC);
        return NULL;
    }
    *E_total = E;

    /* (2) particao identica a do coletivo (preserva o resultado). */
    uint64_t r_ini = (E * (uint64_t) rank)       / (uint64_t) nprocs;
    uint64_t r_fim = (E * (uint64_t) (rank + 1)) / (uint64_t) nprocs;
    uint64_t nrec  = (r_fim > r_ini) ? (r_fim - r_ini) : 0;
    *rec_ini = r_ini;
    *rec_fim = r_fim;

    Aresta *loc = (Aresta *) malloc((size_t) (nrec ? nrec : 1) * sizeof(Aresta));
    if (loc == NULL) {
        fprintf(stderr, "[rank %d] ler_grafo_binario_dist: malloc de %" PRIu64
                " arestas falhou.\n", rank, nrec);
        if (f != NULL) fclose(f);
        MPI_Type_free(&TIPO_REC);
        *E_local = 0;
        return NULL;
    }

    /* (3) transporte: rank 0 le a propria fatia e envia as demais em blocos. */
    if (rank == 0) {
        if (nrec > 0) {
            fseeko(f, (off_t) (r_ini * REC), SEEK_SET);
            if (fread(loc, REC, (size_t) nrec, f) != nrec)
                fprintf(stderr, "[rank 0] ler_grafo_binario_dist: leitura curta da propria fatia.\n");
        }
        Aresta *tmp = (Aresta *) malloc((size_t) CHUNK * sizeof(Aresta));
        if (tmp == NULL) {
            fprintf(stderr, "[rank 0] ler_grafo_binario_dist: malloc do buffer de envio falhou.\n");
            free(loc); fclose(f); MPI_Type_free(&TIPO_REC); *E_local = 0;
            return NULL;
        }
        for (int p = 1; p < nprocs; p++) {
            uint64_t pi  = (E * (uint64_t) p)       / (uint64_t) nprocs;
            uint64_t pf  = (E * (uint64_t) (p + 1)) / (uint64_t) nprocs;
            uint64_t rem = (pf > pi) ? (pf - pi) : 0;
            fseeko(f, (off_t) (pi * REC), SEEK_SET);
            while (rem > 0) {
                uint64_t blk = (rem < CHUNK) ? rem : CHUNK;
                if (fread(tmp, REC, (size_t) blk, f) != blk)
                    fprintf(stderr, "[rank 0] ler_grafo_binario_dist: leitura curta (rank %d).\n", p);
                MPI_Send(tmp, (int) blk, TIPO_REC, p, TAG, MPI_COMM_WORLD);
                rem -= blk;
            }
        }
        free(tmp);
        fclose(f);
    } else {
        /* recebe a propria fatia nos mesmos blocos em que o rank 0 enviou */
        uint64_t rem = nrec, off = 0;
        while (rem > 0) {
            uint64_t blk = (rem < CHUNK) ? rem : CHUNK;
            MPI_Status stt;
            MPI_Recv(loc + off, (int) blk, TIPO_REC, 0, TAG, MPI_COMM_WORLD, &stt);
            off += blk;
            rem -= blk;
        }
    }

    /* (4) V global = maior indice + 1 (identico ao coletivo). */
    uint32_t vmax_local = 0;
    for (uint64_t i = 0; i < nrec; i++) {
        if (loc[i].u > vmax_local) vmax_local = loc[i].u;
        if (loc[i].v > vmax_local) vmax_local = loc[i].v;
    }
    uint32_t vmax_global = 0;
    MPI_Allreduce(&vmax_local, &vmax_global, 1, MPI_UINT32_T, MPI_MAX, MPI_COMM_WORLD);
    *V = (E > 0) ? (vmax_global + 1) : 0;

    fprintf(stderr, "[rank %d] (dist) arestas locais=%" PRIu64
            " registros=[%" PRIu64 ",%" PRIu64 ") V=%" PRIu32 "\n",
            rank, nrec, r_ini, r_fim, *V);

    *E_local = nrec;
    MPI_Type_free(&TIPO_REC);
    return loc;
}

#else
/* Sem -DHAVE_MPI (alvos gcc: teste/sequencial): stubs que retornam NULL. */

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

Aresta *ler_grafo_binario_dist(const char *caminho, uint32_t *V, uint64_t *E_total,
                               int rank, int nprocs, uint64_t *E_local,
                               uint64_t *rec_ini, uint64_t *rec_fim)
{
    (void) caminho; (void) V; (void) E_total; (void) rank; (void) nprocs;
    (void) E_local; (void) rec_ini; (void) rec_fim;
    fprintf(stderr, "ler_grafo_binario_dist: compilado SEM -DHAVE_MPI (use mpicc).\n");
    return NULL;
}

#endif

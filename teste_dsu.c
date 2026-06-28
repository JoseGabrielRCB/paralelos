/*
 * teste_dsu.c — Validação do DSU em um grafo de 5 nós (base 0).
 *
 * Verifica:
 *   1. find/union básicos.
 *   2. fusão redundante retornando 0.
 *   3. contagem de componentes ao longo das fusões.
 *
 * Grafo de teste (5 vértices: 0,1,2,3,4). Aplicaremos uniões e checaremos
 * raízes e número de componentes a cada passo.
 */
#include "grafo.h"
#include <stdio.h>

/* Contadores globais simples de verificação. */
static int g_passou = 0;
static int g_falhou = 0;

/* Macro de asserção: imprime resultado e contabiliza. */
#define CHECK(cond, msg)                                            \
    do {                                                            \
        if (cond) { g_passou++; printf("  [OK]   %s\n", msg); }     \
        else      { g_falhou++; printf("  [FALHA] %s\n", msg); }    \
    } while (0)

int main(void)
{
    DSU dsu;
    const uint32_t N = 5; /* nós 0..4 */

    printf("=== Teste do DSU (grafo de %u nos, base 0) ===\n", N);

    /* Inicialização: 5 conjuntos singleton => 5 componentes. */
    CHECK(dsu_init(&dsu, N) == 1, "dsu_init retorna sucesso");
    CHECK(dsu_num_componentes(&dsu) == 5, "5 componentes apos init");

    /* Cada vértice é raiz de si mesmo logo após init. */
    CHECK(dsu_find(&dsu, 0) == 0, "find(0) == 0 (singleton)");
    CHECK(dsu_find(&dsu, 4) == 4, "find(4) == 4 (singleton)");

    /* União 0-1: deve fundir (retorno 1). Componentes: {0,1}{2}{3}{4} => 4. */
    CHECK(dsu_union(&dsu, 0, 1) == 1, "union(0,1) funde (retorna 1)");
    CHECK(dsu_find(&dsu, 0) == dsu_find(&dsu, 1), "0 e 1 na mesma raiz");
    CHECK(dsu_num_componentes(&dsu) == 4, "4 componentes apos union(0,1)");

    /* União 2-3: funde. Componentes: {0,1}{2,3}{4} => 3. */
    CHECK(dsu_union(&dsu, 2, 3) == 1, "union(2,3) funde (retorna 1)");
    CHECK(dsu_num_componentes(&dsu) == 3, "3 componentes apos union(2,3)");

    /* União 0-2: funde {0,1} com {2,3}. Componentes: {0,1,2,3}{4} => 2. */
    CHECK(dsu_union(&dsu, 0, 2) == 1, "union(0,2) funde (retorna 1)");
    CHECK(dsu_find(&dsu, 1) == dsu_find(&dsu, 3),
          "1 e 3 na mesma raiz apos cadeia de unioes");
    CHECK(dsu_num_componentes(&dsu) == 2, "2 componentes apos union(0,2)");

    /* FUSÃO REDUNDANTE: 1 e 3 já estão juntos => retorno 0, sem mudar nada. */
    CHECK(dsu_union(&dsu, 1, 3) == 0, "union(1,3) redundante (retorna 0)");
    CHECK(dsu_num_componentes(&dsu) == 2,
          "ainda 2 componentes apos fusao redundante");

    /* Vértice 4 permanece isolado e separado do grupo {0,1,2,3}. */
    CHECK(dsu_find(&dsu, 4) != dsu_find(&dsu, 0), "4 isolado do grupo");

    /* União final 4 ao grupo: 1 único componente. */
    CHECK(dsu_union(&dsu, 4, 3) == 1, "union(4,3) funde (retorna 1)");
    CHECK(dsu_num_componentes(&dsu) == 1, "1 componente apos unir tudo");

    dsu_free(&dsu);

    /* Resumo. */
    printf("\n=== Resultado: %d OK, %d FALHA ===\n", g_passou, g_falhou);
    if (g_falhou == 0) {
        printf("Todos os testes do DSU passaram.\n");
        return 0;
    }
    printf("Ha falhas no DSU!\n");
    return 1;
}

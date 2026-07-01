/* Teste unitario do DSU: find/union, fusao redundante e contagem de componentes
 * num grafo de 5 nos (base 0). */
#include "grafo.h"
#include <stdio.h>

static int g_passou = 0;
static int g_falhou = 0;

/* Asserta 'cond', imprime resultado e contabiliza. */
#define CHECK(cond, msg)                                            \
    do {                                                            \
        if (cond) { g_passou++; printf("  [OK]   %s\n", msg); }     \
        else      { g_falhou++; printf("  [FALHA] %s\n", msg); }    \
    } while (0)

int main(void)
{
    DSU dsu;
    const uint32_t N = 5;

    printf("=== Teste do DSU (grafo de %u nos, base 0) ===\n", N);

    /* init -> 5 singletons */
    CHECK(dsu_init(&dsu, N) == 1, "dsu_init retorna sucesso");
    CHECK(dsu_num_componentes(&dsu) == 5, "5 componentes apos init");

    CHECK(dsu_find(&dsu, 0) == 0, "find(0) == 0 (singleton)");
    CHECK(dsu_find(&dsu, 4) == 4, "find(4) == 4 (singleton)");

    /* union funde e reduz o numero de componentes */
    CHECK(dsu_union(&dsu, 0, 1) == 1, "union(0,1) funde (retorna 1)");
    CHECK(dsu_find(&dsu, 0) == dsu_find(&dsu, 1), "0 e 1 na mesma raiz");
    CHECK(dsu_num_componentes(&dsu) == 4, "4 componentes apos union(0,1)");

    CHECK(dsu_union(&dsu, 2, 3) == 1, "union(2,3) funde (retorna 1)");
    CHECK(dsu_num_componentes(&dsu) == 3, "3 componentes apos union(2,3)");

    /* une dois grupos de 2 */
    CHECK(dsu_union(&dsu, 0, 2) == 1, "union(0,2) funde (retorna 1)");
    CHECK(dsu_find(&dsu, 1) == dsu_find(&dsu, 3),
          "1 e 3 na mesma raiz apos cadeia de unioes");
    CHECK(dsu_num_componentes(&dsu) == 2, "2 componentes apos union(0,2)");

    /* fusao redundante: ja juntos -> retorna 0, nada muda */
    CHECK(dsu_union(&dsu, 1, 3) == 0, "union(1,3) redundante (retorna 0)");
    CHECK(dsu_num_componentes(&dsu) == 2,
          "ainda 2 componentes apos fusao redundante");

    CHECK(dsu_find(&dsu, 4) != dsu_find(&dsu, 0), "4 isolado do grupo");

    CHECK(dsu_union(&dsu, 4, 3) == 1, "union(4,3) funde (retorna 1)");
    CHECK(dsu_num_componentes(&dsu) == 1, "1 componente apos unir tudo");

    dsu_free(&dsu);

    printf("\n=== Resultado: %d OK, %d FALHA ===\n", g_passou, g_falhou);
    if (g_falhou == 0) {
        printf("Todos os testes do DSU passaram.\n");
        return 0;
    }
    printf("Ha falhas no DSU!\n");
    return 1;
}

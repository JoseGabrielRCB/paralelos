# Árvore Geradora Mínima — Borůvka (Sequencial e Paralelo com MPI)

Implementação do algoritmo de **Árvore Geradora Mínima (AGM / MST)** de **Borůvka**,
em duas versões:

- **Sequencial** (`sequencial.c`) — transcrição literal da Parte I do enunciado.
- **Paralela** (`paralelo.c`) — transcrição literal da Parte II, no **Estilo B**
  (coletivas MPI: `MPI_Allreduce` com `MPI_Op` customizada + leitura via **MPI-IO**).

A especificação seguida fielmente está em [`ALGORTIMO_REFERENCIA.md`](ALGORTIMO_REFERENCIA.md).

## Estrutura

| Arquivo | Descrição |
|---|---|
| `grafo.h` / `grafo.c` | Estruturas base: `Aresta`, **DSU** (Union-Find com compressão de caminho e union by rank) e I/O do grafo (texto e MPI-IO). |
| `sequencial.c` | Borůvka sequencial (Parte I do enunciado). |
| `paralelo.c` | Borůvka paralelo com MPI (Parte II, Estilo B). |
| `teste_dsu.c` | Teste unitário da estrutura DSU. |
| `ALGORTIMO_REFERENCIA.md` | Pseudocódigo de referência e desvios autorizados. |
| `Makefile` | Alvos de compilação. |

## Compilação

```sh
make teste        # teste do DSU            (gcc  -O2 -Wall)
make sequencial   # programa sequencial     (gcc  -O2 -Wall)
make paralelo     # programa paralelo (MPI) (mpicc -O2 -Wall -DHAVE_MPI)
make clean        # remove executáveis, objetos, logs e saídas
```

Requisitos: `gcc` e uma implementação de MPI (ex.: Open MPI, fornecendo `mpicc`/`mpirun`).

## Execução

```sh
./teste_dsu

./sequencial <arquivo_dados> [arquivo_saida]

mpirun -np <N> ./paralelo <arquivo_dados> [arquivo_saida]
```

A versão paralela deve produzir **o mesmo peso total** que a sequencial para qualquer
número de processos `N`, graças ao desempate determinístico.

## Formato dos dados de entrada

```
Linha 1: V        (número de vértices)
Linha 2: E        (número de arestas)
Linha 3..: u v peso   (uma aresta não direcionada por linha)
```

Índices de vértice em **base 0** (`0 .. V-1`). Pesos em `uint32_t`; a soma da MST é
acumulada em `uint64_t` (desvio **D2**).

### Grafo pequeno de validação (peso total esperado = **37**)

O grafo de 9 vértices e 14 arestas do enunciado serve para validar a corretude
(ver `ALGORTIMO_REFERENCIA.md`).

## Desvios autorizados em relação ao pseudocódigo literal

- **D1 — Representação:** lista de arestas em vez da matriz de pesos `W(u,v)`, que
  seria inviável (~1,96M vértices ≈ 16 TB de RAM). A semântica de `W(u,v)` e do
  `MIN{...}` por componente é preservada.
- **D2 — Tipos:** peso individual em `uint32_t`; soma total em `uint64_t`.
- **D3 — Parada:** o grafo é desconexo, então além do critério literal
  (`raiz(u)==raiz(v)` para todos) há a guarda: se uma fase não fizer nenhuma fusão e
  ainda houver mais de um componente, encerra retornando a **floresta** geradora mínima.

## Saída

1. **Distância (peso) da árvore mínima** — soma dos pesos das arestas selecionadas.
2. O **caminho** da árvore no formato `nó x → aresta → nó y → ...`.

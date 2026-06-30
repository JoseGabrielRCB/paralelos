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

### Eficiência do sequencial em streaming (`SEQ_RAM_EDGES`)

No modo streaming (`.bin` grande), o sequencial relê o arquivo a cada fase de
Borůvka — o I/O é o gargalo. Como o DSU só funde componentes, uma aresta cujos
extremos já estão na mesma componente nunca mais é útil. Definindo
`SEQ_RAM_EDGES=<n>`, assim que as arestas que **ainda cruzam** componentes
couberem em `n` posições na RAM, elas são capturadas e as fases seguintes
**param de reler o disco** (varrem só esse subconjunto). O resultado é idêntico;
serve quando há RAM sobrando:

```sh
# captura ate 80M arestas (~1,3 GB) -> corta drasticamente as releituras
SEQ_RAM_EDGES=80000000 ./sequencial graph.bin saida.txt
```

Sem a variável (padrão), o comportamento é o de só-leitura — memória O(bloco),
seguro em VMs pequenas (ex.: WSL2).

A versão paralela deve produzir **o mesmo peso total** que a sequencial para qualquer
número de processos `N`, graças ao desempate determinístico.

### Cluster por SSH sem disco compartilhado (graph.bin numa máquina só)

Quando o `graph.bin` existe em **apenas uma** das máquinas e **não há NFS** entre
os nós (cenário típico de cluster por SSH sem `sudo`), o leitor binário **padrão**
(`ler_grafo_binario_dist`) faz o **rank 0 ler o arquivo e enviar a fatia de cada
processo pela rede** (MPI_Send). O arquivo só precisa existir na máquina onde roda
o rank 0 (a primeira do `hostfile`). A **partição é a mesma** da leitura coletiva
(rank `i` = registros `[i*E/np, (i+1)*E/np)`), então o **peso total é idêntico**.

```sh
# graph.bin existe só na 1ª máquina do hostfile; nada é copiado para as outras.
mpirun -np 8 --hostfile hosts \
  --mca btl_tcp_if_include 172.26.1.0/24 \
  ./paralelo graph.bin saida.txt
```

Se **houver** disco compartilhado (NFS) com o `graph.bin` visível em todos os nós,
exporte `GRAFO_MPIIO=1` para usar a leitura **coletiva MPI-IO** (cada rank lê sua
fatia em paralelo, mais rápida nesse caso):

```sh
GRAFO_MPIIO=1 mpirun -np 8 --hostfile hosts ./paralelo graph.bin saida.txt
```

## Formato dos dados de entrada

O leitor é escolhido **pela extensão** do arquivo passado em `argv[1]`:

### Texto (qualquer extensão != `.bin`) — usado nos testes de corretude

```
Linha 1: V        (número de vértices)
Linha 2: E        (número de arestas)
Linha 3..: u v peso   (uma aresta não direcionada por linha)
```

### Binário (`*.bin`) — arquivo de execução real (`graph.bin`, 12,8 GB)

Sequência de registros contíguos de **16 bytes**, **sem cabeçalho**:

| campo | tipo | bytes |
|---|---|---|
| `u`    | `uint32` little-endian | 4 |
| `v`    | `uint32` little-endian | 4 |
| `peso` | `double` (IEEE-754) little-endian | 8 |

Logo `E = tamanho_do_arquivo / 16` e, como não há cabeçalho, `V = (maior índice de
vértice) + 1` (no paralelo, via `MPI_Allreduce(MAX)`). O layout casa exatamente com
`struct Aresta`, então a leitura é cópia direta de bytes (sem parsing). No paralelo,
a partição é aritmética exata por registro — sem alinhamento de bordas de linha.

Índices de vértice em **base 0** (`0 .. V-1`). Pesos em `double`; a soma da MST
também é acumulada em `double` (desvio **D2**).

### Grafo pequeno de validação (peso total esperado = **37**)

O grafo de 9 vértices e 14 arestas do enunciado serve para validar a corretude
(ver `ALGORTIMO_REFERENCIA.md`).

## Desvios autorizados em relação ao pseudocódigo literal

- **D1 — Representação:** lista de arestas em vez da matriz de pesos `W(u,v)`, que
  seria inviável (~10M vértices ≈ 800 TB de RAM). A semântica de `W(u,v)` e do
  `MIN{...}` por componente é preservada.
- **D2 — Tipos:** peso individual em `double` (o `graph.bin` traz pesos em ponto
  flutuante); soma total também em `double`.
- **D3 — Parada:** o grafo é desconexo, então além do critério literal
  (`raiz(u)==raiz(v)` para todos) há a guarda: se uma fase não fizer nenhuma fusão e
  ainda houver mais de um componente, encerra retornando a **floresta** geradora mínima.

## Saída

1. **Distância (peso) da árvore mínima** — soma dos pesos das arestas selecionadas.
2. O **caminho** da árvore no formato `nó x → aresta → nó y → ...`.

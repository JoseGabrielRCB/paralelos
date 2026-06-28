# Algoritmo de referência (transcrição fiel do enunciado)

> Transcrição literal dos pseudocódigos do `arvore_geradora_minima.pdf` (onde aparecem como
> imagem). **Esta é a especificação a ser seguida.** Não alterar a lógica; desvios só com
> anotação e justificativa, conforme a política do projeto.

---

## Parte I — Algoritmo Sequencial (pseudocódigo literal)

```
algoritmo Árvore Geradora Mínima
    entrada: Um grafo não direcionado ponderado G = (V, E).
    saída: F, uma floresta estendida mínima de G.

    Inicialize uma floresta F para (V, E') onde E' = {}.

    concluído := false
    enquanto não concluído faça
        Encontre os componentes conectados de F e atribua a cada vértice seu componente
        Inicialize a aresta de menor peso para cada componente como "Nenhuma"
        para cada aresta uv em E, onde uv estão em diferentes componentes de F:
            deixe wx ser a aresta de menor peso para a componente de u
            se for-preferido(uv, wx) então
                Defina uv como a aresta de menor peso para a componente de u
            deixe yz ser a aresta de menor peso para a componente de v
            se for-preferido-sobre(uv, yz) então
                Defina uv como a aresta de menor peso para o componente de v
        se todos os componentes tiverem a aresta de menor peso definida como "Nenhum" então
            // nenhuma outra árvore pode ser mesclada -- nós terminamos
            complete := true
        else
            complete := false
            para cada componente cuja aresta de menor peso não é "Nenhuma" faça
                Adicione sua aresta de menor peso em E'

função for-preferido-sobre(aresta1, aresta2)
    return (aresta2 is "None") ou
           (peso(aresta1) < peso(aresta2)) ou
           (peso(aresta1) = peso(aresta2) e regra-desempate(aresta1, aresta2))

função regra-desempate(aresta1, aresta2)
    a regra de desempate retorna true se e somente se aresta1
    for-preferido-sobre aresta2 no caso de um empate.
```

**Notas de leitura (não alteram a lógica):** este é o algoritmo de **Borůvka**. Cada iteração
do `enquanto` é uma fase: acha-se, para cada componente, a aresta de menor peso que o liga a
outro componente (com desempate determinístico via `for-preferido-sobre`), e essas arestas são
adicionadas de uma vez a `E'`. O `regra-desempate` deve dar uma ordem total (ex.: menor peso;
em empate, menor par (u,v)) para o resultado ser reproduzível.

---

## Parte II — Algoritmo Paralelo (pseudocódigo literal)

```
ALGORITMO PARALELO
Entrada: Grafo não direcionado ponderado G = (V, E) representado por uma matriz de pesos W(u,v).
Saída: Árvore Geradora de Custo Mínimo de G.

Início
    para cada u ∈ V faça em paralelo
        raiz(i) := i
    fim_paralelo
    finalizado := Falso
    enquanto not (finalizado) faça
        para cada u ∈ V faça em paralelo
            vertice+proximo(u) := v tal que raiz(u) seja diferente de raiz(v)
                                  e W(u,v) := MIN{W(u,w) tal que w ∈ V}
        fim_paralelo
        para cada componente k da floresta Fi faça em paralelo
            escolha um vertice u tal que W(u, vertice+proximo(u)) seja o menor vertice de k
        fim_paralelo
        combine os novos vertices e crie a nova floresta Fi+1
        para cada vertice u ∈ V faça em paralelo
            procure pela raiz(u)
            se raiz(u) == raiz(v) para todo u, v ∈ V então
                finalizado := True
        fim_paralelo
    fim_enquanto
Fim
```

**Estruturas citadas no enunciado (Parte II / seção "Explicando o Algoritmo Paralelo"):**
- `raiz` — dá a raiz da árvore (componente) à qual um nó pertence (= Union-Find / DSU).
- `vertice+proximo(u)` — o vértice de OUTRA árvore que está mais próximo (aresta de menor peso) de `u`.
- A floresta `Fi` contém várias árvores; cada árvore é representada pela relação de seus pais.

**Saída exigida (duas partes):**
1. A **distância da árvore mínima** = somatória dos pesos das arestas da árvore.
2. O **caminho** da árvore, no formato: `nó x → aresta → nó y → aresta → nó z → ...`

---

## Desvios autorizados sobre a Parte II (com anotação obrigatória no código)

- **D1 (representação):** a "matriz de pesos W(u,v)" é inviável (~1,96M vértices ≈ 16 TB de RAM).
  Usar **lista de arestas**, preservando a semântica de `W(u,v)` = peso da aresta (u,v) e o
  `MIN{...}` por componente. Comentar a justificativa.
- **D3 (parada):** o critério literal é `raiz(u)==raiz(v)` para todos (1 componente). Como o
  grafo é **desconexo** (maior componente conexo ≈ 99,6% dos nós), implementar o laço literal e
  adicionar a guarda: se uma fase não fizer nenhuma fusão e houver >1 componente, encerrar
  retornando a **floresta** geradora mínima, comentando que a condição literal não foi atendida.
- **D2 (tipos):** peso em `uint32_t`; soma total em `uint64_t`.
- **Estilo MPI (Estilo B):** o "faça em paralelo" é realizado por processos MPI que particionam
  as arestas; a combinação dos mínimos por componente usa `MPI_Allreduce` + `MPI_Op` customizada
  (mesmo desempate determinístico). Leitura/distribuição via **MPI-IO**.

---

## Formato dos dados de entrada (verificado)

Arquivo `dados_entrada_sequencial.txt`:
- Linha 1: número de vértices V = **1965206**
- Linha 2: número de arestas E = **2769244** (cabeçalho do próprio arquivo; o arquivo de
  características cita 2766607 — confiar SEMPRE no cabeçalho do arquivo de dados)
- Da linha 3 em diante: `u v peso` (um por linha)

Características verificadas:
- Índices de vértice são **base 0**: variam de **0** a **1965205** (portanto DSU/arrays de
  tamanho V, índices 0..V-1).
- Pesos vão além de 4×10⁹ (ex.: 4185718553) → **não** cabem em `int` com sinal; usar `uint32_t`
  e somar em `uint64_t`.
- Grafo **não direcionado** e **desconexo**.

---

## Grafo pequeno para teste de corretude (resultado conhecido: peso total = 37)

Grafo do próprio enunciado (9 vértices, 14 arestas). Arestas `u v peso`:

```
9
14
0 1 4
0 7 8
1 2 8
1 7 11
2 3 7
2 8 2
2 5 4
3 4 9
3 5 14
4 5 10
5 6 2
6 7 1
6 8 6
7 8 7
```

MST esperada (peso total = **37**): arestas {6-7:1, 2-8:2, 5-6:2, 0-1:4, 2-5:4, 2-3:7, 0-7:8, 3-4:9}.
Use este grafo para validar o sequencial e, depois, comparar a saída do paralelo (deve dar 37
para qualquer número de processos).
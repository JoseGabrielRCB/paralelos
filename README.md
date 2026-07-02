# Árvore Geradora Mínima — Borůvka (Sequencial e Paralelo com MPI)

Algoritmo de **Árvore Geradora Mínima (AGM / MST)** de **Borůvka** em duas versões
que dão **o mesmo peso total** para qualquer número de processos:

- **`sequencial.c`** — roda numa máquina só (`gcc`).
- **`paralelo.c`** — roda em várias máquinas via **MPI** (`mpicc`/`mpirun`).

Este README é focado em **como rodar** e **como usar o Makefile** (inclusive num
cluster por SSH). No fim há um resumo de como o código funciona.

---

## Requisitos

- `gcc` (para o sequencial e o teste).
- Uma implementação de **MPI** — ex.: **Open MPI**, que fornece `mpicc` e `mpirun`
  (para o paralelo).
- Para o cluster: acesso **SSH** às máquinas e o **mesmo caminho** (`REMOTE_PATH`)
  existindo em todas elas (o home é local, não compartilhado).

---

## O que mudar no Makefile para rodar

Todas as configurações ficam no bloco **`CONFIGURACAO`** no topo do `Makefile`.
Edite estas variáveis conforme o seu ambiente:

| Variável | O que é | Troque por |
|---|---|---|
| `RGM` | Seu login nas máquinas do cluster | O seu usuário (ex.: `rgm12345`) |
| `MAE` | Nó que roda o **rank 0** (vira a 1ª linha do hostfile) | O nome da sua máquina "mãe" |
| `MAQUINAS` | Os **demais** nós (separados por espaço; **não** repita o `MAE`) | A lista das suas máquinas |
| `REMOTE_PATH` | Caminho **idêntico** em todas as máquinas onde o binário é copiado | Ex.: `/home/local/$(RGM)/paralelos` |
| `GRAFO` | Arquivo de entrada (fica na `MAE` quando se usa o leitor distribuído) | Caminho do seu grafo (`.bin` ou `.txt`) |
| `SAIDA` | Arquivo de **saída** (a lista de arestas da MST) | O nome que quiser (ex.: `saida.txt`) |
| `NP` | Nº de processos padrão do `make run` | Quantos processos quer (ou passe na linha) |
| `IFACE` | Rede/interface comum dos nós (evita `docker0`/`virbr0`) | A sua sub-rede (ex.: `172.26.1.0/24`) ou `eth0` |

> `GRAFO` e `SAIDA` também podem ser passados na hora: `make run-seq GRAFO=g37.txt SAIDA=out.txt`.

---

## Rodar no cluster — passo a passo

Rode a partir da máquina **mãe** (`MAE`), dentro de `REMOTE_PATH`.

### 1. `make keys` — libera SSH sem senha (só na primeira vez)

```sh
make keys
```

Gera uma chave SSH (se ainda não houver) e a copia (`ssh-copy-id`) para a `MAE` e
todas as `MAQUINAS`. A partir daí o `mpirun` e o `scp` conectam sem pedir senha.

### 2. `make` — compila o paralelo **e distribui** para os nós

```sh
make
```

É o alvo padrão (`all = compila + envia`): compila o `paralelo` com
`mpicc -DHAVE_MPI` e **copia o binário** para `REMOTE_PATH` em cada máquina
(o home é local, então o executável precisa existir no mesmo caminho em todas).

### 3. `make run` — executa o paralelo

```sh
make run       # usa NP (padrão 4)
make run 8     # atalho: roda com 8 processos (NP=8)
```

Gera o `hosts` (hostfile) a partir de `MAE` + `MAQUINAS` e chama:

```sh
mpirun -np <NP> --hostfile hosts \
  --mca btl_tcp_if_include <IFACE> --mca oob_tcp_if_include <IFACE> \
  ./paralelo <GRAFO> <SAIDA>
```

Fixar a `IFACE` nos dois canais (dados e controle) evita que o MPI tente usar
`docker0`/`virbr0` e trave a conexão.

### 4. `make run-seq` — confere o resultado numa máquina só

```sh
make run-seq                                  # usa GRAFO e SAIDA do Makefile
make run-seq GRAFO=/tmp/graph.bin SAIDA=seq.txt
```

Compila (se preciso) e roda `./sequencial <GRAFO> <SAIDA>`. Serve para **conferir**
que o paralelo dá o mesmo **peso total** que o sequencial.

---

## Todos os alvos do Makefile

| Comando | O que faz |
|---|---|
| `make keys` | Gera/instala a chave SSH em todas as máquinas (1ª vez). |
| `make` / `make all` | Compila o paralelo **e** distribui o binário (`compila` + `envia`). |
| `make compila` | Só compila o `paralelo` (`mpicc -DHAVE_MPI`). |
| `make envia` | Só distribui o binário `paralelo` para os nós. |
| `make run` / `make run 8` | Roda o paralelo com `NP` processos (`8` = atalho p/ `NP=8`). |
| `make seq` | Compila só o `sequencial` (`gcc`). |
| `make run-seq` | Roda o `sequencial` com `GRAFO`/`SAIDA`. |
| `make teste` | Compila o teste do DSU (`teste_dsu`). Ver nota abaixo. |
| `make clean` | Remove executáveis, hostfile, objetos, logs e saídas. |

> **Nota:** `teste_dsu.c` não faz parte do repositório (é mantido localmente), então
> `make teste` só funciona se o arquivo estiver presente na pasta.

---

## Rodar sem o Makefile (local, para testar)

Numa máquina só, sem cluster:

```sh
# compilar
gcc   -O2 -Wall -o sequencial sequencial.c grafo.c
mpicc -O2 -Wall -DHAVE_MPI -o paralelo   paralelo.c   grafo.c

# rodar
./sequencial <arquivo_dados> [arquivo_saida]
mpirun -np <N> ./paralelo <arquivo_dados> [arquivo_saida]
```

Para testar o paralelo com mais processos que núcleos na mesma máquina, use
`mpirun --oversubscribe -np <N> ...`.



## Formato dos dados de entrada

O leitor é escolhido **pela extensão** do arquivo:

- **Texto** (extensão ≠ `.bin`): linha 1 = `V`, linha 2 = `E`, demais linhas
  `u v peso` (uma aresta por linha). Índices em base 0.
- **Binário** (`*.bin`): registros contíguos de **16 bytes**, sem cabeçalho —
  `u` (`uint32`), `v` (`uint32`), `peso` (`double`), little-endian. Assim
  `E = tamanho/16` e `V = maior índice + 1`.

Arquivos `.bin` grandes: o sequencial entra em **streaming** automaticamente acima
de 2 GB (ou com `--stream`); no cluster, o paralelo usa por padrão o leitor
**distribuído** (só a `MAE` precisa ter o arquivo). Com NFS, `GRAFO_MPIIO=1` usa a
leitura coletiva MPI-IO.

---

## Como o código funciona (resumo)

- **`grafo.c`/`grafo.h`** — base compartilhada: a struct `Aresta` (16 bytes), o
  **DSU (Union-Find)** com compressão de caminho e união por rank, e os leitores de
  I/O. O DSU é **ponderado**: um vetor `peso[raiz]` acumula a soma dos pesos da MST
  dentro de cada componente (`peso[nova_raiz] = peso[A] + peso[B] + peso_da_aresta`)
  — a **soma parcial**. O peso total é a soma dos `peso[raiz]`, usada como
  verificação cruzada (linha `Soma parcial (DSU)`).
- **`sequencial.c`** — Borůvka em uma máquina. Em cada fase acha a menor aresta de
  saída de cada componente e aplica as uniões. Modo **em-RAM** (grafo inteiro em
  memória) ou **streaming** (relê o `.bin` em blocos) para arquivos gigantes.
- **`paralelo.c`** — Borůvka com MPI no **Estilo B**: o DSU é **replicado** em todos
  os ranks; cada rank varre só sua fatia de arestas, um `MPI_Allreduce` com `MPI_Op`
  customizada combina as melhores arestas, e todos aplicam as **mesmas uniões na
  mesma ordem**. Por isso o peso total independe do número de processos.

O desempate entre arestas de peso igual segue uma **ordem total** (menor peso; em
empate, menor par `(min,max)`), idêntica nas duas versões — é o que garante o
**mesmo resultado** em seq e paralelo, para qualquer `N`.

## Saída

1. No terminal: **`Peso total da MST: <soma>`** e a linha `Soma parcial (DSU)`.
2. No arquivo `SAIDA`: a MST como lista de arestas `u → (peso) → v`, uma por linha.

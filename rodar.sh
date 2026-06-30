#!/usr/bin/env bash
# rodar.sh — lança o Borůvka PARALELO no cluster usando os NOMES das máquinas
# (via hostfile) e fixando a rede, para a conexão MPI "simplesmente funcionar"
# sem ter de digitar as flags --mca toda vez.
#
# Uso:
#   ./rodar.sh <arquivo_dados> [arquivo_saida]
#
# Exemplos:
#   ./rodar.sh /tmp/graph.bin
#   NP=16 ./rodar.sh /tmp/graph.bin resultado.txt
#
# Variáveis de ambiente (todas opcionais):
#   HOSTFILE   arquivo de hosts            (padrão: ./hostfile, ao lado do script)
#   REDE       sub-rede OU interface dos nós (padrão: 172.26.1.0/24; pode ser "eth0")
#   NP         nº de processos             (padrão: soma dos slots do hostfile)
#   MPIIO=1    usa leitura coletiva MPI-IO (só se o arquivo estiver visível em todos
#              os nós, ex.: pasta compartilhada). Sem isso, usa o leitor distribuído
#              (rank 0 lê e envia as fatias pela rede; arquivo só no 1º nó do hostfile).
set -euo pipefail

DIR="$(cd "$(dirname "$0")" && pwd)"
HOSTFILE="${HOSTFILE:-$DIR/hostfile}"
REDE="${REDE:-172.26.1.0/24}"
ENTRADA="${1:?uso: $0 <arquivo_dados> [arquivo_saida]}"
SAIDA="${2:-saida.txt}"

if [ ! -f "$HOSTFILE" ]; then
    echo "ERRO: hostfile '$HOSTFILE' nao encontrado." >&2
    echo "Crie um com uma maquina por linha, ex.:" >&2
    echo "  l1m20u24 slots=4" >&2
    echo "  l1m21u24 slots=4" >&2
    exit 1
fi

if [ ! -x "$DIR/paralelo" ]; then
    echo "ERRO: binario '$DIR/paralelo' nao existe. Compile com 'make paralelo'." >&2
    exit 1
fi

# nº de processos: usa NP se definido; senão soma os slots declarados no hostfile.
if [ -n "${NP:-}" ]; then
    NPROC="$NP"
else
    NPROC="$(awk '{for(i=1;i<=NF;i++) if($i ~ /^slots=/){sub(/slots=/,"",$i); s+=$i}} END{print (s>0?s:1)}' "$HOSTFILE")"
fi

# 1º nó do hostfile = onde roda o rank 0 (e onde o arquivo precisa estar no modo distribuído).
RANK0="$(awk 'NF && $1 !~ /^#/ {print $1; exit}' "$HOSTFILE")"

echo ">> hostfile : $HOSTFILE   (rank 0 = $RANK0)"
echo ">> rede     : $REDE"
echo ">> processos: $NPROC"
echo ">> entrada  : $ENTRADA"
echo ">> saida    : $SAIDA"
[ -n "${MPIIO:-}" ] && echo ">> leitura  : MPI-IO coletiva (GRAFO_MPIIO=1)"

# GRAFO_MPIIO precisa chegar aos processos em TODOS os nós -> -x repassa a variável.
EXTRA=()
if [ -n "${MPIIO:-}" ]; then EXTRA+=(-x GRAFO_MPIIO=1); fi

exec mpirun --hostfile "$HOSTFILE" -np "$NPROC" \
    --mca oob_tcp_if_include "$REDE" \
    --mca btl_tcp_if_include "$REDE" \
    "${EXTRA[@]}" \
    "$DIR/paralelo" "$ENTRADA" "$SAIDA"

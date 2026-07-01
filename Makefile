# Makefile — Projeto AGM (Borůvka): sequencial (gcc) + paralelo (MPI) no cluster.
#
# Fluxo tipico no cluster (rode a partir do no MAE, dentro de REMOTE_PATH):
#   make keys        # 1x: gera chave SSH e copia p/ todas as maquinas (sem senha)
#   make             # compila o paralelo e DISTRIBUI o binario p/ todos os nos
#   make run         # roda o paralelo (NP processos). Ex.: make run 8  -> NP=8
#   make run-seq     # roda o sequencial (1 maquina) p/ conferir o peso
#
# Outros alvos: make compila | envia | seq | teste | clean

# ----------------------------- CONFIGURACAO --------------------------------
RGM         := rgm48936
# MAE = no que roda o rank 0 (vira a 1a linha do hostfile).
MAE         := l1m19u24
# MAQUINAS = demais nos (acrescente os que quiser; NAO repita o MAE aqui).
MAQUINAS    := l1m21u24 l1m22u24 l1m18u24 l1m14u24 l1m28u24 l1m30u24 l1m31u24
# REMOTE_PATH = MESMO caminho em todas as maquinas (home e' local, nao compartilhado).
REMOTE_PATH := /home/local/$(RGM)/paralelos

parBIN      := paralelo
parFONTE    := paralelo.c
seqBIN      := sequencial
seqFONTE    := sequencial.c
testeBIN    := teste_dsu
testeFONTE  := teste_dsu.c
# COMUM = fonte compartilhado (DSU + I/O).
COMUM       := grafo.c
HDR         := grafo.h

HOSTFILE    := hosts
NP          := 4
# GRAFO = arquivo de entrada (no MAE, quando se usa o leitor distribuido).
GRAFO       := /tmp/graph.bin
# SAIDA = 2o argumento do programa (arquivo de SAIDA, NAO um numero).
SAIDA       := saida.txt
# IFACE = rede comum dos nos (evita docker0/virbr0). Pode ser "eth0" ou o CIDR.
IFACE       := 172.26.1.0/24

CC          := gcc
MPICC       := mpicc
CFLAGS      := -O2 -Wall -D_FILE_OFFSET_BITS=64

CHAVE       := $(HOME)/.ssh/id_rsa
SSH_OPTS    := -o StrictHostKeyChecking=no
SCP_OPTS    := $(SSH_OPTS) -o BatchMode=yes

# Permite "make run 8" -> NP=8 (o "8" vira um alvo no-op).
ifeq (run,$(firstword $(MAKECMDGOALS)))
RUN_NP := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
ifneq ($(RUN_NP),)
    NP := $(RUN_NP)
    $(eval $(RUN_NP):;@true)
endif
endif

SHELL := /bin/bash
.ONESHELL:
.DEFAULT_GOAL := all
.PHONY: all keys compila envia run seq run-seq teste clean

all: compila envia

# ------------------------------- SSH KEYS ----------------------------------
keys:
	@set -e
	test -f $(CHAVE) || ssh-keygen -t rsa -b 4096 -N "" -f $(CHAVE)
	for m in $(MAE) $(MAQUINAS); do
	  echo "--- $$m ---"
	  ssh-copy-id $(SSH_OPTS) -i $(CHAVE).pub "$(RGM)@$$m" || echo "  falhou: $$m"
	done

# ------------------------------- COMPILAR ----------------------------------
# IMPORTANTE: o paralelo precisa de grafo.c E de -DHAVE_MPI (ativa a implementacao
# MPI-IO / distribuida em grafo.c; sem a macro fica o stub que retorna NULL).
compila: $(parBIN)

$(parBIN): $(parFONTE) $(COMUM) $(HDR)
	$(MPICC) $(CFLAGS) -DHAVE_MPI -o $(parBIN) $(parFONTE) $(COMUM)
	@echo "OK: make run   (ou: mpirun -np N --hostfile $(HOSTFILE) ./$(parBIN) <grafo> <saida>)"

# ------------------------------- DISTRIBUIR --------------------------------
# Home local (nao compartilhado) => o binario precisa existir no MESMO caminho
# em cada no. Copia p/ todos os MAQUINAS (no MAE e' apenas um cp local).
envia: $(parBIN)
	@set -e
	mkdir -p $(REMOTE_PATH)
	[ "$$(readlink -f $(parBIN))" = "$$(readlink -f $(REMOTE_PATH)/$(parBIN))" ] || cp -f $(parBIN) $(REMOTE_PATH)/
	for m in $(MAQUINAS); do
	  echo "--- $$m ---"
	  ssh $(SCP_OPTS) "$(RGM)@$$m" "mkdir -p $(REMOTE_PATH)"
	  scp $(SCP_OPTS) $(parBIN) "$(RGM)@$$m:$(REMOTE_PATH)/"
	done

# ------------------------------- HOSTFILE ----------------------------------
$(HOSTFILE): Makefile
	@set -e
	: > $(HOSTFILE)
	echo "$(MAE) slots=1" >> $(HOSTFILE)
	for m in $(MAQUINAS); do echo "$$m slots=1" >> $(HOSTFILE); done

# --------------------------------- RUN -------------------------------------
# Fixa a interface/rede nos dois canais (oob = controle, btl = dados) para a
# conexao MPI nao cair em docker0/virbr0. (Open MPI: nada de lamboot.)
run: $(HOSTFILE)
	@set -e
	mpirun -np $(NP) --hostfile $(HOSTFILE) \
	  --mca btl_tcp_if_include $(IFACE) --mca oob_tcp_if_include $(IFACE) \
	  ./$(parBIN) $(GRAFO) $(SAIDA)

# ------------------------------ SEQUENCIAL ---------------------------------
seq: $(seqBIN)

$(seqBIN): $(seqFONTE) $(COMUM) $(HDR)
	$(CC) $(CFLAGS) -o $(seqBIN) $(seqFONTE) $(COMUM)

run-seq: $(seqBIN)
	./$(seqBIN) $(GRAFO) $(SAIDA)

# --------------------------------- TESTE -----------------------------------
teste: $(testeBIN)

$(testeBIN): $(testeFONTE) $(COMUM) $(HDR)
	$(CC) $(CFLAGS) -o $(testeBIN) $(testeFONTE) $(COMUM)

# --------------------------------- CLEAN -----------------------------------
clean:
	rm -f $(parBIN) $(seqBIN) $(testeBIN) $(HOSTFILE)

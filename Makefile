# Makefile — Projeto Árvore Geradora Mínima (Borůvka)
#
# Alvos:
#   make teste       -> compila e prepara o teste do DSU      (gcc  -O2 -Wall)
#   make sequencial  -> compila o programa sequencial         (gcc  -O2 -Wall)
#   make paralelo    -> compila o programa paralelo (MPI)     (mpicc -O2 -Wall)
#   make clean       -> remove executaveis e objetos
#
# Fase 0: apenas a BASE (DSU + teste). Os programas `sequencial` e `paralelo`
# serao adicionados em fases posteriores; seus alvos ja estao prontos com os
# compiladores corretos e avisam caso o fonte ainda nao exista.

CC      := gcc
MPICC   := mpicc
CFLAGS  := -O2 -Wall
COMMON  := grafo.c

.PHONY: all teste sequencial paralelo clean

# Por padrao, constroi o que existe na Fase 0: o teste do DSU.
all: teste

# --- Teste do DSU (gcc) ---
teste: teste_dsu.c $(COMMON) grafo.h
	$(CC) $(CFLAGS) -o teste_dsu teste_dsu.c $(COMMON)
	@echo "OK: execute ./teste_dsu"

# --- Programa sequencial (gcc) ---
sequencial: sequencial.c $(COMMON) grafo.h
	$(CC) $(CFLAGS) -o sequencial sequencial.c $(COMMON)
	@echo "OK: execute ./sequencial <arquivo_dados> [arquivo_saida]"

# --- Programa paralelo (mpicc / Open MPI) ---
# -DHAVE_MPI ativa a implementacao MPI-IO de ler_grafo_mpiio em grafo.c
# (nos alvos com gcc, sem essa macro, fica o stub).
paralelo: paralelo.c $(COMMON) grafo.h
	$(MPICC) $(CFLAGS) -DHAVE_MPI -o paralelo paralelo.c $(COMMON)
	@echo "OK: mpirun -np <N> ./paralelo <arquivo_dados> [arquivo_saida]"

# --- Limpeza ---
# Remove executaveis, objetos, logs por rank e saidas geradas.
clean:
	rm -f teste_dsu sequencial paralelo *.o log_rank*.txt

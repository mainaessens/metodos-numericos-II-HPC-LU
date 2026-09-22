# ============================================================================
#  Makefile — Factorización LU Paralela
#  Proyecto 6 · Métodos Numéricos II y Computación Científica · FACET-UNT 2026
#
#  Uso:
#      make            compila los dos binarios
#      make test       corrida rápida de validación serial vs MPI
#      make clean      borra binarios y objetos
#
#  En el cluster, si -march=native da problemas (nodo de compilación distinto
#  al de ejecución), compilar con:  make ARCH=
# ============================================================================

CC      ?= gcc
MPICC   ?= mpicc

# Banderas de optimización. -O3 habilita vectorización automática del lazo
# interno de la eliminación, que es donde está todo el tiempo.
# IMPORTANTE: el serial y el MPI se compilan con LAS MISMAS banderas.
# Si no, el speedup que midan no significa nada.
ARCH    ?= -march=native
OPT     ?= -O3 $(ARCH)
WARN    ?= -Wall -Wextra -Wno-unused-parameter
STD     ?= -std=c11

CFLAGS  ?= $(OPT) $(WARN) $(STD) -Iinclude
LDLIBS  ?= -lm

SRCDIR   = src
BINDIR   = bin
OBJDIR   = build

COMMON_SRC = $(SRCDIR)/matgen.c $(SRCDIR)/verify.c
COMMON_OBJ = $(OBJDIR)/matgen.o $(OBJDIR)/verify.o

.PHONY: all clean test dirs

all: dirs $(BINDIR)/lu_serial $(BINDIR)/lu_mpi

dirs:
	@mkdir -p $(BINDIR) $(OBJDIR) results

$(OBJDIR)/%.o: $(SRCDIR)/%.c include/lu.h | dirs
	$(CC) $(CFLAGS) -c $< -o $@

$(BINDIR)/lu_serial: $(SRCDIR)/lu_serial.c $(COMMON_OBJ) include/lu.h | dirs
	$(CC) $(CFLAGS) $(SRCDIR)/lu_serial.c $(COMMON_OBJ) -o $@ $(LDLIBS)

# El MPI se compila con mpicc, que es un wrapper de gcc que agrega los
# includes y las librerías de MPI. Los .o comunes se pueden reutilizar
# porque no dependen de MPI.
$(BINDIR)/lu_mpi: $(SRCDIR)/lu_mpi.c $(COMMON_OBJ) include/lu.h | dirs
	$(MPICC) $(CFLAGS) $(SRCDIR)/lu_mpi.c $(COMMON_OBJ) -o $@ $(LDLIBS)

# ---------------------------------------------------------------------------
#  Validación rápida: el resultado paralelo tiene que coincidir con el serial
#  para cualquier número de procesos, incluidos los que NO dividen a n.
# ---------------------------------------------------------------------------
test: all
	@echo ""
	@echo "=== SERIAL (referencia) ==============================="
	@$(BINDIR)/lu_serial -n 300 -k diagdom -s 7
	@echo "=== MPI con 1,2,3,4,5,7 procesos ======================"
	@for np in 1 2 3 4 5 7; do \
		echo "--- p = $$np ---"; \
		mpirun --oversubscribe -n $$np $(BINDIR)/lu_mpi -n 300 -k diagdom -s 7 2>/dev/null \
		  | grep -E "residuo|T total"; \
	done
	@echo "=== reparto en bloques vs ciclico (p=4) ==============="
	@mpirun --oversubscribe -n 4 $(BINDIR)/lu_mpi -n 600 -d cyclic --csv
	@mpirun --oversubscribe -n 4 $(BINDIR)/lu_mpi -n 600 -d block  --csv
	@echo "=== pivoteo: la matriz con diagonal nula =============="
	@echo "con pivoteo (debe andar):"
	@$(BINDIR)/lu_serial -n 200 -k zerodiag | grep -E "residuo"
	@echo "sin pivoteo (debe fallar):"
	@-$(BINDIR)/lu_serial -n 200 -k zerodiag --nopivot
	@echo ""

clean:
	rm -rf $(BINDIR) $(OBJDIR)

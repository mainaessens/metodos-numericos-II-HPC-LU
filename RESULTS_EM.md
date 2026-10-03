# Resultados — corridas de Ezequiel

Corridas de la implementación en dos sistemas: una PC personal y el cluster
de la cátedra. Los datos crudos, los resúmenes y las figuras están en
`results/<sistema>/` y `figuras/<sistema>/`.

---

## Sistemas

| | PC | Cluster |
|---|---|---|
| CPU | AMD Ryzen 7 7735HS | AMD EPYC 7B12 (2 sockets × 64 núcleos) |
| Núcleos usados | 8 físicos / 16 hilos | 48 asignados (sin hyperthreading) |
| Cache L3 | 16 MiB, una sola, compartida por todos los núcleos | 16 MiB por grupo de núcleos (512 MiB en total) |
| RAM | 6.6 GiB visibles en WSL2 | 192 GB asignados |
| Entorno | Ubuntu 24.04 en WSL2 | Ubuntu 24.04, VS Code en el navegador |
| Compilador | gcc 13.3.0 (`-O3 -march=native`) | gcc 13.3.0 (`-O3 -march=native`) |
| MPI | Open MPI 4.1.6 | Open MPI 4.1.6 |

Detalle completo de cada máquina en `results/pc/info_pc.txt` y
`results/cluster/info_cluster.txt`. El `lscpu` del cluster muestra la máquina
física entera (128 núcleos, 2 nodos NUMA), no solo los 48 núcleos asignados.

---

## Qué se corrió

En todos los casos: matriz diagonal dominante (`-k diagdom`), reparto cíclico,
semilla 42 y 5 repeticiones por punto. Se reporta la **mediana** de las
repeticiones. T₁ es la versión serial pura, no MPI con p = 1.

| Sistema | Experimento | N | p |
|---|---|---|---|
| PC | Escalabilidad fuerte | 1000, 2000, 4000 | serial, 2, 4, 8 |
| PC | Precisión | 100, 200, 500, 1000 · 4 tipos de matriz | serial y 4 |
| Cluster | Escalabilidad fuerte | 1000, 2000, 4000, 8000 | serial, 2, 4, 8, 16, 32 |

La precisión se corrió solo en la PC: el residuo depende del algoritmo y de la
aritmética IEEE 754, no del hardware.

Comandos usados:

```bash
# PC
PROCS="2 4 8" N_STRONG="1000 2000 4000" REPS=5 ./scripts/run_experiments.sh strong
./scripts/run_experiments.sh precision

# Cluster
PROCS="2 4 8 16 32" N_STRONG="1000 2000 4000 8000" REPS=5 ./scripts/run_experiments.sh strong
```

---

## Dónde están los archivos

| Archivo | Contenido |
|---|---|
| `results/<sistema>/raw.csv` | Todas las corridas, sin procesar |
| `results/<sistema>/resumen.md` | Medianas, speedup, eficiencia y GFLOP/s por punto |
| `results/<sistema>/info_<sistema>.txt` | Características de la máquina |
| `figuras/<sistema>/fig1_speedup.png` | Speedup según p, una curva por N |
| `figuras/<sistema>/fig2_eficiencia.png` | Eficiencia según p |
| `figuras/<sistema>/fig3_descomposicion.png` | Cómputo vs comunicación, para el N que elige el script |
| `figuras/pc/fig6_residuo.png` | Residuo relativo según N, por tipo de matriz |

---

## Speedup (mediana)

### PC

| N | T₁ serial | p = 2 | p = 4 | p = 8 |
|---|---|---|---|---|
| 1000 | 0.11 s | 1.45 | 1.81 | 2.44 |
| 2000 | 1.83 s | 1.49 | 2.10 | 2.04 |
| 4000 | 18.0 s | 1.41 | 1.77 | 1.88 |

### Cluster

| N | T₁ serial | p = 2 | p = 4 | p = 8 | p = 16 | p = 32 |
|---|---|---|---|---|---|---|
| 1000 | 0.07 s | 1.78 | 2.61 | 3.40 | 2.94 | 2.71 |
| 2000 | 1.13 s | 1.53 | 5.76 | 8.75 | 11.26 | 11.90 |
| 4000 | 10.4 s | 1.40 | 4.33 | 10.87 | 15.72 | 24.53 |
| 8000 | 82.8 s | 1.37 | 3.74 | 5.68 | 10.74 | 34.15 |

---

## Observaciones

**PC: el límite es la memoria, no la comunicación.** El speedup se estanca
cerca de 2. Con N = 4000 y p = 2 la comunicación es solo ~3 % del tiempo, y
aun así el speedup es 1.4. El serial rinde ~6 GFLOP/s con N = 1000 (la matriz,
de 8 MB, entra en la L3) y ~2 GFLOP/s con N = 4000 (128 MB, no entra). Los 8
procesos comparten una sola L3 y el mismo bus de memoria, así que agregar
procesos no agrega ancho de banda.

**Cluster: speedup superlineal cuando el bloque de cada proceso entra en la
cache.** Hay casos con eficiencia mayor al 100 %:

| Caso | Datos por proceso | Eficiencia |
|---|---|---|
| N = 2000, p = 4 | 500 filas × 2000 × 8 B = 8 MB | 144 % |
| N = 4000, p = 8 | 500 filas × 4000 × 8 B = 16 MB | 136 % |
| N = 8000, p = 32 | 250 filas × 8000 × 8 B = 16 MB | 107 % |

En los tres casos, las filas de cada proceso entran en una L3 de 16 MiB,
mientras que el serial no. Cada proceso trabaja desde la cache en lugar de la
RAM, y por eso se gana más que p veces. Con p = 2 los bloques no entran, y el
speedup queda cerca de 1.4 para todos los N, igual que en la PC.

**Con N chico domina la comunicación.** Con N = 1000 el speedup del cluster
llega a 3.4 con p = 8 y después baja: hay 2N colectivas por factorización y
no se achican con p, mientras que el cómputo por proceso sí. En la PC, con
N = 1000 y p = 8, la comunicación es el 57 % del tiempo.

**Precisión.** El residuo relativo queda entre 1e-16 y 5e-15 para las matrices
bien condicionadas (diagonal dominante, aleatoria y diagonal nula). La matriz de
diagonal nula se resuelve gracias al pivoteo; sin pivoteo falla en el paso 0.
Hilbert da un residuo ~1e-18, pero la solución calculada no sirve
(κ(A) ≈ 1e20 y elementos de x del orden de 1e18): un residuo chico no
garantiza una solución exacta si el problema está mal condicionado.

---

## Notas sobre las figuras

- **fig1:** la línea de speedup ideal se ve quebrada porque el eje x está en
  escala logarítmica y el eje y no. Los datos son correctos.
- **fig3:** el script elige automáticamente N = 1000, que es el caso donde más
  pesa la comunicación. Como no se usó `--split-idle`, la comunicación incluye
  también el tiempo de espera por desbalance, y la barra de "ocioso" queda vacía.
- **fig5 (escalabilidad débil):** no se corrió ese experimento. El script la
  genera igual con un solo punto tomado de escalabilidad fuerte, así que se
  descartó.

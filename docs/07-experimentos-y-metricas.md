# 07 · Experimentos y métricas

> Prerrequisitos: [03](03-costo-y-por-que-paralelizar.md),
> [05](05-lu-en-paralelo.md), [06](06-el-codigo-explicado.md).

---

## 1. Las métricas

### Speedup

```
           T₁
   S(p) = ─────
           T(p)
```

`T₁` es el tiempo de la **versión serial pura** (`bin/lu_serial`), **no** el
de la versión MPI corrida con un proceso.

> Esto es lo primero que se revisa en un informe de HPC. Usar la versión MPI
> con `p = 1` como `T₁` infla el speedup porque arrastra el mismo overhead
> arriba y abajo de la fracción: estás comparando el programa contra sí mismo,
> no contra el algoritmo secuencial.

El ideal es `S(p) = p` (speedup lineal). Lo normal es sublineal.

### Eficiencia

```
           S(p)
   E(p) = ──────
             p
```

Es el número honesto. `S(16) = 6` suena bien hasta que lo traducís: es
**37 % de eficiencia**, o sea que 10 de los 16 procesadores estuvieron
trabajando para nada.

### GFLOP/s

```
             (2/3)·n³
   GFLOP/s = ──────────── / 10⁹
                 T
```

Se usa el **conteo de operaciones del algoritmo** `(2/3)n³`, no la cantidad de
instrucciones que realmente ejecutó la implementación. Es la convención
estándar (la misma que usa HPL / TOP500) y es lo que hace comparables
implementaciones distintas.

### Residuo relativo

```
        ‖A·x − b‖∞
   r = ──────────────
        ‖A‖∞ · ‖x‖∞
```

Esperable: `O(ε) ≈ 10⁻¹⁶`. Aceptable hasta `~10⁻¹³`. Por encima de eso, hay
que explicar por qué (condicionamiento) o buscar el bug.

### Escalabilidad fuerte vs. débil

| | Escalabilidad **fuerte** | Escalabilidad **débil** |
|---|---|---|
| Qué se fija | el **problema** (`n` constante) | el **trabajo por proceso** |
| Qué varía | `p` | `p` **y** `n` juntos |
| Pregunta que responde | "¿resuelvo *este* problema más rápido?" | "¿puedo resolver problemas *más grandes*?" |
| Métrica | speedup, eficiencia | tiempo, que debería quedar plano |
| Límite | Amdahl | Gustafson |

**La trampa de la escalabilidad débil en LU:** como el trabajo crece como
`n³`, para mantener el trabajo por proceso constante hay que tomar

```
   n(p) = n₀ · p^(1/3)
```

y **no** `n(p) = n₀·p`, que es el error habitual (ese multiplicaría el trabajo
por proceso por `p²` y daría una curva desastrosa sin que haya nada mal).

| p | `n = 1000·p^(1/3)` | trabajo total | trabajo por proceso |
|---|---|---|---|
| 1 | 1 000 | 6.7 × 10⁸ | 6.7 × 10⁸ |
| 2 | 1 260 | 1.3 × 10⁹ | 6.7 × 10⁸ |
| 4 | 1 587 | 2.7 × 10⁹ | 6.7 × 10⁸ |
| 8 | 2 000 | 5.3 × 10⁹ | 6.7 × 10⁸ |
| 16 | 2 520 | 1.1 × 10¹⁰ | 6.7 × 10⁸ |

---

## 2. El protocolo de medición

Esto no es burocracia: es lo que hace que los resultados signifiquen algo.

### (a) Mínimo 5 repeticiones, y se reporta la MEDIANA

El cluster no es dedicado. Hay otros jobs, ruido del sistema operativo,
contención de red. Eso produce corridas ocasionalmente **muy lentas**, pero
nunca corridas anormalmente rápidas: la distribución de tiempos tiene una cola
larga hacia arriba.

El **promedio** se deja arrastrar por esa cola. La **mediana** no.

Reportamos además el mínimo y el máximo observados, para que se vea cuánto
ruido había. Si el máximo es 3× la mediana, eso es información: el cluster
estaba cargado.

> La consigna del TP1 lo pide explícitamente: *"No haga una única ejecución
> para calcular el speed up. Corra el código múltiples veces y saque
> estadísticas (como la mediana)."*

### (b) Qué entra y qué no entra en el cronómetro

| Dentro | Fuera |
|---|---|
| El lazo de factorización | Generación de la matriz |
| | Alocación de memoria |
| | Copia de restauración entre repeticiones |
| | Sustituciones triangulares |
| | Cálculo del residuo |
| | Escritura de resultados |

Y siempre un `MPI_Barrier` antes de arrancar, para que todos larguen parejo.

### (c) Mismas banderas de compilación en los dos binarios

Si el serial se compila con `-O3` y el MPI con `-O0`, el speedup que midan es
ficción. El `Makefile` usa la misma variable `CFLAGS` para los dos.

### (d) Compilar en el mismo tipo de nodo donde se va a correr

`-march=native` genera instrucciones para la CPU **del nodo que compila**. Si
el nodo de login tiene otra CPU que los de cómputo, el binario puede fallar
con *"Illegal instruction"*. Por eso `job.slurm` compila dentro del job. Si
hay problemas: `make ARCH=`.

---

## 3. La batería de experimentos

Todo está en `scripts/run_experiments.sh`.

| # | Experimento | N | p | Responde |
|---|---|---|---|---|
| **E1** | Validación | 100, 500 | 1…8 | ¿El paralelo da lo mismo que el serial? |
| **E2** | Escalabilidad fuerte | 1000, 2000, 4000 | 1…32 | ¿Cómo escala el cómputo? |
| **E3** | Escalabilidad débil | `1000·p^(1/3)` | 1…16 | ¿Sirve para problemas más grandes? |
| **E4** | Descomposición del tiempo | 4000 | 1…32 | ¿Cómo escala la comunicación? |
| **E5** | Cíclico vs bloques | 2000 | 2…32 | ¿Cuánto importa el reparto? |
| **E6** | Precisión | 100…1000 | 4 | ¿Cuál es la precisión del residuo? |
| **E7** | Contra LAPACK *(opcional)* | 4000 | 1 | ¿Cuánto nos falta para el estado del arte? |

### Cómo correrlos

```bash
# todo (en el cluster, desde un job):
./scripts/run_experiments.sh all

# uno solo:
./scripts/run_experiments.sh strong

# ajustando parámetros:
PROCS="1 2 4 8 16 32" REPS=7 ./scripts/run_experiments.sh strong
```

Los resultados se **acumulan** en `results/raw.csv` (no se pisan), con todas
las repeticiones crudas. La estadística la hace `scripts/plots.py`, así que
siempre se puede volver sobre los datos originales.

---

## 4. Las figuras del informe

`python3 scripts/plots.py` genera seis figuras en `figuras/`:

| Figura | Qué muestra | Qué hay que escribir debajo |
|---|---|---|
| `fig1_speedup.png` | Speedup vs p, una curva por N | Dónde se aparta de la recta ideal y por qué; cómo mejora al crecer N |
| `fig2_eficiencia.png` | Eficiencia vs p | El número honesto. A partir de qué p deja de valer la pena |
| `fig3_descomposicion.png` | Barras apiladas: cómputo / comunicación / ocioso | Cómo crece la fracción de comunicación con p |
| `fig4_reparto.png` | Cíclico vs bloques | La evidencia experimental del argumento de balance de carga |
| `fig5_debil.png` | Escalabilidad débil | Qué tan plana queda la curva y qué significa |
| `fig6_residuo.png` | Residuo vs N por tipo de matriz | Que el residuo se mantiene en ε; y que eso no implica solución exacta |

Además genera `results/resumen.md` con la tabla de speedup y eficiencia, lista
para pegar en el informe.

**Cada figura necesita un párrafo que la explique.** Una gráfica sin
interpretación no suma nada: la consigna de la cátedra pide explícitamente
*"generar resultados y explicarlos"*.

---

## 5. Cómo interpretar lo que van a ver

### El speedup se aplana

**Esperable.** Tres causas simultáneas:

1. **Amdahl**: la parte serial (búsqueda del pivote, difusión) no se paraleliza.
2. **Comunicación**: `O(n)` colectivas de `O(log p)` etapas, y el costo crece
   con `p`.
3. **Desbalance**: en el paso `k` solo trabajan los procesos con filas `> k`.

La figura 3 cuantifica cuánto pesa cada una. Eso es lo que convierte una
observación en un análisis.

### El speedup baja después de cierto `p`

**También esperable**, y es un resultado, no un fracaso. Ocurre cuando `n/p`
se hace chico y la comunicación domina. Identificar ese punto y explicarlo con
el modelo de costos es uno de los mejores resultados que puede tener el
informe.

### Speedup superlineal (`S(p) > p`)

Puede pasar, y hay que explicarlo, no esconderlo. Suele ser **efecto de
cache**: al repartir la matriz, la porción de cada proceso entra en niveles de
cache más rápidos que la matriz completa en la corrida serial. No viola nada:
es que `T₁` está medido con una jerarquía de memoria peor.

### Curvas de speedup que se separan por N

**Es el resultado central.** La curva de `N = 4000` va a estar por encima de la
de `N = 1000` para todo `p > 1`. Esa separación **es** la verificación
experimental de que `T_com/T_comp ~ p/n`.

### El residuo no cambia con p

Correcto y esperable. Confirma que la paralelización no introdujo error
numérico. Si cambiara, habría un bug.

---

## 6. Errores frecuentes

| Error | Consecuencia |
|---|---|
| `T₁` = versión MPI con p=1 | Speedup inflado. Lo detectan enseguida |
| Una sola corrida por punto | Las curvas quedan con dientes de sierra sin sentido |
| Promedio en vez de mediana | Una corrida contaminada arruina el punto |
| Escalabilidad débil con `n = n₀·p` | Curva desastrosa que no significa nada |
| Medir con `n = 500` y `p = 32` | 15 filas por proceso: se mide latencia, no el algoritmo |
| Banderas de compilación distintas | El speedup no significa nada |
| Incluir la generación de la matriz | Distorsión creciente con `n` |
| Usar corridas con `--split-idle` para el speedup | Las barreras agregan overhead |
| Gráficas sin párrafo explicativo | La consigna pide explicarlas |

---

## 7. Preguntas típicas de defensa

1. ¿Por qué `T₁` es el serial puro y no el MPI con p=1?
2. ¿Por qué la mediana y no el promedio?
3. Diferencia entre escalabilidad fuerte y débil.
4. ¿Por qué en escalabilidad débil `n` crece como `p^(1/3)`?
5. ¿Por qué el speedup se aplana? Nombrá las tres causas.
6. ¿Puede haber speedup superlineal? ¿Cómo lo explicarías?
7. ¿Qué entra y qué no entra en el tiempo medido?
8. Si la eficiencia con 16 procesos es 40 %, ¿qué significa eso en recursos?

---

**Siguiente:** [08 · Guía del cluster](08-guia-del-cluster.md)

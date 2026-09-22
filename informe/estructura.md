# Estructura del informe

Formato tomado del documento `Formato de informe.docx` del material de la
cátedra. **Confirmar con la docente que sigue vigente para este proyecto antes
de maquetar.**

---

## Portada

- Nombre de la materia: Métodos Numéricos II y Computación Científica
- Carrera
- Integrantes (los tres)
- Título del trabajo: *Factorización LU paralela con MPI: análisis de
  escalabilidad y precisión numérica*
- Año: 2026

---

## Resumen — hasta 200 palabras

Tiene que decir de qué trata **y adelantar los resultados**. Se escribe último.

Esqueleto, para completar con los números reales:

> Se implementó la factorización LU con pivoteo parcial de matrices densas en
> paralelo, utilizando el paradigma de pasaje de mensajes (MPI) sobre un modelo
> de memoria distribuida. La matriz se reparte entre los procesos mediante
> distribución cíclica por filas, elegida por su balance de carga sostenido
> frente a la distribución en bloques contiguos. Se evaluó la escalabilidad
> fuerte y débil sobre \<hardware\> para matrices de hasta N = \<N\>,
> alcanzando un speedup de \<S\> con \<p\> procesos (\<E\> % de eficiencia).
> La descomposición del tiempo de ejecución muestra que la fracción dedicada a
> comunicación crece de \<x\> % con 2 procesos a \<y\> % con \<p\>, lo que
> explica el aplanamiento de la curva de speedup y es consistente con el modelo
> teórico T_com/T_comp ~ p/N. La verificación numérica arrojó un residuo
> relativo ‖Ax−b‖∞/(‖A‖∞‖x‖∞) del orden de \<r\>, idéntico al de la
> implementación serial de referencia para todo número de procesos, lo que
> confirma que la paralelización no introduce error adicional.

---

## Introducción — hasta 1 carilla

**Qué va:**

1. **El problema.** Qué es la factorización LU, por qué `Ax = b` denso aparece
   en todos lados, por qué se factoriza en vez de invertir. (→ [doc 01](../docs/01-que-es-la-factorizacion-lu.md))
2. **Por qué merece HPC.** El costo `(2/3)n³`, la tabla de tiempos y memoria
   vs `n`. (→ [doc 03](../docs/03-costo-y-por-que-paralelizar.md))
3. **Estado del arte.** ScaLAPACK y la distribución 2D block-cyclic; los
   algoritmos por bloques y BLAS3; CALU y el pivoteo por torneo; SLATE y las
   variantes para GPU del Exascale Computing Project; HPL como benchmark del
   TOP500. (→ [doc 05 §7](../docs/05-lu-en-paralelo.md))
4. **Objetivos concretos del trabajo.** Enumerados, para poder responderlos
   uno por uno en las conclusiones:
   - Implementar LU con pivoteo parcial en serial y en paralelo con MPI.
   - Verificar la corrección numérica mediante el residuo `‖Ax−b‖`.
   - Medir la escalabilidad fuerte y débil, y el speedup alcanzado.
   - Cuantificar cómo escalan por separado el cómputo y la comunicación al
     incrementar N.
   - Comparar experimentalmente dos estrategias de distribución de datos.

---

## Datos y Métodos — máximo 1½ carillas

**Qué va:**

1. **Matrices utilizadas.** Los cuatro tipos (`diagdom`, `random`, `hilbert`,
   `zerodiag`), cómo se generan y para qué sirve cada uno. Mencionar el
   generador posicional sin estado y **por qué** hace falta (permite que cada
   proceso genere sus filas sin que nadie aloque `n²`, y garantiza que serial y
   paralelo resuelvan el mismo sistema). (→ [doc 06 §2](../docs/06-el-codigo-explicado.md))
2. **El algoritmo.** Right-looking con pivoteo parcial, almacenamiento
   in-place, representación de `P` como vector de enteros.
3. **La paralelización.** **Acá va la figura de distribución de datos**
   (cíclica vs bloques): es la figura más importante del informe. Las cuatro
   operaciones del paso `k` y qué colectiva MPI usa cada una.
4. **El modelo de costos.** `T_comp ≈ (2/3)n³/p`, `T_com ≈ n·log(p)·α +
   (n²/2)·log(p)·β`, y el cociente `~ p/n`.
5. **Hardware y software.** Modelo de CPU, cores por nodo, RAM, red,
   compilador y banderas, versión de MPI, gestor de colas. (→ [doc 08 §6](../docs/08-guia-del-cluster.md))
6. **Protocolo de medición.** 5+ repeticiones, mediana, qué entra y qué no
   entra en el cronómetro, `T₁` = serial puro, barrera antes de arrancar.

---

## Resultados — máximo 2 carillas

Las figuras con su interpretación. **Cada una necesita un párrafo.**

| Figura | Qué escribir |
|---|---|
| `fig1_speedup.png` | Dónde se aparta de la recta ideal. Cómo la curva mejora al crecer N. El p a partir del cual deja de valer la pena |
| `fig2_eficiencia.png` | El número honesto. Traducirlo a recursos |
| `fig3_descomposicion.png` | Cómo crece la fracción de comunicación con p. Conectarlo con el aplanamiento de fig1 |
| `fig4_reparto.png` | La evidencia experimental del argumento de balance de carga |
| `fig5_debil.png` | Qué tan plana queda la curva. Qué dice sobre resolver problemas más grandes |
| `fig6_residuo.png` | El residuo se mantiene en ε. Y la distinción residuo / error con Hilbert |

Más la tabla de `results/resumen.md`.

**Obligatorio:** el link al repositorio con el código y los datos:

> El código fuente, los scripts de ejecución y los datos crudos de todas las
> corridas están disponibles en https://github.com/mainaessens/mnii-hpc-lu

---

## Conclusiones — 1 carilla

**Responder objetivo por objetivo**, en el mismo orden en que se enunciaron en
la introducción. Después:

- **El límite encontrado y por qué.** A partir de qué `n/p` deja de escalar, y
  qué lo causa según la descomposición del tiempo.
- **Qué se confirmó del modelo teórico.** La separación de las curvas de
  speedup por N es la verificación experimental de `T_com/T_comp ~ p/n`.
- **Trabajo futuro.** Distribución 2D block-cyclic; algoritmo por bloques con
  BLAS3; pivoteo por torneo (CALU). Con una frase de por qué cada uno
  ayudaría.

---

## Referencias

Normas APA. La lista está en [`referencias.md`](referencias.md).

---

## Checklist antes de entregar

- [ ] El resumen adelanta resultados concretos, con números
- [ ] Todos los objetivos de la introducción están respondidos en las conclusiones
- [ ] Cada figura tiene su párrafo de interpretación
- [ ] La descripción del hardware está completa (sin corchetes sin llenar)
- [ ] El link al repositorio está y funciona
- [ ] Las referencias están en APA y todas se citan en el texto
- [ ] Se respetan los límites de extensión de cada sección
- [ ] Las tres personas leyeron el informe entero

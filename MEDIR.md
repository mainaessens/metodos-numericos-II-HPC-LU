# Cómo correr las mediciones (PC y cluster)

Guía para repetir las corridas en otra máquina con **los mismos parámetros que
las de Ezequiel** (ver `RESULTS_EM.md`), más el experimento E4 con
`--split-idle`, que faltaba.

## Experimentos acordados

| # | Experimento | ¿Se hace? | PC | Cluster |
|---|---|---|---|---|
| E1 | Validación | sí (`medir.sh` chequea el residuo antes de empezar) | ✔ | ✔ |
| E2 | Escalabilidad fuerte | sí | N = 1000, 2000, 4000 · p = 2, 4, 8 | N = 1000…8000 · p = 2…32 |
| E3 | Escalabilidad débil | **no** | | |
| E4 | Descomposición del tiempo (`--split-idle`) | sí | N = 4000 · p = 2, 4, 8 | N = 4000 · p = 2…32 |
| E5 | Cíclico vs bloques | **no** | | |
| E6 | Precisión | sí | ✔ | ✔ (sirve para mostrar que no depende del hardware) |
| E7 | Contra LAPACK | **no** | | |

Todo con matriz diagonal dominante, reparto cíclico, semilla 42, 5
repeticiones, mediana. T₁ = serial puro.

Todo lo anterior lo hace un solo comando:

```bash
./scripts/medir.sh <sistema> <pc|cluster>
```

y deja `results/<sistema>/` (raw.csv, info de la máquina, log, resumen) y
`figuras/<sistema>/`. Las corridas con `--split-idle` **no** entran en el
speedup: `plots.py` las separa solo (las reconoce porque tienen `t_idle > 0`).

---

## En la PC (Windows → WSL)

1. Abrir **Ubuntu (WSL)**. Trabajar en el home de Linux, **no** en
   `/mnt/c/...` ni en OneDrive.
2. Una sola vez, instalar lo necesario:
   ```bash
   sudo apt update
   sudo apt install -y build-essential openmpi-bin libopenmpi-dev python3-pandas python3-matplotlib git
   ```
3. Clonar (o actualizar) el repo:
   ```bash
   cd ~
   git clone https://github.com/mainaessens/metodos-numericos-ii-hpc-lu.git
   cd metodos-numericos-ii-hpc-lu
   # si ya estaba clonado:  git pull origin main
   ```
4. Antes de medir: notebook **enchufada**, modo de energía **"Máximo
   rendimiento"**, cerrar navegador/Teams/etc. No usar la PC mientras corre.
5. Correr (tarda ~15 min):
   ```bash
   chmod +x scripts/*.sh
   ./scripts/medir.sh pc-maia pc
   ```

## En el cluster (lab.ccad.unc.edu.ar)

1. Entrar a https://lab.ccad.unc.edu.ar/ con la cuenta de Google.
2. Entorno de trabajo: **Python Base**. Recursos: **48 núcleos** (los mismos
   que usó Ezequiel; como mínimo 32).
3. En VS Code web, abrir una terminal y clonar:
   ```bash
   cd ~
   git clone https://github.com/mainaessens/metodos-numericos-ii-hpc-lu.git
   cd metodos-numericos-ii-hpc-lu
   # si ya estaba clonado:  git pull origin main
   chmod +x scripts/*.sh
   ```
4. Correr en segundo plano para que no se corte si se cierra la pestaña
   (tarda ~25 min, sobre todo por N = 8000):
   ```bash
   nohup ./scripts/medir.sh cluster-maia cluster > medir.out 2>&1 &
   tail -f medir.out          # Ctrl+C sale del tail, NO corta la corrida
   ```
5. Al terminar, aparece `Fin:` en `medir.out`.

## Subir los resultados

**Opción fácil (PC):** copiar las carpetas al repo de Windows y que lo suba
otro (o Claude) desde ahí:

```bash
DEST="/mnt/c/Users/naess/OneDrive/Documentos/5to año/Métodos numéricos II/mnii-hpc-lu"
mkdir -p "$DEST/results" "$DEST/figuras"
cp -r results/pc-maia "$DEST/results/"
cp -r figuras/pc-maia "$DEST/figuras/"
```

**Opción git:** desde la máquina donde se corrió (PC o cluster):

```bash
git checkout -b results-maia          # la primera vez; después: git checkout results-maia
git add results/pc-maia figuras/pc-maia          # o cluster-maia
git commit -m "resultados de Maia: PC (o cluster)"
git push origin results-maia
```

Después se mergea `results-maia` a `main`, igual que se hizo con
`results-ezequiel`.

Si en el cluster `git push` pide credenciales y no las tienen ahí, alternativa:
descargar las carpetas `results/cluster-maia` y `figuras/cluster-maia` desde
el explorador de VS Code web (clic derecho → *Download*) y subirlas desde la PC.

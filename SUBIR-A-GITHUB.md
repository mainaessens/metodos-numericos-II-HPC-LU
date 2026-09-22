# Cómo subir esto a GitHub

Todo el proyecto está acá, listo para hacer el primer push. Son tres pasos.

---

## 1. Crear el repositorio vacío

1. Entrá a **https://github.com/new**
2. **Repository name:** `mnii-hpc-lu`
3. **Public**
4. **NO** marques "Add a README file", ni .gitignore, ni licencia — tiene que
   quedar completamente vacío
5. **Create repository**

---

## 2. Subirlo

Abrí una terminal (Git Bash, PowerShell o la terminal de VS Code) **dentro de
esta carpeta** y pegá:

```bash
git init -b main
git add .
git commit -m "Estructura inicial: LU serial y MPI, scripts, guia de estudio"
git remote add origin https://github.com/mainaessens/mnii-hpc-lu.git
git push -u origin main
```

Si te pide credenciales, GitHub ya no acepta contraseña: usá un **Personal
Access Token** (Settings → Developer settings → Personal access tokens →
Tokens (classic) → Generate new token, con el permiso `repo`). O instalá
[GitHub CLI](https://cli.github.com/) y corré `gh auth login` una vez.

---

## 3. Agregar a las otras dos personas

En el repo → **Settings → Collaborators → Add people**.

Con el repo público igual pueden clonarlo sin ser colaboradores, pero para
hacer `push` necesitan el permiso.

---

## Flujo de trabajo del equipo

Para evitar pisarse:

```bash
# antes de empezar a trabajar, siempre:
git pull

# al terminar algo:
git add .
git commit -m "que hiciste"
git push
```

Si van a tocar el mismo archivo al mismo tiempo, conviene trabajar en ramas:

```bash
git checkout -b mi-rama
# ... trabajar, commitear ...
git push -u origin mi-rama
# y después abrir un Pull Request en GitHub
```

---

## Qué NO se sube (está en `.gitignore`)

- `bin/`, `build/`, `*.o` — binarios compilados, se regeneran con `make`
- `logs/` — salidas de los jobs del cluster

## Qué SÍ se sube, aunque parezca que no

- **`results/*.csv`** — los datos crudos de todas las corridas. La consigna
  pide entregar el repo "con el código **y los datos usados**"
- **`figuras/*.png`** — las figuras que van al informe

---

## Lo primero que hay que completar

- `README.md` → la sección **Equipo**, con los tres nombres
- `docs/08-guia-del-cluster.md` → los datos reales del cluster apenas los tengan
- `scripts/job.slurm` → partición, nodos y tiempo límite reales

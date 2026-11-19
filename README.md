<div align="center">

<img width="200" alt="PAIPAI-logo" src="logo.png"/>

# Package for Alloy Interstitial Predictions using Artificial Intelligence<br>(PAIPAI)

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://spdx.org/licenses/MIT.html)
![C++](https://img.shields.io/badge/C%2B%2B-17-blue)
![Python](https://img.shields.io/badge/python-3.12-blue)
![Platforms](https://img.shields.io/badge/platform-linux%20%7C%20macos-lightgrey)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21480012.svg)](https://doi.org/10.5281/zenodo.21480012)
[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.commatsci.2026.114752-blue)](https://doi.org/10.1016/j.commatsci.2026.114752)

`PAIPAI` searches for low-energy configurations of multi-component alloys with interstitials and defects by combining Monte Carlo sampling with machine learning interatomic potentials (MLIPs). This repository is a fork of [siyazhu/PAIPAI](https://github.com/siyazhu/PAIPAI) v2.0 by Siya Zhu with bug fixes and two new options, `--seed` and `--pack-every`.

<p>
  <a href="https://github.com/dogusariturk/PAIPAI/issues/new?labels=bug">Report a Bug</a> |
  <a href="https://github.com/dogusariturk/PAIPAI/issues/new?labels=enhancement">Request a Feature</a>
</p>

</div>

---

## Key Features

- Search for low-energy structures with two kinds of workers: fast workers screen trial structures with short relaxations, and slow workers refine the promising ones (`search` mode)
- Sample finite-temperature ensembles with sequential Markov-chain Monte Carlo (`finiteT` mode), and repeat a run with `--seed`
- Propose metal swaps, interstitial swaps, cluster interstitial moves, and local interstitial hops, with adjustable weights
- Handle bulk, grain-boundary, and slab structures with interstitial solutes such as B, C, N, and O
- Move the interstitial sites with the cell when a relaxation changes the cell shape
- Relax and evaluate structures with MLIPs through [MaterialsFramework](https://github.com/dogusariturk/MaterialsFramework) on CPU or GPU
- Start from an earlier run's state (`--resume-state`) or continue an interrupted run (`--continue`); PAIPAI checks the state against the run's cutoffs and requeues tasks the interrupted run left unfinished
- Pack accepted states into `tar.gz` archives during a run (`--pack-every`)
- Stop the whole run as soon as any worker exits
- Build `struc.in` input files from a metal-only POSCAR with `findinter`
- Post-process accepted states with `warrencowley` (Warren-Cowley short-range order) and `packmc` (CONTCAR archive and summary table)

---

## Installation

### Prerequisites

- CMake ≥ 3.15
- C++17 compiler (GCC or Clang)
- Python 3.12 and [uv](https://docs.astral.sh/uv/)
- [MaterialsFramework](https://github.com/dogusariturk/MaterialsFramework) with an MLIP extra (GRACE for the default model)
- (Optional) NVIDIA GPU with CUDA for `--device cuda`

### Build Instructions

**1. Clone the repository**

```bash
git clone https://github.com/dogusariturk/PAIPAI.git
cd PAIPAI
```

**2. Create the Python environment for the workers**

```bash
uv venv --python 3.12
uv pip install "materialsframework[grace]"
```

**3. Build and install**

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$HOME/.local" \
  -DPython3_EXECUTABLE="$PWD/.venv/bin/python"
cmake --build build -j 8
cmake --install build
paipai --help
```

> [!NOTE]
> This installs `paipai`, `mc_paipai`, `findinter`, `warrencowley`, `packmc`, and the worker scripts into `~/.local/bin`. Add `~/.local/bin` to your `PATH` if it isn't there already. The `paipai` wrapper runs the workers with `.venv/bin/python`; set `PAIPAI_PYTHON` to use a different Python.

---

## Quickstart

### Input structure

```bash
findinter --input POSCAR --inter B --internum 8 --gn 100 \
  --output struc.in --site-poscar interstitial_sites.vasp
```

> [!NOTE]
> On its first run in a directory, `findinter` writes a default `radii.dat` and exits. Check the radii, then run the same command again. See [docs/findinter.md](docs/findinter.md) for the site-search parameters.

### Ground-state search

```bash
paipai --input struc.in --mode search --device cpu \
  --fast 15 --slow 15 --steps 2000 --temp 10 --p-cluster-inter 0
```

A low Monte Carlo temperature works best for structure search.

### Finite-temperature sampling

```bash
paipai --input struc.in --mode finiteT --device cuda --ngpu 1 \
  --fast 0 --slow 1 --steps 200000 --temp 700 --p-cluster-inter 0 \
  --seed 1 --pack-every 60
```

To start a chain from the end of a search, replace `--input struc.in` with `--resume-state <search run>/mcprocess`, which picks up the latest accepted state. To continue a stopped run, rerun it in the same `--root` with `--continue` in place of `--input` or `--resume-state`.

### Analysis

```bash
warrencowley --inter-metal-cutoff 2.5 --metal-metal-cutoff 3.0
packmc
```

> [!NOTE]
> Both tools read `mcprocess/`. If the run used `--pack-every`, unpack the archives first:
>
> ```bash
> for f in mcprocess_packed/*.tar.gz; do tar -xzf "$f" -C mcprocess; done
> ```

### More options

`paipai --help` lists every option. Design notes are in [docs/](docs/), and the [upstream README](https://github.com/siyazhu/PAIPAI/blob/v2.0-dev/README.md) has the full option tables and a description of each move type.

---

## License

Distributed under the MIT License. See [LICENSE](LICENSE) for details. `include/json.hpp` is third-party code with its own license notice.

---

## Citation

PAIPAI was written by Siya Zhu. If you use it in your research, please cite:

> Zhu, S., & Arróyave, R. (2026). Ground-state structure search of defective high-entropy alloys using machine-learning potentials and Monte Carlo sampling. *Computational Materials Science*, 270, 114752. https://doi.org/10.1016/j.commatsci.2026.114752

> Zhu, S. (2026). siyazhu/PAIPAI: PAIPAI V2.0 (v2.0). Zenodo. https://doi.org/10.5281/zenodo.21480012

```bibtex
@article{zhu_2026_paipai,
  author  = {Zhu, Siya and Arróyave, Raymundo},
  title   = {Ground-state structure search of defective high-entropy alloys using machine-learning potentials and Monte Carlo sampling},
  journal = {Computational Materials Science},
  volume  = {270},
  pages   = {114752},
  year    = {2026},
  doi     = {10.1016/j.commatsci.2026.114752},
}

@software{zhu_2026_21480012,
  author    = {Zhu, Siya},
  title     = {siyazhu/PAIPAI: PAIPAI V2.0},
  version   = {v2.0},
  month     = jul,
  year      = 2026,
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.21480012},
  url       = {https://doi.org/10.5281/zenodo.21480012},
}
```

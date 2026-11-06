# Intermodulation Calculation
This script calculates the critical intermodulation frequencies in a given list of frequencies. See [this page](https://www.shure.com/en-US/performance-production/louder/all-about-wireless-intermodulation-distortion) for a good introduction on the topic

The default list is not particularly well optimized.

## Finding an optimized frequency set
`optimize.py` searches for a frequency set that is as intermodulation free as possible.
No dependencies beyond the Python standard library (>= 3.9) are needed at runtime.
Build the compiled search core once (needs a C compiler, e.g. gcc or MSVC Build Tools, and setuptools):

```bash
python setup.py build_ext --inplace
```

Alternatively, `pip install .` installs it including the command `im-optimize` (same arguments as `optimize.py`).

```bash
python optimize.py --start 400.0 --stop 406.0 --spacing 0.3 --count 16
```

Optional: `--threshold` (default 0.035 MHz), `--step` grid (default 0.01 MHz), `--time-limit` of the exact search (default 60 s, negative = unlimited), `--no-exact`, `--iterations`, `--restarts`, `--seed`.

How it works:
- IM products (same IM3/IM5/IM7 definition as `main.py`) within the threshold of a carrier are penalized, weighted IM3 (1000) >> IM5 (30) >> IM7 (1).
- Simulated annealing finds a good set quickly, its cost is the upper bound for an exact branch & bound search (`imcore.c`).
- The result is reported as *proven optimal* if the exact search finishes within the time limit (or a set without any critical IM product is found). Otherwise it is the best set found.
- Without the compiled core, only the annealing result is returned.

Examples for 400–406 MHz:

| Carriers | Spacing | Default list in `main.py` | Optimized | Proven optimal |
|---|---|---|---|---|
| 16 | 0.3 MHz | 10 IM5, 7 IM7 | IM free | yes (< 1 s) |
| 17 | 0.3 MHz | – | 2 IM7 | no (60 s) |
| 20 | 0.2 MHz | – | 1 IM5, 10 IM7 | no (30 s) |

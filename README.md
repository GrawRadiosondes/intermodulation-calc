# Intermodulation Calculation
This script calculates the critical intermodulation frequencies in a given list of frequencies. See [this page](https://www.shure.com/en-US/performance-production/louder/all-about-wireless-intermodulation-distortion) for a good introduction on the topic

The default list is not particularly well optimized.

## Finding an optimized frequency set
`optimize.py` searches for a frequency set that is as intermodulation free as possible:

```bash
python optimize.py --start 400.0 --stop 406.0 --spacing 0.3 --count 16
```

Optional: `--threshold` (default 0.035 MHz), `--step` grid (default 0.01 MHz), `--iterations`, `--restarts`, `--seed`.

It uses simulated annealing on the frequency grid, penalizing IM products (same IM3/IM5/IM7 definition as `main.py`) within the threshold of a carrier, weighted IM3 >> IM5 >> IM7.
Example for the band above: the default list has 10 critical IM5 / 7 IM7 hits, the optimized set 0 IM5 / 7 IM7 (no IM3 in both).
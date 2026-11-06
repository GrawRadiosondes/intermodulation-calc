"""Find a frequency set that is as intermodulation free as possible.

Usage example:
    python optimize.py --start 400.0 --stop 406.0 --spacing 0.3 --count 16

All frequencies are placed on a grid (default 10 kHz) starting at --start. Because every
two-tone IM product k*f1 - (k-1)*f2 of grid frequencies lies on the same grid again, the whole
search works with integer grid indices, which makes it exact and fast.

Objective (lower is better): for every IM product (IM3/IM5/IM7, as in main.py) that lands within
--threshold of any carrier, a penalty weighted by its order is added. IM3 hits are far more expensive
than IM5 hits, which are more expensive than IM7 hits. A simulated annealing
(over "move one carrier to another free grid position" steps) finds a good set quickly. Its cost is
then the upper bound for an exact branch & bound search in the compiled core imcore.c, which either
improves the set or proves that no better set exists. Build the core once with:
    python setup.py build_ext --inplace
"""
import argparse
import math
import random
from collections import Counter

try:
    import imcore  # compiled exact search, build with: python setup.py build_ext --inplace
except ImportError:
    imcore = None

# (k, order name, penalty weight). Product = k*f1 - (k-1)*f2, see main.py
IM_TYPES = [(2, 'IM3', 1000), (3, 'IM5', 30), (4, 'IM7', 1)]


class IMState:
    """Carrier set on an integer grid with incrementally maintained IM product counts."""

    def __init__(self, channels, threshold):
        self.channels = list(channels)
        self.threshold = threshold
        self.products = {name: Counter() for _, name, _ in IM_TYPES}
        for i in range(len(self.channels)):
            self._update_products(i, +1, skip_lower=True)

    def _update_products(self, i, sign, skip_lower=False):
        """Add (sign=+1) or remove (sign=-1) all products in which channel i takes part."""
        a = self.channels[i]
        for j, b in enumerate(self.channels):
            if j == i or (skip_lower and j < i):
                continue
            for k, name, _ in IM_TYPES:
                counter = self.products[name]
                counter[k * a - (k - 1) * b] += sign
                counter[k * b - (k - 1) * a] += sign

    def hits(self, name):
        """Number of (product, carrier) pairs closer than the threshold for one IM order."""
        counter = self.products[name]
        t = self.threshold
        return sum(counter[c + d] for c in self.channels for d in range(-t, t + 1))

    def cost(self):
        return sum(weight * self.hits(name) for _, name, weight in IM_TYPES)

    def move(self, i, new_pos):
        self._update_products(i, -1)
        self.channels[i] = new_pos
        self._update_products(i, +1)


def is_free(pos, channels, skip, spacing):
    return all(abs(pos - c) >= spacing for j, c in enumerate(channels) if j != skip)


def random_start(n_grid, count, spacing, rng):
    """Random feasible placement: distribute the slack randomly between minimally spaced carriers."""
    slack = n_grid - 1 - (count - 1) * spacing
    cuts = sorted(rng.randint(0, slack) for _ in range(count))
    return [cuts[i] + i * spacing for i in range(count)]


def anneal(n_grid, count, spacing, threshold, iterations, rng):
    state = IMState(random_start(n_grid, count, spacing, rng), threshold)
    cost = state.cost()
    best_cost, best = cost, list(state.channels)
    t_start, t_end = 300.0, 0.3
    for it in range(iterations):
        temp = t_start * (t_end / t_start) ** (it / iterations)
        i = rng.randrange(count)
        old = state.channels[i]
        # mix of small local nudges and long jumps
        if rng.random() < 0.5:
            new = old + rng.randint(-3 * spacing, 3 * spacing)
        else:
            new = rng.randrange(n_grid)
        if new == old or not 0 <= new < n_grid or not is_free(new, state.channels, i, spacing):
            continue
        state.move(i, new)
        new_cost = state.cost()
        if new_cost <= cost or rng.random() < math.exp((cost - new_cost) / temp):
            cost = new_cost
            if cost < best_cost:
                best_cost, best = cost, list(state.channels)
                if best_cost == 0:
                    break
        else:
            state.move(i, old)
    return best_cost, sorted(best)


def optimize(start, stop, spacing, count, threshold=0.035, step=0.01,
             iterations=20000, restarts=3, seed=None, time_limit=60.0, exact=True):
    """Return (frequency list, hits per IM order, proven optimal) of the best set found. Frequencies in MHz.

    Simulated annealing provides a good set quickly; its cost is then used as upper bound for the
    exact branch & bound search of the compiled core (imcore), which either improves the set or
    proves that no better one exists (within time_limit seconds CPU time, None = unlimited).
    """
    rng = random.Random(seed)
    n_grid = int(round((stop - start) / step)) + 1
    spacing_idx = int(math.ceil(spacing / step - 1e-9))
    threshold_idx = int(math.floor(threshold / step + 1e-9))
    if (count - 1) * spacing_idx > n_grid - 1:
        raise ValueError('Band too narrow for %d carriers with %.3f MHz spacing' % (count, spacing))

    best_cost, best = None, None
    for _ in range(restarts):
        cost, channels = anneal(n_grid, count, spacing_idx, threshold_idx, iterations, rng)
        if best_cost is None or cost < best_cost:
            best_cost, best = cost, channels
        if best_cost == 0:
            break

    proven = best_cost == 0
    if exact and not proven:
        if imcore is None:
            print('Note: compiled core not built, skipping exact search '
                  '(run: python setup.py build_ext --inplace)')
        else:
            weights = tuple(weight for _, _, weight in IM_TYPES)
            cost, channels, proven, _ = imcore.search(
                n_grid, count, spacing_idx, threshold_idx, weights, best_cost,
                -1.0 if time_limit is None else time_limit)
            if channels is not None:
                best_cost, best = cost, channels

    state = IMState(best, threshold_idx)
    assert state.cost() == best_cost
    hits = {name: state.hits(name) for _, name, _ in IM_TYPES}
    freqs = [round(start + c * step, 6) for c in best]
    return freqs, hits, proven


def main():
    p = argparse.ArgumentParser(description='Find an intermodulation free (as far as possible) frequency set.')
    p.add_argument('--start', type=float, required=True, help='lowest allowed frequency [MHz]')
    p.add_argument('--stop', type=float, required=True, help='highest allowed frequency [MHz]')
    p.add_argument('--spacing', type=float, required=True, help='minimum carrier spacing [MHz]')
    p.add_argument('--count', type=int, required=True, help='number of frequencies')
    p.add_argument('--threshold', type=float, default=0.035,
                   help='min. distance of an IM product to a carrier [MHz] (default 0.035)')
    p.add_argument('--step', type=float, default=0.01, help='frequency grid [MHz] (default 0.01)')
    p.add_argument('--iterations', type=int, default=20000, help='annealing steps per restart')
    p.add_argument('--restarts', type=int, default=3, help='number of annealing restarts')
    p.add_argument('--seed', type=int, default=None, help='random seed for reproducible results')
    p.add_argument('--time-limit', type=float, default=60.0,
                   help='CPU time limit of the exact search [s] (default 60, negative = unlimited)')
    p.add_argument('--no-exact', action='store_true', help='annealing only, skip the exact search')
    args = p.parse_args()

    time_limit = None if args.time_limit < 0 else args.time_limit
    freqs, hits, proven = optimize(args.start, args.stop, args.spacing, args.count, args.threshold,
                                   args.step, args.iterations, args.restarts, args.seed,
                                   time_limit, not args.no_exact)
    print('### Optimized Frequency List: ', freqs)
    print('### Critical IM hits: ', ', '.join('%s: %d' % kv for kv in hits.items()))
    print('### Proven optimal' if proven else '### Best set found, NOT proven optimal (time limit reached)')


if __name__ == '__main__':
    main()

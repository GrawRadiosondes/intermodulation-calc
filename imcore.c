/*
 * Exact branch & bound search for intermodulation free frequency sets (compiled core of optimize.py).
 *
 * Carriers are integer grid indices 0..n_grid-1 and are placed in ascending order (depth first).
 * Placing a carrier can only add IM hits, never remove them, so the partial cost is a lower bound
 * for every completion. A branch is cut as soon as its partial cost reaches the best known cost.
 * If the search space is exhausted, the best set found is proven optimal.
 *
 * Build: python setup.py build_ext --inplace
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define N_ORDERS 3
static const int K[N_ORDERS] = {2, 3, 4};  /* product = k*a - (k-1)*b  ->  IM3, IM5, IM7 */

typedef struct {
    int n, count, spacing, t, off, lo_store, hi_store;
    long w[N_ORDERS];
    int *pos, *best_pos, *occ, *cnt[N_ORDERS];
    long best;
    long long nodes;
    int found, aborted, interrupted;
    clock_t deadline;
    int has_deadline;
} Search;

static inline int stored(const Search *s, int p) { return p >= s->lo_store && p <= s->hi_store; }

static inline int near_carriers(const Search *s, int p) {
    int sum = 0;
    for (int d = -s->t; d <= s->t; d++) sum += s->occ[s->off + p + d];
    return sum;
}

static void update_products(Search *s, int m, int x, int sign) {
    for (int j = 0; j < m; j++) {
        int a = s->pos[j];
        for (int o = 0; o < N_ORDERS; o++) {
            int k = K[o];
            int p1 = k * x - (k - 1) * a, p2 = k * a - (k - 1) * x;
            if (stored(s, p1)) s->cnt[o][s->off + p1] += sign;
            if (stored(s, p2)) s->cnt[o][s->off + p2] += sign;
        }
    }
}

static void dfs(Search *s, int m, long cost) {
    if (s->aborted) return;
    if ((++s->nodes & 0xFFFF) == 0) {
        if (s->has_deadline && clock() > s->deadline) { s->aborted = 1; return; }
        if (PyErr_CheckSignals() != 0) { s->aborted = s->interrupted = 1; return; }
    }
    if (m == s->count) {
        if (cost < s->best) {
            s->best = cost;
            s->found = 1;
            memcpy(s->best_pos, s->pos, sizeof(int) * s->count);
        }
        return;
    }
    int lo = m ? s->pos[m - 1] + s->spacing : 0;
    int hi = s->n - 1 - (s->count - 1 - m) * s->spacing;
    for (int x = lo; x <= hi && s->best > 0 && !s->aborted; x++) {
        long c = cost;
        /* products of already placed carriers landing near x */
        for (int o = 0; o < N_ORDERS; o++) {
            int sum = 0;
            for (int d = -s->t; d <= s->t; d++) sum += s->cnt[o][s->off + x + d];
            c += s->w[o] * sum;
        }
        if (c >= s->best) continue;
        s->occ[s->off + x] = 1;
        /* new products of x with every placed carrier, checked against all carriers incl. x */
        for (int j = 0; j < m && c < s->best; j++) {
            int a = s->pos[j];
            for (int o = 0; o < N_ORDERS; o++) {
                int k = K[o];
                int p1 = k * x - (k - 1) * a, p2 = k * a - (k - 1) * x;
                if (stored(s, p1)) c += s->w[o] * near_carriers(s, p1);
                if (stored(s, p2)) c += s->w[o] * near_carriers(s, p2);
            }
        }
        if (c < s->best) {
            update_products(s, m, x, +1);
            s->pos[m] = x;
            dfs(s, m + 1, c);
            update_products(s, m, x, -1);
        }
        s->occ[s->off + x] = 0;
    }
}

static PyObject *search(PyObject *self, PyObject *args, PyObject *kwargs) {
    static char *kwlist[] = {"n_grid", "count", "spacing", "threshold", "weights",
                             "upper_bound", "time_limit", NULL};
    Search s;
    memset(&s, 0, sizeof s);
    long upper_bound;
    double time_limit = -1;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "iiii(lll)l|d", kwlist, &s.n, &s.count, &s.spacing,
                                     &s.t, &s.w[0], &s.w[1], &s.w[2], &upper_bound, &time_limit))
        return NULL;
    if (s.n < 1 || s.count < 1 || s.count > s.n || s.spacing < 1 || s.t < 0) {
        PyErr_SetString(PyExc_ValueError, "invalid search parameters");
        return NULL;
    }
    /* products are only stored if they can come within t of a carrier: [-t, n-1+t] */
    s.lo_store = -s.t;
    s.hi_store = s.n - 1 + s.t;
    s.off = 2 * s.t + 1;
    size_t size = (size_t)s.n + 4 * (size_t)s.t + 2;
    s.pos = calloc(s.count, sizeof(int));
    s.best_pos = calloc(s.count, sizeof(int));
    s.occ = calloc(size, sizeof(int));
    for (int o = 0; o < N_ORDERS; o++) s.cnt[o] = calloc(size, sizeof(int));
    PyObject *result = NULL;
    if (!s.pos || !s.best_pos || !s.occ || !s.cnt[0] || !s.cnt[1] || !s.cnt[2]) {
        PyErr_NoMemory();
        goto cleanup;
    }
    s.best = upper_bound;
    if (time_limit >= 0) {
        s.has_deadline = 1;
        s.deadline = clock() + (clock_t)(time_limit * CLOCKS_PER_SEC);
    }
    if ((s.count - 1) * (long)s.spacing <= s.n - 1) dfs(&s, 0, 0);
    if (s.interrupted) goto cleanup;  /* KeyboardInterrupt is already set */

    PyObject *positions;
    if (!s.found) {
        Py_INCREF(Py_None);
        positions = Py_None;
    } else {
        positions = PyList_New(s.count);
        if (!positions) goto cleanup;
        for (int i = 0; i < s.count; i++) PyList_SET_ITEM(positions, i, PyLong_FromLong(s.best_pos[i]));
    }
    /* proven: search exhausted (or cost 0 reached, which is a trivial lower bound) */
    int proven = !s.aborted || s.best == 0;
    result = Py_BuildValue("(lNOL)", s.best, positions, proven ? Py_True : Py_False, s.nodes);

cleanup:
    free(s.pos);
    free(s.best_pos);
    free(s.occ);
    for (int o = 0; o < N_ORDERS; o++) free(s.cnt[o]);
    return result;
}

static PyMethodDef methods[] = {
    {"search", (PyCFunction)(void (*)(void))search, METH_VARARGS | METH_KEYWORDS,
     "search(n_grid, count, spacing, threshold, weights, upper_bound, time_limit=-1)\n"
     "-> (best_cost, positions or None, proven, nodes)\n\n"
     "Exact search for a carrier set with cost < upper_bound. positions is None if none was found.\n"
     "proven is True if no better set exists. time_limit < 0 means no limit."},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef module = {PyModuleDef_HEAD_INIT, "imcore",
                                    "Compiled branch & bound core for optimize.py", -1, methods};

PyMODINIT_FUNC PyInit_imcore(void) { return PyModule_Create(&module); }

import numpy as np


def init_lattice(L, seed=None):
    if seed is not None:
        np.random.seed(seed)
    return np.random.choice([-1, 1], size=(L, L))


def neighbours(i, j, L):
    return [
        ((i + 1) % L, j),
        ((i - 1) % L, j),
        (i, (j + 1) % L),
        (i, (j - 1) % L),
    ]


def wolff_single_cluster_update(lattice, beta, J=1.0):
    L = lattice.shape[0]
    i0 = np.random.randint(0, L)
    j0 = np.random.randint(0, L)

    seed_spin = lattice[i0, j0]
    p_add = 1.0 - np.exp(-2.0 * beta * J)

    stack = [(i0, j0)]
    in_cluster = np.zeros((L, L), dtype=bool)
    in_cluster[i0, j0] = True

    lattice[i0, j0] = -seed_spin
    cluster_size = 1

    while stack:
        i, j = stack.pop()
        for ni, nj in neighbours(i, j, L):
            if (not in_cluster[ni, nj]) and lattice[ni, nj] == seed_spin:
                if np.random.rand() < p_add:
                    in_cluster[ni, nj] = True
                    lattice[ni, nj] = -seed_spin
                    stack.append((ni, nj))
                    cluster_size += 1

    return cluster_size


def compute_energy(lattice, J=1.0):
    L = lattice.shape[0]
    E = 0.0
    for i in range(L):
        for j in range(L):
            E -= J * lattice[i, j] * lattice[(i + 1) % L, j]
            E -= J * lattice[i, j] * lattice[i, (j + 1) % L]
    return E


def compute_magnetization(lattice):
    return np.sum(lattice)


def susceptibility(m_vals, beta, N):
    M = np.asarray(m_vals)
    M_abs = np.abs(M)
    return beta * (np.mean(M**2) - np.mean(M_abs)**2) / N


def binder_cumulant(m_vals):
    M = np.asarray(m_vals)
    M_abs = np.abs(M)
    m2 = np.mean(M_abs**2)
    m4 = np.mean(M_abs**4)
    return 1.0 - m4 / (3.0 * m2**2)


def autocorrelation_time(vals, max_lag=None):
    x = np.asarray(vals, dtype=float)
    N = x.size
    x -= np.mean(x)
    var = np.var(x)

    if var == 0.0:
        return 0.0

    if max_lag is None:
        max_lag = N // 10

    acf = np.zeros(max_lag)
    for tau in range(1, max_lag):
        acf[tau] = np.sum(x[:N - tau] * x[tau:]) / ((N - tau) * var)

    return 0.5 + np.sum(acf[1:])


def run_wolff_simulation(L, beta, n_updates, J=1.0, seed=None):
    lattice = init_lattice(L, seed)
    M_list = []

    for _ in range(n_updates):
        wolff_single_cluster_update(lattice, beta, J)
        M_list.append(compute_magnetization(lattice))

    return np.array(M_list)


def run_wolff_with_measurements(L, beta, n_updates, J=1.0, seed=None):
    lattice = init_lattice(L, seed)
    N = L * L

    m_list = []
    e_list = []

    for _ in range(n_updates):
        wolff_single_cluster_update(lattice, beta, J)
        m_list.append(compute_magnetization(lattice))
        e_list.append(compute_energy(lattice, J))

    return {
        "mean_M": np.mean(m_list),
        "mean_E": np.mean(e_list),
        "chi": susceptibility(m_list, beta, N),
        "U4": binder_cumulant(m_list),
        "tau_M": autocorrelation_time(m_list),
        "tau_E": autocorrelation_time(e_list),
        "m_list": np.array(m_list),
        "e_list": np.array(e_list),
    }


def run_linear_annealing(L, beta_start, beta_end, steps, J=1.0, seed=None):
    lattice = init_lattice(L, seed)
    betas = np.linspace(beta_start, beta_end, steps)

    E_list = []
    M_list = []

    for beta in betas:
        wolff_single_cluster_update(lattice, beta, J)
        E_list.append(compute_energy(lattice, J) / (L * L))
        M_list.append(np.abs(compute_magnetization(lattice)) / (L * L))

    return {
        "beta": betas,
        "E": np.array(E_list),
        "M": np.array(M_list),
    }


if __name__ == "__main__":
    L = 32
    T = 2.27
    beta = 1.0 / T
    n_updates = 5000

    results = run_wolff_with_measurements(
        L=L,
        beta=beta,
        n_updates=n_updates,
        J=1.0,
        seed=1234,
    )

    print("Mean magnetization:", results["mean_M"])
    print("Mean energy:", results["mean_E"])
    print("Susceptibility:", results["chi"])
    print("Binder cumulant:", results["U4"])
    print("Autocorrelation time (M):", results["tau_M"])
    print("Autocorrelation time (E):", results["tau_E"])

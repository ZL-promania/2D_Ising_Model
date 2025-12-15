import numpy as np
import matplotlib.pyplot as plt


def main():
    size = 30          # size N x N
    T = 2.0            # k_B = 1
    n_sweeps = 1000    # Metropolis steps

    ising = get_one_sample(size_of_sample=size,
                           temperature=T, n_sweeps=n_sweeps)

    plot_spins(ising)

    energy_s = local_energies(ising)
    plt.figure()
    plt.hist(energy_s.ravel(), bins=50, density=True, alpha=0.7)
    plt.xlabel("local energy")
    plt.ylabel("probability density")
    plt.title("Local energy distribution")
    plt.show()


def get_one_sample(size_of_sample, temperature, n_sweeps=1000, J=1.0):
    """get random initial spin +1/-1"""

    S = np.random.choice([-1, 1], size=(size_of_sample, size_of_sample))
    print("size of the sample:", S.shape)

    energy = total_energy(S, J=J)
    print("initial energy:", energy)

    for sweep in range(n_sweeps):
        S, energy = metropolis_sweep(
            S, temperature, J=J, current_energy=energy)
        if (sweep + 1) % 100 == 0:
            print(f"epoch #{sweep + 1}, current energy: {energy}")

    print(f"epoch #{n_sweeps}, final energy: {energy}")
    return S


def metropolis_sweep(S, T, J=1.0, current_energy=None):
    """randomly flip one spin, accept the flip based on Boltzmann factor"""

    N = S.shape[0]
    if current_energy is None:
        current_energy = total_energy(S, J=J)

    beta = 1.0 / T

    i = np.random.randint(0, N)
    j = np.random.randint(0, N)

    s_ij = S[i, j]

    nn_sum = (
        S[(i + 1) % N, j] +
        S[(i - 1) % N, j] +
        S[i, (j + 1) % N] +
        S[i, (j - 1) % N]
    )
    # sum of neighbor spins

    dE = 2.0 * J * s_ij * nn_sum  # change of energy

    if dE <= 0.0 or np.random.rand() < np.exp(-beta * dE):  # accept
        S[i, j] = -s_ij
        current_energy += dE

    return S, current_energy


def total_energy(S, J=1.0):
    """total energy, periodic B.C."""
    E = 0.0
    E -= J * np.sum(S * np.roll(S, shift=1, axis=0))
    E -= J * np.sum(S * np.roll(S, shift=1, axis=1))
    return E


def local_energies(S, J=1.0):
    """local energy, periodic B.C."""
    nn_sum = (
        np.roll(S, shift=1, axis=0) +
        np.roll(S, shift=-1, axis=0) +
        np.roll(S, shift=1, axis=1) +
        np.roll(S, shift=-1, axis=1)
    )
    E_local = -0.5 * J * S * nn_sum  # divide double counting
    return E_local


def plot_spins(S):
    """plot an Ising state"""
    plt.figure()
    plt.imshow(S, cmap="coolwarm", interpolation="nearest")
    plt.colorbar(label="spin")
    plt.title("Ising spin configuration")
    plt.axis("off")
    plt.show()


if __name__ == "__main__":
    main()

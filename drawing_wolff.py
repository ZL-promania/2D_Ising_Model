import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def magnetization_autocorrelation(M_series):
    M = np.asarray(M_series, dtype=np.float64)
    M = M - M.mean()
    n = len(M)

    f = np.fft.fft(M, n=2*n)
    acf = np.fft.ifft(f * np.conjugate(f)).real
    acf = acf[:n]

    acf /= acf[0]
    return acf

def integrated_autocorrelation_time(M_series, max_lag=None, use_positive_window=True):
    rho = magnetization_autocorrelation(M_series)
    n = len(rho)

    if max_lag is None:
        if use_positive_window:
            positive_idx = np.where(rho > 0)[0]
            if len(positive_idx) == 0:
                max_lag = 0
            else:
                max_lag = int(positive_idx[-1])
        else:
            max_lag = n - 1
    else:
        max_lag = min(max_lag, n - 1)

    tau_int = 0.5 + np.sum(rho[1:max_lag+1])
    return tau_int, rho

def run_wolff_time_series(L, T, n_therm, n_steps, J=1.0, seed=None):
    if seed is not None:
        np.random.seed(seed)

    beta = 1.0 / T
    lattice = init_lattice(L)

    for _ in range(n_therm):
        wolff_single_cluster_update(lattice, beta, J)

    M_list = []
    for _ in range(n_steps):
        wolff_single_cluster_update(lattice, beta, J)
        M_list.append(measure_magnetization(lattice))

    return np.array(M_list)

def compute_tau_vs_T(L, T_list, n_therm, n_steps, J=1.0, seed_base=1234,
                     max_lag=None, use_positive_window=True):
    tau_list = []

    for idx, T in enumerate(T_list):
        seed = seed_base + idx
        M_series = run_wolff_time_series(L, T, n_therm, n_steps, J=J, seed=seed)
        tau_int, rho = integrated_autocorrelation_time(
            M_series, max_lag=max_lag, use_positive_window=use_positive_window
        )
        tau_list.append(tau_int)
        print(f"T = {T:.3f}, tau_int = {tau_int:.3f}")

    return np.array(tau_list)

def plot_tau_vs_T(L, T_list, tau_list):
    plt.figure(figsize=(7,5))
    plt.plot(T_list, tau_list, marker='o')
    plt.xlabel("Temperature T")
    plt.ylabel(r"Integrated autocorrelation time $\tau_{\mathrm{int}}$")
    plt.title(fr"Wolff 2D Ising: $\tau_{{\mathrm{{int}}}}$ vs $T$ (L={L})")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()




def compute_tau_vs_L(T, L_list, n_therm_per_spin, n_steps_per_spin,
                     J=1.0, seed_base=2025,
                     max_lag=None, use_positive_window=True):
    tau_list = []

    for idx, L in enumerate(L_list):
        n_therm = int(n_therm_per_spin * L * L)
        n_steps = int(n_steps_per_spin * L * L)

        seed = seed_base + idx
        M_series = run_wolff_time_series(L, T, n_therm, n_steps, J=J, seed=seed)

        
        tau_int, rho = integrated_autocorrelation_time(
            M_series, max_lag=max_lag, use_positive_window=use_positive_window
        )
        tau_list.append(tau_int)
        print(f"L = {L:3d}, n_therm = {n_therm}, n_steps = {n_steps}, tau_int = {tau_int:.3f}")

    return np.array(tau_list)


def plot_tau_vs_L(T, L_list, tau_list, loglog=True):
    plt.figure(figsize=(7,5))
    if loglog:
        plt.loglog(L_list, tau_list, marker='o')
    else:
        plt.plot(L_list, tau_list, marker='o')
    plt.xlabel("System size L")
    plt.ylabel(r"Integrated autocorrelation time $\tau_{\mathrm{int}}$")
    plt.title(fr"Wolff 2D Ising: $\tau_{{\mathrm{{int}}}}$ vs $L$ (T={T:.3f})")
    plt.grid(alpha=0.3, which="both")
    plt.tight_layout()
    plt.show()




def run_ising_M_vs_T(L_list, T_list, n_therm_per_spin, n_steps_per_spin, J=1.0):
    M_data = {}
    M_std = {}

    for L in L_list:
        M_data[L] = {}
        M_std[L] = {}

        n_therm = int(n_therm_per_spin * L * L)
        n_steps = int(n_steps_per_spin * L * L)

        for idx, T in enumerate(T_list):
            beta = 1.0 / T
            seed = 1000 + 10*L + idx

            M_series = run_wolff_time_series(L, T, n_therm, n_steps, J=J, seed=seed)
            absM = np.abs(M_series) / (L * L)

            M_data[L][T] = absM.mean()
            M_std[L][T] = absM.std() / np.sqrt(len(absM))

            print(f"L={L:3d}, T={T:.3f},  <|M|>={M_data[L][T]:.4f}")

    return M_data, M_std


def plot_raw_M_vs_T(L_list, T_list, M_data, M_std):
    plt.figure(figsize=(7,5))

    for L in L_list:
        M_list = np.array([M_data[L][T] for T in T_list])
        err_list = np.array([M_std[L][T] for T in T_list])
        plt.errorbar(T_list, M_list, yerr=err_list, marker='o', label=f"L={L}")

    plt.xlabel("Temperature T")
    plt.ylabel(r"$\langle |M| \rangle$")
    plt.title("Raw magnetization curves (no scaling)")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_scaled_M(L_list, T_list, M_data, beta=1/8, nu=1, Tc=2.269185):
    plt.figure(figsize=(7,5))

    for L in L_list:
        x = (np.array(T_list) - Tc) * (L ** (1/nu))
        y = np.array([M_data[L][T] for T in T_list]) * (L ** (beta/nu))

        plt.plot(x, y, marker='o', label=f"L={L}")

    plt.xlabel(r"$(T - T_c)L^{1/\nu}$")
    plt.ylabel(r"$\langle |M| \rangle L^{\beta/\nu}$")
    plt.title("Finite-size scaling collapse of magnetization")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()




def make_ising_gif(filename="ising_wolff.gif", 
                   L=64, T=2.27, frames=100, interval=150):
    beta = 1.0 / T
    lattice = init_lattice(L)

    for _ in range(200):
        wolff_single_cluster_update(lattice, beta)

    fig, ax = plt.subplots(figsize=(6,6))
    img = ax.imshow(lattice, cmap="gray", vmin=-1, vmax=1, animated=True)
    ax.axis("off")

    def update(frame):
        wolff_single_cluster_update(lattice, beta)
        img.set_data(lattice)
        return [img]

    ani = animation.FuncAnimation(
        fig, update, frames=frames, interval=interval, blit=True
    )

    ani.save(filename, writer="pillow")
    plt.close(fig)
    print(f"GIF saved as {filename}")

make_ising_gif("ising_demo.gif", L=64, T=2.27, frames=120, interval=120)




if __name__ == "__main__":
    L = 32
    T_crit = 2.27
    J = 1.0

    #swithches to control which analysis to run
    RUN_TAU_VS_T = True
    RUN_TAU_VS_L = True
    RUN_M_VS_T   = False
    RUN_GIF      = False

    if RUN_TAU_VS_T:
        T_list = np.linspace(1.5, 3.5, 9)
        n_therm = 5000
        n_steps = 50000

        tau_T = compute_tau_vs_T(
            L, T_list, n_therm, n_steps,
            max_lag=None, use_positive_window=True
        )
        plot_tau_vs_T(L, T_list, tau_T)

    if RUN_TAU_VS_L:
        L_list = np.array([24, 32, 48, 64])

        n_therm_per_spin = 0.5
        n_steps_per_spin = 2.0

        tau_L = compute_tau_vs_L(
            T_crit, L_list,
            n_therm_per_spin,
            n_steps_per_spin
        )
        plot_tau_vs_L(T_crit, L_list, tau_L, loglog=True)


    if RUN_M_VS_T:
        T_list = np.concatenate([
            np.linspace(2.0, 2.6, 13),
            np.arange(2.7, 7.0 + 1e-9, 0.5)
        ])

        L_list = [16, 24, 32]

        n_therm_per_spin = 0.5
        n_steps_per_spin = 2.0

        M_data, M_std = run_ising_M_vs_T(
            L_list, T_list,
            n_therm_per_spin,
            n_steps_per_spin
        )

        plot_raw_M_vs_T(L_list, T_list, M_data, M_std)
        plot_scaled_M(L_list, T_list, M_data)


    if RUN_GIF:
        make_ising_gif(
            "ising_demo.gif",
            L=64,
            T=2.27,
            frames=120,
            interval=120
        )




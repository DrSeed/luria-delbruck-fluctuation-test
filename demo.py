# Self-contained demo: random mutation vs induced resistance, and three estimators of the rate.
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from fluctuation import grow_culture, induced_culture, m_p0, m_lea_coulson, m_mss_mle

N0, GENS, MU, CULTURES = 1000, 17, 1e-8, 60


def main():
    os.makedirs("figures", exist_ok=True)
    os.makedirs("results", exist_ok=True)
    rng = np.random.default_rng(1943)
    nt = N0 * 2 ** GENS
    true_m = MU * N0 * (2 ** GENS - 1)
    random_counts = np.array([grow_culture(N0, GENS, MU, rng) for _ in range(CULTURES)])
    induced_counts = np.array([induced_culture(nt, random_counts.mean() / nt, rng) for _ in range(CULTURES)])

    rows = []
    for name, counts in [("random mutation (Luria-Delbruck)", random_counts), ("induced on the plate (Poisson)", induced_counts)]:
        rows.append({"hypothesis": name, "mean": round(counts.mean(), 2), "variance": round(counts.var(ddof=1), 1),
                     "variance_to_mean": round(counts.var(ddof=1) / counts.mean(), 1), "max": int(counts.max()),
                     "cultures_with_zero": int((counts == 0).sum())})
    table = pd.DataFrame(rows)
    est = pd.DataFrame([{"method": n, "m_estimate": round(f(random_counts), 3),
                         "mutation_rate_per_division": f"{f(random_counts) / nt:.2e}"}
                        for n, f in [("p0 method", m_p0), ("Lea-Coulson median", m_lea_coulson), ("MSS maximum likelihood", m_mss_mle)]])
    est.loc[len(est)] = {"method": "truth", "m_estimate": round(true_m, 3), "mutation_rate_per_division": f"{MU:.2e}"}
    table.to_csv("results/summary.csv", index=False)
    est.to_csv("results/rate_estimates.csv", index=False)

    reps = []
    for rep in range(40):
        c = np.array([grow_culture(N0, GENS, MU, rng) for _ in range(CULTURES)])
        reps.append({"p0": m_p0(c), "Lea-Coulson": m_lea_coulson(c), "MSS-MLE": m_mss_mle(c)})
    reps = pd.DataFrame(reps)

    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))
    bins = np.concatenate([[-0.5, 0.5], np.logspace(0, np.log10(max(random_counts.max(), 2) + 1), 15)])
    axes[0].hist(induced_counts, bins=bins, alpha=0.6, label="induced (Poisson)", color="#1f77b4")
    axes[0].hist(random_counts, bins=bins, alpha=0.6, label="random mutation", color="#c0504d")
    axes[0].set_xscale("symlog")
    axes[0].set_xlabel("resistant colonies per culture")
    axes[0].set_ylabel("cultures")
    axes[0].set_title("Same mean, wildly different spread")
    axes[0].legend()
    axes[1].bar(table["hypothesis"].str.split(" ").str[0], table["variance_to_mean"], color=["#c0504d", "#1f77b4"])
    axes[1].axhline(1, color="k", ls="--", lw=1, label="Poisson (variance = mean)")
    axes[1].set_yscale("log")
    axes[1].set_ylabel("variance / mean")
    axes[1].set_title("The 1943 test: jackpots mean random mutation")
    axes[1].legend()
    axes[2].boxplot([reps[c].dropna() for c in reps.columns], showfliers=False)
    axes[2].set_xticks([1, 2, 3], reps.columns)
    axes[2].axhline(true_m, color="k", ls="--", lw=1, label=f"true m = {true_m:.2f}")
    axes[2].set_ylabel("estimated mutations per culture (m)")
    axes[2].set_title("40 repeat experiments of 60 cultures")
    axes[2].legend()
    fig.tight_layout()
    fig.savefig("figures/demo.png", dpi=120)
    print(table.to_string(index=False))
    print(est.to_string(index=False))
    print(reps.describe().round(3))


if __name__ == "__main__":
    main()

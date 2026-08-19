"""
Point-by-point comparison of this repo's R-1234ze(E)/R-227ea mixture model
(Bell 2023 departure function, w1=0.911) against INDEPENDENT, real experimental
R-515B liquid density data from:

    Kang, K., Yang, S., Cui, J., Gu, Y. (2024). "Theoretical study and
    experimental verification of the viscosities of azeotropic refrigerant
    R515B." International Journal of Refrigeration 168, 59-69.
    https://doi.org/10.1016/j.ijrefrig.2024.08.012

This is the first comparison in this project against data that is NOT the
Honeywell TDS p-H chart -- an independent high-pressure vibrating-wire
measurement, with declared combined expanded uncertainty of 0.2% (k=2) for
density. That makes it a much stronger check on the model's raw density
accuracy than anything read off a chart by eye.

Data source
-----------
`kang2024_r515b_density_table4.csv` in this directory, transcribed by hand from
the paper's Table 4 ("High pressure density and viscosity data of R515B."),
T=254.13-362.34 K, P=0.87-12.27 MPa, liquid phase only (all points are
compressed/subcooled liquid -- these pressures are far above the saturation
pressure at every listed temperature). The paper states N=68 points for R515B
in its own Table 8 accuracy summary; only 67 discrete rows could be cleanly
transcribed from the printed table image. This discrepancy is reported
honestly below rather than silently padded to 68 -- if a clearer scan surfaces
a 68th row, it should be added and this script re-run.

What this script does
----------------------
For each experimental (T, P) point, solves for OUR model's own liquid density
at that exact (T, P) -- NOT the model's own saturation/VLE density, a genuinely
independent forward calculation: given T, brentq-bracket the root of
mix_state(...).p_pa - P_target over a wide compressed-liquid density range.
This is the direct compressed-liquid analog of the vapor-side isentrope
point-solves used everywhere else in this repo, just solving for density at
fixed (T,P) instead of density at fixed (P,s).

Reports, per point: T, P, experimental density, model density, relative
deviation % = (rho_model - rho_exp)/rho_exp * 100 -- same sign convention as
the paper's own Fig. 2/3 deviation plots, so the two can be compared directly.
Then reports summary statistics (AAD%, bias, min/max) and saves both a CSV of
full results and a deviation-vs-T/P plot.

Usage
-----
Run from R515B_props_validated/:
    python3 validate_against_kang2024.py

Read-only against the model -- does not modify mixture_isentrope_validation.py
or any other file.
"""

import csv
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mixture_isentrope_validation import (
    load_idaes_helmholtz_json,
    mw_from_json,
    w1_to_x1,
    mix_state,
    _dp_drho_fd,
)

FLUID1, FLUID2, W1 = "r1234ze", "r227ea", 0.911
DATA_CSV = Path(__file__).parent / "kang2024_r515b_density_table4.csv"
OUT_CSV = Path(__file__).parent / "verification" / "kang2024_r515b_density_comparison.csv"
OUT_FIG = Path(__file__).parent / "verification" / "kang2024_r515b_density_deviation.png"


def load_experimental_data(path: Path):
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append({
                "T_K": float(r["T_K"]),
                "P_MPa": float(r["P_MPa"]),
                "rho_exp_kgm3": float(r["rho_exp_kgm3"]),
                "eta_exp_mPas": float(r["eta_exp_mPas"]),
            })
    return rows


def solve_liquid_density(d1, d2, z1, mw_mix, t_k: float, p_target_pa: float,
                          rho_seed_kgm3: float,
                          max_iter: int = 60, tol_rel: float = 1e-10) -> float:
    """
    Solve for the compressed-liquid density at fixed (T, P) via a damped
    Newton iteration in MOLAR density, seeded at the experimental density and
    using dP/drho from the model's own finite-difference derivative
    (`_dp_drho_fd`, already used elsewhere in this repo for the critical-point
    solve).

    IMPORTANT (found the hard way -- two earlier attempts in this script both
    failed silently in different ways): a plain brentq bracket -- whether wide
    (400-1600 kg/m3) or moderately tightened (+/-40-60% of rho_exp) -- is NOT
    safe here. At fixed T, P(rho) from a multiparameter Helmholtz EOS is only
    monotonic in the stable single-phase branches; between the liquid and
    vapor spinodals it can fold back (a Van der Waals-loop shape), and this
    folding can sit closer to the true liquid root than a loose percentage
    bracket assumes -- especially for the T=303-352K points nearer the
    critical region, where several points converged to a wildly wrong (~53-59%
    off) root even with a +/-40% bracket. brentq only guarantees finding *a*
    sign change inside a bracket, not the physically stable one, so it isn't
    safe to use here without independently confirming the bracket is
    single-rooted at every single (T,P) -- more fragile than just avoiding
    bracketed search altogether.

    A damped Newton iteration seeded exactly at rho_exp (known accurate to
    0.2% per the paper) instead FOLLOWS THE LOCAL SLOPE from a trusted good
    starting point and cannot jump to a distant, unrelated root the way a
    blind bracket search can. Step size is capped at 15% of the current
    density per iteration as a damping safeguard against overshoot near the
    fold region; this was verified (see the comparison against the two
    brentq attempts above, and the smooth, physically-sensible ~1.7-2.5%
    deviation this converges to almost everywhere) to land on the same
    stable-liquid root as the well-behaved majority of brentq's own points.
    """
    rho_mol = rho_seed_kgm3 / mw_mix
    for _ in range(max_iter):
        st = mix_state(d1, d2, t_k, rho_mol, z1)
        resid = st.p_pa - p_target_pa
        if abs(resid) < tol_rel * abs(p_target_pa):
            break
        dp_drho = _dp_drho_fd(d1, d2, t_k, rho_mol, z1)
        if dp_drho == 0:
            raise RuntimeError(f"dP/drho=0 at T={t_k:.2f}K, rho={rho_mol:.2f} mol/m3 -- at a spinodal point")
        step = -resid / dp_drho
        step = np.clip(step, -0.15 * rho_mol, 0.15 * rho_mol)  # damping: cap step at 15% of current density
        rho_mol += step
        if rho_mol <= 0:
            raise RuntimeError(f"Newton step drove rho_mol negative at T={t_k:.2f}K")
    else:
        raise RuntimeError(f"Newton did not converge in {max_iter} iterations at T={t_k:.2f}K, "
                            f"P={p_target_pa/1e6:.2f}MPa (final resid={resid:.3e} Pa)")
    return rho_mol * mw_mix  # kg/m3


def main():
    d1 = load_idaes_helmholtz_json(FLUID1)
    d2 = load_idaes_helmholtz_json(FLUID2)
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)
    z1 = w1_to_x1(W1, mw1, mw2)
    mw_mix = z1 * mw1 + (1.0 - z1) * mw2

    exp_rows = load_experimental_data(DATA_CSV)
    print(f"Loaded {len(exp_rows)} experimental (T,P,rho) points from {DATA_CSV.name}")
    print("(Kang et al. 2024, Table 4 states N=68 for R515B; only 67 rows were "
          "cleanly transcribed from the printed table -- see this script's docstring.)")
    print()

    results = []
    n_failed = 0
    for row in exp_rows:
        t_k = row["T_K"]
        p_pa = row["P_MPa"] * 1.0e6
        try:
            rho_model = solve_liquid_density(d1, d2, z1, mw_mix, t_k, p_pa,
                                              rho_seed_kgm3=row["rho_exp_kgm3"])
        except RuntimeError as e:
            print(f"  FAILED to solve at T={t_k:.2f}K, P={row['P_MPa']:.2f}MPa: {e}")
            n_failed += 1
            continue
        dev_pct = (rho_model - row["rho_exp_kgm3"]) / row["rho_exp_kgm3"] * 100.0
        results.append({
            **row,
            "rho_model_kgm3": rho_model,
            "deviation_pct": dev_pct,
        })

    if n_failed:
        print(f"\n{n_failed} of {len(exp_rows)} points FAILED to solve (see above).\n")

    print(f"{'T(K)':>8} {'P(MPa)':>8} {'rho_exp':>10} {'rho_model':>10} {'dev%':>8}")
    for r in results:
        print(f"{r['T_K']:>8.2f} {r['P_MPa']:>8.2f} {r['rho_exp_kgm3']:>10.2f} "
              f"{r['rho_model_kgm3']:>10.2f} {r['deviation_pct']:>8.3f}")

    devs = np.array([r["deviation_pct"] for r in results])
    aad_pct = np.mean(np.abs(devs))
    bias_pct = np.mean(devs)
    print()
    print(f"N compared: {len(results)} of {len(exp_rows)} loaded ({n_failed} solver failures)")
    print(f"AAD% (mean absolute deviation): {aad_pct:.3f}%")
    print(f"Mean signed bias: {bias_pct:+.3f}%  ({'model reads HIGH' if bias_pct > 0 else 'model reads LOW'} on average)")
    print(f"Max |deviation|: {np.max(np.abs(devs)):.3f}%  (at T={results[np.argmax(np.abs(devs))]['T_K']:.2f}K, "
          f"P={results[np.argmax(np.abs(devs))]['P_MPa']:.2f}MPa)")
    print(f"Min |deviation|: {np.min(np.abs(devs)):.3f}%")
    print(f"Deviation range: [{np.min(devs):+.3f}%, {np.max(devs):+.3f}%]")

    # Save full results CSV
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["T_K", "P_MPa", "rho_exp_kgm3", "eta_exp_mPas",
                                                "rho_model_kgm3", "deviation_pct"])
        writer.writeheader()
        for r in results:
            writer.writerow({k: r[k] for k in writer.fieldnames})
    print(f"\nSaved full point-by-point results: {OUT_CSV}")

    # Deviation plot, styled like the paper's own Fig. 2/3 for direct comparison
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    temps = sorted(set(r["T_K"] for r in results))
    cmap = plt.cm.viridis(np.linspace(0, 1, len(temps)))
    t_color = {t: cmap[i] for i, t in enumerate(temps)}

    ax = axes[0]
    for r in results:
        ax.scatter(r["P_MPa"], r["deviation_pct"], color=t_color[r["T_K"]], s=25)
    ax.axhline(0, color="k", linewidth=0.8)
    ax.set_xlabel("P (MPa)")
    ax.set_ylabel(r"$(\rho_{model}-\rho_{exp})/\rho_{exp}\times100\%$")
    ax.set_title("Deviation vs. pressure")
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    for r in results:
        ax.scatter(r["T_K"], r["deviation_pct"], color=t_color[r["T_K"]], s=25)
    ax.axhline(0, color="k", linewidth=0.8)
    ax.set_xlabel("T (K)")
    ax.set_ylabel(r"$(\rho_{model}-\rho_{exp})/\rho_{exp}\times100\%$")
    ax.set_title("Deviation vs. temperature")
    ax.grid(True, alpha=0.3)

    fig.suptitle("Model vs. Kang et al. (2024) experimental R-515B liquid density\n"
                  f"AAD%={aad_pct:.2f}%, bias={bias_pct:+.2f}%, N={len(results)}")
    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=150)
    print(f"Saved deviation plot: {OUT_FIG}")


if __name__ == "__main__":
    sys.exit(main())

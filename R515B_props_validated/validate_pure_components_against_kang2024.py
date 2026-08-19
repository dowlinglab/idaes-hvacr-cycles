"""
Extension of validate_against_kang2024.py to the PURE-COMPONENT tables in the
same paper (Kang et al. 2024): Table 2 (R227ea) and Table 3 (R1234ze(E)),
high-pressure liquid density, same T/P grid and same vibrating-wire method as
the R515B mixture data already checked in validate_against_kang2024.py.

Purpose
-------
The mixture comparison (validate_against_kang2024.py) found a consistent
~1.9% liquid-density underprediction. This project's standing open question
(flagged since 2026-08-12, "Pure R-227ea's real saturated-liquid density near
25C still unconfirmed") is WHERE that bias comes from: the pure-component
Helmholtz EOS layer (r1234ze.json / r227ea.json, used directly, no mixing
rule), or the Bell (2023) mixing/departure-function layer on top of it. The
project's earlier pure-fluid checks (2026-08-12) were a single external
data point each (one for R-1234ze(E), matched to ~0.003%; none at all for
R-227ea, since no source was found at the time). This script uses Kang et
al.'s own real measurements of BOTH pure fluids -- same lab, same method, same
T/P grid as the mixture data -- to settle this properly across a wide range
instead of one spot value.

Method
------
Identical density solve to validate_against_kang2024.py (damped Newton in
molar density, seeded at the experimental value, using the model's own
_dp_drho_fd) -- just with z1 fixed to the pure-component limit (z1=1.0 for
R-1234ze(E), z1=0.0 for R227ea) instead of the R515B mixture composition.
Confirmed beforehand that mix_state evaluates cleanly at both pure limits
(no divide-by-zero in the corresponding-states mixing formulas).

Usage
-----
Run from R515B_props_validated/:
    python3 validate_pure_components_against_kang2024.py

Read-only against the model.
"""

import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mixture_isentrope_validation import (
    load_idaes_helmholtz_json,
    mw_from_json,
    mix_state,
    _dp_drho_fd,
)

FLUID1, FLUID2 = "r1234ze", "r227ea"
HERE = Path(__file__).parent

CASES = [
    {"label": "R1234ze(E) (pure)", "z1": 1.0, "csv": HERE / "kang2024_r1234ze_density_table3.csv",
     "out_csv": HERE / "verification" / "kang2024_r1234ze_density_comparison.csv"},
    {"label": "R227ea (pure)", "z1": 0.0, "csv": HERE / "kang2024_r227ea_density_table2.csv",
     "out_csv": HERE / "verification" / "kang2024_r227ea_density_comparison.csv"},
]
OUT_FIG = HERE / "verification" / "kang2024_pure_components_density_deviation.png"


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


def solve_liquid_density(d1, d2, z1, mw, t_k, p_target_pa, rho_seed_kgm3,
                          max_iter=60, tol_rel=1e-10):
    """Same damped-Newton method as validate_against_kang2024.py (see that
    script's docstring for why a bracketed brentq search is unsafe here)."""
    rho_mol = rho_seed_kgm3 / mw
    for _ in range(max_iter):
        st = mix_state(d1, d2, t_k, rho_mol, z1)
        resid = st.p_pa - p_target_pa
        if abs(resid) < tol_rel * abs(p_target_pa):
            break
        dp_drho = _dp_drho_fd(d1, d2, t_k, rho_mol, z1)
        if dp_drho == 0:
            raise RuntimeError(f"dP/drho=0 at T={t_k:.2f}K")
        step = -resid / dp_drho
        step = np.clip(step, -0.15 * rho_mol, 0.15 * rho_mol)
        rho_mol += step
        if rho_mol <= 0:
            raise RuntimeError(f"Newton step drove rho_mol negative at T={t_k:.2f}K")
    else:
        raise RuntimeError(f"Newton did not converge at T={t_k:.2f}K, P={p_target_pa/1e6:.2f}MPa")
    return rho_mol * mw


def main():
    d1 = load_idaes_helmholtz_json(FLUID1)
    d2 = load_idaes_helmholtz_json(FLUID2)
    mw1 = mw_from_json(d1)
    mw2 = mw_from_json(d2)

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    summary = []

    for case, mw in zip(CASES, [mw1, mw2]):
        exp_rows = load_experimental_data(case["csv"])
        results = []
        n_failed = 0
        for row in exp_rows:
            t_k = row["T_K"]
            p_pa = row["P_MPa"] * 1.0e6
            try:
                rho_model = solve_liquid_density(d1, d2, case["z1"], mw, t_k, p_pa,
                                                   rho_seed_kgm3=row["rho_exp_kgm3"])
            except RuntimeError as e:
                print(f"  [{case['label']}] FAILED at T={t_k:.2f}K, P={row['P_MPa']:.2f}MPa: {e}")
                n_failed += 1
                continue
            dev_pct = (rho_model - row["rho_exp_kgm3"]) / row["rho_exp_kgm3"] * 100.0
            results.append({**row, "rho_model_kgm3": rho_model, "deviation_pct": dev_pct})

        devs = np.array([r["deviation_pct"] for r in results])
        aad = np.mean(np.abs(devs))
        bias = np.mean(devs)
        print(f"=== {case['label']}: N={len(results)}/{len(exp_rows)} (failed={n_failed}) "
              f"AAD%={aad:.3f} bias={bias:+.3f}% range=[{devs.min():+.3f}%,{devs.max():+.3f}%] ===")
        summary.append({"label": case["label"], "N": len(results), "AAD%": aad, "bias%": bias,
                         "min%": float(devs.min()), "max%": float(devs.max())})

        case["out_csv"].parent.mkdir(parents=True, exist_ok=True)
        with open(case["out_csv"], "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["T_K", "P_MPa", "rho_exp_kgm3", "eta_exp_mPas",
                                                    "rho_model_kgm3", "deviation_pct"])
            writer.writeheader()
            for r in results:
                writer.writerow({k: r[k] for k in writer.fieldnames})
        print(f"  Saved: {case['out_csv']}")

        axes[0].scatter([r["P_MPa"] for r in results], devs, label=case["label"], s=20, alpha=0.7)
        axes[1].scatter([r["T_K"] for r in results], devs, label=case["label"], s=20, alpha=0.7)

    for ax, xlabel, title in zip(axes, ["P (MPa)", "T (K)"], ["Deviation vs. pressure", "Deviation vs. temperature"]):
        ax.axhline(0, color="k", linewidth=0.8)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r"$(\rho_{model}-\rho_{exp})/\rho_{exp}\times100\%$")
        ax.set_title(title)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    fig.suptitle("Pure-component model vs. Kang et al. (2024) experimental liquid density")
    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=150)
    print(f"\nSaved deviation plot: {OUT_FIG}")

    print("\n=== Summary ===")
    for s in summary:
        print(s)


if __name__ == "__main__":
    sys.exit(main())

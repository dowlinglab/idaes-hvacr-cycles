"""
Author: Shilpa Narasimhan
Support: Claude AI
Date Created: 10/16/2026

What this does
    Computes TEWI (Total Equivalent Warming Impact) for the cold-storage
    freezer, using Zhuorui's annual electricity results, and makes a
    2 x 2 figure (rows: no limit / 21 C limit, columns: single stage / cascade).

Formula (per system, over the lifetime n)
    TEWI = leaks + end_of_life + electricity
    leaks       = GWP * m * L * n
    end_of_life = GWP * m * (1 - alpha)
    electricity = E * beta * n

Inputs
    from_zhuorui/ColdStorage_annual_results.xlsx  (sheets Annual_No_limit, Annual_21C_limit)

Outputs (folder tewi_outputs/)
    tewi_results_<tag>.csv
    tewi_figure_<tag>.png and .pdf

Run
    cd ~/idaes-hvacr-cycles/Cascade
    conda activate myidaesenv
    python compute_tewi.py
"""

import os

import matplotlib.pyplot as plt
import pandas as pd


# =====================================================================
# 1. SETTINGS  (change things here, nowhere else)
# =====================================================================

# Lifetime, leak rate, recovery: IPCC 2019 Refinement, Vol. 3, Ch. 7,
# Table 7.9, "Medium and large commercial refrigeration".
LIFETIME_YR = 15          # n, high end of 7-15 yr
LEAK_RATE = 0.10          # L, per year, low end of 10-35 %/yr (sensitivity: 0.35)
RECOVERY = 0.70           # alpha, high end of 0-70 %

# Charge: 2 kg/kW for HFC/HFO, 1 kg/kW for CO2 (Bellos and Tzivanidis 2019, Appl. Sci. 9(4), 790)
CAPACITY_KW = 60.0                 # Zhuorui's design capacity at 35 C
SPECIFIC_CHARGE_HFC = 2.0          # kg per kW
SPECIFIC_CHARGE_CO2 = 1.0          # kg per kW
CASCADE_HOT_CHARGE_FRACTION = 1.0  # 1.0 = scenario A, 0.25 = B, 0.10 = C

# GWP, 100 yr, IPCC AR6.  R515B = 0.911 * 1.37 (R1234ze(E)) + 0.089 * 3600 (R227ea)
GWP = {
    "R134a": 1530.0,
    "R1234yf": 0.501,
    "R515B": 0.911 * 1.37 + 0.089 * 3600.0,
    "CO2": 1.0,
}

# Grid emission factor, kg CO2e per kWh, eGRID2023 subregion CO2e output rates
# (lb/MWh * 0.45359237 / 1000)
BETA = {
    "Chicago, IL": 916.1 * 0.45359237 / 1000,     # RFCW
    "Atlanta, GA": 846.0 * 0.45359237 / 1000,     # SRSO
    "San Diego, CA": 430.0 * 0.45359237 / 1000,   # CAMX
    "Tucson, AZ": 706.2 * 0.45359237 / 1000,      # AZNM
}

# Zhuorui's file and sheet names
DATA_FILE = os.path.join("from_zhuorui", "ColdStorage_annual_results.xlsx")
SHEETS = {
    "No limit": "Annual_No_limit",
    "21 C limit": "Annual_21C_limit",
}

# Figure
USE_TEX = False   # True only if LaTeX + cm-super are installed
OUT_DIR = "tewi_outputs"

CITIES = ["Chicago, IL", "Atlanta, GA", "San Diego, CA", "Tucson, AZ"]
FLUIDS = ["R134a", "R1234yf", "R515B"]
CYCLES = ["Single stage", "Cascade (with CO2)"]

# names shown on the figure
FLUID_LABEL = {"R134a": "R-134a", "R515B": "R-515B", "R1234yf": "R-1234yf"}
CO2_LABEL = "R-744"      # CO2 in the cascade


# =====================================================================
# 2. LOAD ZHUORUI'S ELECTRICITY
# =====================================================================

def load_electricity():
    """One table with columns: Case, City, Cycle, Refrigerant, E_kWh_per_yr."""
    tables = []
    for case, sheet in SHEETS.items():
        df = pd.read_excel(DATA_FILE, sheet_name=sheet)
        df = df[["City", "Cycle", "Refrigerant", "Annual chiller electricity (kWh/yr)"]]
        df = df.rename(columns={"Annual chiller electricity (kWh/yr)": "E_kWh_per_yr"})
        df.insert(0, "Case", case)
        tables.append(df)
    return pd.concat(tables, ignore_index=True)


# =====================================================================
# 3. TEWI FOR ONE ROW
# =====================================================================

def refrigerant_charges(cycle, fluid):
    """Returns a list of (refrigerant, charge in kg) for this system."""
    hfc_charge = SPECIFIC_CHARGE_HFC * CAPACITY_KW
    if cycle == "Single stage":
        return [(fluid, hfc_charge)]
    co2_charge = SPECIFIC_CHARGE_CO2 * CAPACITY_KW
    return [(fluid, hfc_charge * CASCADE_HOT_CHARGE_FRACTION), ("CO2", co2_charge)]


def tewi_row(row):
    """Leaks, end of life, electricity and total, all in tonnes CO2e."""
    leaks_kg = 0.0
    end_of_life_kg = 0.0
    for refrigerant, m in refrigerant_charges(row["Cycle"], row["Refrigerant"]):
        leaks_kg += GWP[refrigerant] * m * LEAK_RATE * LIFETIME_YR
        end_of_life_kg += GWP[refrigerant] * m * (1.0 - RECOVERY)

    electricity_kg = row["E_kWh_per_yr"] * BETA[row["City"]] * LIFETIME_YR

    return pd.Series({
        "leaks_t": leaks_kg / 1000,
        "end_of_life_t": end_of_life_kg / 1000,
        "electricity_t": electricity_kg / 1000,
        "TEWI_t": (leaks_kg + end_of_life_kg + electricity_kg) / 1000,
    })


# =====================================================================
# 4. FIGURE
# =====================================================================

def make_figure(results, tag):
    plt.rcParams.update({
        "text.usetex": USE_TEX,
        "mathtext.fontset": "cm",
        "font.family": "serif",
        "font.serif": ["cmr10"],
        "axes.formatter.use_mathtext": True,
        "axes.unicode_minus": False,
        "font.size": 10,
        "ytick.direction": "in",
        "ytick.right": True,
        "hatch.color": "white",
        "hatch.linewidth": 0.8,
    })

    dark, mid, light = "#3a3a3a", "#8c8c8c", "#d0d0d0"
    bar_width = 0.26

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.0), sharey=True)

    panel_list = [
        ("No limit", "Single stage", "Single stage, no condensing limit"),
        ("No limit", "Cascade (with CO2)", "Cascade with " + CO2_LABEL + r", no condensing limit"),
        ("21 C limit", "Single stage", r"Single stage, $21^{\circ}\mathrm{C}$ condensing limit"),
        ("21 C limit", "Cascade (with CO2)", "Cascade with " + CO2_LABEL + r", $21^{\circ}\mathrm{C}$ condensing limit"),
    ]

    for ax, (case, cycle, title) in zip(axes.flat, panel_list):
        for i, city in enumerate(CITIES):
            for j, fluid in enumerate(FLUIDS):
                x = i + (j - 1) * (bar_width + 0.02)
                pick = ((results["Case"] == case) & (results["Cycle"] == cycle)
                        & (results["City"] == city) & (results["Refrigerant"] == fluid))
                r = results[pick].iloc[0]

                # electricity at the bottom, then leaks, then end of life
                ax.bar(x, r["electricity_t"], bar_width, color=dark, edgecolor="white", linewidth=0.6)
                ax.bar(x, r["leaks_t"], bar_width, bottom=r["electricity_t"],
                       color=mid, edgecolor="white", linewidth=0.6, hatch="////")
                ax.bar(x, r["end_of_life_t"], bar_width, bottom=r["electricity_t"] + r["leaks_t"],
                       color=light, edgecolor="white", linewidth=0.6)

                # refrigerant name under each bar
                ax.text(x, -20, FLUID_LABEL[fluid], rotation=90, ha="center", va="top", fontsize=6.5, color="#52514e")

        ax.set_title(title, fontsize=9.5)
        ax.set_xticks(range(len(CITIES)))
        ax.set_xticklabels([c.split(",")[0] for c in CITIES])
        ax.tick_params(axis="x", pad=30, length=0)
        ax.tick_params(axis="y", labelleft=True)
        ax.set_ylabel(r"TEWI over %d yr [t CO$_2$e]" % LIFETIME_YR)
        ax.minorticks_on()
        ax.tick_params(axis="x", which="minor", bottom=False, top=False)
        ax.grid(axis="y", which="major", color="#b0b0b0", linestyle="--", linewidth=0.6)
        ax.grid(axis="y", which="minor", color="#dcdcdc", linestyle=":", linewidth=0.5)
        ax.set_axisbelow(True)

    handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=dark, edgecolor="white"),
        plt.Rectangle((0, 0), 1, 1, facecolor=mid, edgecolor="white", hatch="////"),
        plt.Rectangle((0, 0), 1, 1, facecolor=light, edgecolor="white"),
    ]
    labels = [
        r"Electricity ($E\,\beta\,n$)",
        r"Leaks ($\mathrm{GWP}\,m\,L\,n$)",
        r"End of life ($\mathrm{GWP}\,m\,(1-\alpha)$)",
    ]
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.995))
    fig.tight_layout(rect=(0, 0, 1, 0.94))

    fig.savefig(os.path.join(OUT_DIR, f"tewi_figure_{tag}.png"), dpi=300)
    fig.savefig(os.path.join(OUT_DIR, f"tewi_figure_{tag}.pdf"))
    plt.close(fig)


# =====================================================================
# 5. RUN
# =====================================================================

if __name__ == "__main__":
    os.makedirs(OUT_DIR, exist_ok=True)

    # tag goes into the file names, so different settings never overwrite each other
    tag = f"L{int(round(LEAK_RATE * 100))}_n{LIFETIME_YR}_cascade{int(round(CASCADE_HOT_CHARGE_FRACTION * 100))}"

    electricity = load_electricity()
    results = pd.concat([electricity, electricity.apply(tewi_row, axis=1)], axis=1)

    out_csv = os.path.join(OUT_DIR, f"tewi_results_{tag}.csv")
    results.to_csv(out_csv, index=False, float_format="%.3f")
    print(f"Saved {out_csv} ({len(results)} rows)")

    make_figure(results, tag)
    print(f"Saved {OUT_DIR}/tewi_figure_{tag}.png and .pdf")

    # quick look: Chicago, single stage R134a, no limit
    check = results[(results["Case"] == "No limit") & (results["City"] == "Chicago, IL")
                    & (results["Cycle"] == "Single stage") & (results["Refrigerant"] == "R134a")]
    print(check[["leaks_t", "end_of_life_t", "electricity_t", "TEWI_t"]].round(1).to_string(index=False))

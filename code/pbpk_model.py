"""
================================================================================
LAYER 2: 7-COMPARTMENT PBPK  — COMPATIBILITY SHIM
================================================================================
The model now lives in `core.py` (the single source of truth). This module
re-exports it so existing imports keep working. DO NOT redefine the model here.
================================================================================
"""

from __future__ import annotations

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from core import (  # noqa: F401  (re-export)
    get_physiology,
    calibrate_kp,
    pbpk_odes,
    extract_pk,
    extract_pk_pbpk,
    recommend_dose,
    compute_dose_recommendation,
    simulate_pbpk,
    run_pbpk_analysis,
)


if __name__ == '__main__':
    from pathlib import Path
    from drug_properties import analyze_drug

    OUT = Path(__file__).resolve().parent.parent / 'figures'
    OUT.mkdir(exist_ok=True)

    print('VALIDATING 7-COMPARTMENT PBPK (imports from core.py)')
    # Anchors are dose-NORMALISED to 500 mg (see clean_dataset.py):
    #   Earth Cmax ≈ 7.1 µg/mL (pooled), Kovachevich single-study 500 mg = 5.13
    lit_cmax = {500: 5.13, 1000: 9.41}
    for dose in (500, 1000):
        drug = analyze_drug('Paracetamol', mission_days=30, dose_mg=dose, verbose=False)
        r = run_pbpk_analysis(drug, mission_days=30, verbose=True)
        lit = lit_cmax[dose]
        err = (r['pk_earth']['Cmax'] - lit) / lit * 100
        print(f'  vs single-study lit Cmax≈{lit} ({dose} mg): error {err:+.1f}%')

    drug = analyze_drug('Paracetamol', mission_days=30, verbose=False)
    result = run_pbpk_analysis(drug, mission_days=30, verbose=False)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.patch.set_facecolor('#0a0a1a')
    titles = ['Arterial (systemic)', 'Liver', 'Peripheral tissue']
    pairs = [
        (result['C_earth'], result['C_space']),
        (result['C_liver_e'], result['C_liver_s']),
        (result['C_tissue_e'], result['C_tissue_s']),
    ]
    for ax, title, (ce, cs) in zip(axes, titles, pairs):
        ax.set_facecolor('#0d1117')
        ax.plot(result['t'], ce, '#2196F3', lw=2.5, label='Earth')
        ax.plot(result['t'], cs, '#F44336', lw=2.5, label='Space')
        ax.set(xlabel='Time (h)', ylabel='Conc (mg/L)', title=title)
        ax.legend(facecolor='#1a1a2e', labelcolor='white', fontsize=9)
        ax.tick_params(colors='white')
        for lbl in (ax.xaxis.label, ax.yaxis.label, ax.title):
            lbl.set_color('white')
        ax.grid(True, alpha=0.2)
    plt.suptitle('7-compartment PBPK (Vd-calibrated) — Paracetamol Earth vs Space',
                 color='white', fontweight='bold')
    plt.tight_layout()
    plt.savefig(OUT / 'pbpk_layer2.png', dpi=150, bbox_inches='tight', facecolor='#0a0a1a')
    print('figure saved')

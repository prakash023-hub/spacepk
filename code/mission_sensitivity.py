"""
Mission-phase sensitivity analysis — SpacePK Framework
Tracks Cmax, AUC, and dose adjustment across mission timeline.
"""

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from drug_properties import DRUG_LIST, analyze_drug
from iss_drug_catalog import ISS_DRUG_CATALOG
from pbpk_model import compute_dose_recommendation, run_pbpk_analysis

MISSION_DAYS = [0, 1, 3, 7, 14, 30, 60, 90, 180]
DEMO_DRUGS = [d for d in DRUG_LIST if ISS_DRUG_CATALOG.get(d, {}).get('space_data')]
if len(DEMO_DRUGS) < 4:
    DEMO_DRUGS = ['Paracetamol', 'Ibuprofen', 'Ciprofloxacin', 'Promethazine']


def run_mission_sensitivity(drugs=None, body_weight=75, out_dir=None):
    drugs = drugs or DEMO_DRUGS
    if out_dir is None:
        out_dir = Path(__file__).resolve().parent.parent / 'figures'
    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True)

    print('SPACE PK — MISSION PHASE SENSITIVITY ANALYSIS')
    print('=' * 60)

    all_results = {}
    for drug in drugs:
        if drug not in DRUG_LIST:
            continue
        rows = []
        for day in MISSION_DAYS:
            profile = analyze_drug(drug, mission_days=day, body_weight_kg=body_weight, verbose=False)
            sim = run_pbpk_analysis(profile, mission_days=day, body_weight=body_weight, verbose=False)
            pk_e = sim['pk_earth']
            pk_s = sim['pk_space']
            dose = compute_dose_recommendation(profile, pk_e, pk_s)
            rows.append({
                'mission_day': day,
                'phase': profile['space_phase'],
                'Cmax_earth': pk_e['Cmax'],
                'Cmax_space': pk_s['Cmax'],
                'AUC_earth': pk_e['AUC'],
                'AUC_space': pk_s['AUC'],
                'dose_factor': dose['adjustment_factor'] if dose else 1.0,
                'space_dose_mg': dose['space_dose_mg'] if dose else profile['dose_mg'],
            })
        all_results[drug] = rows
        print(f"  {drug}: day-180 dose factor ×{rows[-1]['dose_factor']:.2f}")

    _plot_sensitivity(all_results, out_dir)
    _plot_heatmap(all_results, out_dir)
    print(f"\n✅ Sensitivity plots saved to {out_dir}")
    return all_results


def _plot_sensitivity(all_results, out_dir):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.patch.set_facecolor('#0a0a1a')
    colors = ['#2196F3', '#F44336', '#4CAF50', '#FF9800', '#9C27B0', '#00BCD4']

    for ax, (drug, rows) in zip(axes.flat, all_results.items()):
        ax.set_facecolor('#0d1117')
        days = [r['mission_day'] for r in rows]
        cmax_chg = [(r['Cmax_space'] - r['Cmax_earth']) / r['Cmax_earth'] * 100
                    if r['Cmax_earth'] > 0 else 0 for r in rows]
        dose_fac = [r['dose_factor'] for r in rows]

        ax2 = ax.twinx()
        ax.plot(days, cmax_chg, 'o-', color=colors[0], lw=2, label='Cmax change (%)')
        ax2.plot(days, dose_fac, 's--', color=colors[1], lw=2, label='Dose factor')
        ax.axhline(0, color='white', alpha=0.2, lw=1)
        ax.axvline(3, color='#FFC107', alpha=0.4, ls=':', label='Acute→Adapt')
        ax.axvline(14, color='#FF5722', alpha=0.4, ls=':', label='Adapt→Chronic')

        ax.set(xlabel='Mission day', ylabel='Cmax change vs Earth (%)', title=drug)
        ax2.set_ylabel('Recommended dose factor', color=colors[1])
        ax.tick_params(colors='white')
        ax2.tick_params(colors=colors[1])
        ax.xaxis.label.set_color('white')
        ax.yaxis.label.set_color('white')
        ax.title.set_color('white')
        ax.grid(True, alpha=0.2)

    plt.suptitle('Mission-Phase PK Sensitivity — Earth vs Spaceflight',
                 color='white', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(out_dir / 'mission_sensitivity.png', dpi=150, bbox_inches='tight', facecolor='#0a0a1a')
    plt.close()


def _plot_heatmap(all_results, out_dir):
    drugs = list(all_results.keys())
    matrix = np.array([[r['dose_factor'] for r in all_results[d]] for d in drugs])

    fig, ax = plt.subplots(figsize=(12, 5))
    fig.patch.set_facecolor('#0a0a1a')
    ax.set_facecolor('#0d1117')
    im = ax.imshow(matrix, aspect='auto', cmap='RdYlGn_r', vmin=0.8, vmax=1.6)
    ax.set_xticks(range(len(MISSION_DAYS)))
    ax.set_xticklabels(MISSION_DAYS, color='white')
    ax.set_yticks(range(len(drugs)))
    ax.set_yticklabels(drugs, color='white')
    ax.set_xlabel('Mission day', color='white')
    ax.set_title('Space Dose Adjustment Factor (Cmax-matched)', color='white', fontweight='bold')

    for i in range(len(drugs)):
        for j in range(len(MISSION_DAYS)):
            ax.text(j, i, f'{matrix[i, j]:.2f}', ha='center', va='center',
                    color='white', fontsize=9, fontweight='bold')

    cbar = plt.colorbar(im, ax=ax)
    cbar.ax.yaxis.set_tick_params(color='white')
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')
    plt.tight_layout()
    plt.savefig(out_dir / 'mission_dose_heatmap.png', dpi=150, bbox_inches='tight', facecolor='#0a0a1a')
    plt.close()


if __name__ == '__main__':
    run_mission_sensitivity()

"""
SpacePK — Interactive Flight Surgeon Decision Support
Streamlit app: cheminformatics → PBPK → dose recommendation
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import streamlit as st

CODE_DIR = Path(__file__).resolve().parent / 'code'
sys.path.insert(0, str(CODE_DIR))

from drug_properties import DRUG_LIST, analyze_drug  # noqa: E402
from pbpk_model import run_pbpk_analysis  # noqa: E402

st.set_page_config(
    page_title='SpacePK',
    page_icon='🛸',
    layout='wide',
    initial_sidebar_state='expanded',
)

st.markdown("""
<style>
    .main { background-color: #0a0a1a; }
    .stMetric { background: #0d1117; padding: 12px; border-radius: 8px; }
    h1, h2, h3, p, label { color: #e6edf3 !important; }
</style>
""", unsafe_allow_html=True)

st.title('🛸 SpacePK — Spaceflight Pharmacokinetics Framework')
st.caption('Integrated RDKit → 7-compartment PBPK → Bayesian PopPK | NASA-relevant dose adjustment')

with st.sidebar:
    st.header('Mission parameters')
    drug = st.selectbox('Drug', DRUG_LIST, index=0)
    mission_days = st.slider('Mission day', 0, 180, 30, help='0 = Earth baseline')
    body_weight = st.number_input('Body weight (kg)', 50, 120, 75)
    st.divider()
    st.markdown('**Mission phases**')
    if mission_days == 0:
        st.info('🌍 Earth baseline')
    elif mission_days <= 3:
        st.warning('⚡ Acute — fluid shift, faster absorption')
    elif mission_days <= 14:
        st.info('🔄 Adaptation — mixed physiology')
    else:
        st.success('🛰 Chronic — delayed absorption, reduced plasma volume')

profile = analyze_drug(drug, mission_days=mission_days, body_weight_kg=body_weight, verbose=False)
result = run_pbpk_analysis(profile, mission_days=mission_days, body_weight=body_weight, verbose=False)
pk_e, pk_s = result['pk_earth'], result['pk_space']
dose = result['dose_recommendation']

tab1, tab2, tab3, tab4 = st.tabs(['📊 Properties', '📈 PBPK Curves', '💊 Dose', '🔬 Research'])

with tab1:
    c1, c2, c3, c4 = st.columns(4)
    c1.metric('MW', f"{profile['MW']:.1f} g/mol")
    c2.metric('LogP', f"{profile['logP']:.2f}")
    c3.metric('BCS Class', profile['BCS'])
    c4.metric('QED', f"{profile['QED']:.3f}")

    c5, c6, c7, c8 = st.columns(4)
    c5.metric('F (Earth)', f"{profile['F_earth']:.2f}")
    c6.metric('F (Space)', f"{profile['F_space']:.2f}")
    c7.metric('Vd Earth', f"{profile['Vd_earth']:.1f} L")
    c8.metric('Vd Space', f"{profile['Vd_space']:.1f} L")

    mods = profile['space_modifiers']
    st.subheader(f"Space modifiers — {profile['space_phase']} phase (day {mission_days})")
    m1, m2, m3, m4 = st.columns(4)
    m1.metric('ka factor', f"×{mods['ka_factor']}")
    m2.metric('F factor', f"×{mods['F_factor']}")
    m3.metric('Vd factor', f"×{mods['Vd_factor']}")
    m4.metric('ke factor', f"×{mods['ke_factor']}")

with tab2:
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    fig.patch.set_facecolor('#0a0a1a')
    titles = ['Systemic (Venous)', 'Liver', 'Peripheral Tissue']
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
        ax.legend(facecolor='#1a1a2e', labelcolor='white', fontsize=8)
        ax.tick_params(colors='white')
        ax.grid(True, alpha=0.2)
    plt.suptitle(f'{drug} — Mission day {mission_days}', color='white', fontweight='bold')
    plt.tight_layout()
    st.pyplot(fig)
    plt.close()

    st.subheader('PK comparison')
    st.dataframe({
        'Parameter': ['Cmax (mg/L)', 'Tmax (h)', 'AUC (mg·h/L)', 't½ (h)'],
        'Earth': [pk_e['Cmax'], pk_e['Tmax'], pk_e['AUC'], pk_e['t12']],
        'Space': [pk_s['Cmax'], pk_s['Tmax'], pk_s['AUC'], pk_s['t12']],
    }, use_container_width=True, hide_index=True)

with tab3:
    if dose:
        st.subheader('Flight dose recommendation')
        d1, d2, d3 = st.columns(3)
        d1.metric('Earth dose', f"{dose['earth_dose_mg']:.0f} mg")
        d2.metric('Recommended space dose', f"{dose['space_dose_mg']:.0f} mg")
        d3.metric('Adjustment', f"{dose['adjustment_pct']:+.1f}%")
        st.info(dose['rationale'])
        if dose['adjustment_factor'] > 1.15:
            st.warning('⚠️ >15% dose increase suggested — monitor for toxicity at peak exposure.')
        elif dose['adjustment_factor'] < 0.85:
            st.warning('⚠️ Dose reduction suggested — verify efficacy at trough concentrations.')
    else:
        st.error('Could not compute dose recommendation for this profile.')

with tab4:
    st.markdown("""
    ### Novelty — SpacePK integrated framework

    **What is new in this work:**
    - Single pipeline linking **RDKit cheminformatics → mechanistic PBPK → Bayesian PopPK**
    - **Mission-phase modifiers** (acute / adaptation / chronic) from HDT and spaceflight literature
    - **Multi-drug extensibility** for ISS-relevant compounds
    - **Actionable dose output** with uncertainty (Bayesian layer)

    **Key references:** Gandia 2003, Kovachevich 2009, Polyakov 2021, Leach 1981, Dello Russo 2022

    **Run full pipeline:**
    ```bash
    cd ~/spacepk/code
    python drug_properties.py
    python pbpk_model.py
    python bayesian_popPK.py
    python mission_sensitivity.py
    ```

    **GitHub:** [prakash023-hub/spacepk](https://github.com/prakash023-hub/spacepk)
    """)

    if st.button('Generate mission sensitivity plots (all drugs)'):
        with st.spinner('Running sensitivity analysis…'):
            from mission_sensitivity import run_mission_sensitivity
            run_mission_sensitivity()
        st.success('Saved to figures/mission_sensitivity.png and mission_dose_heatmap.png')
        st.image(str(Path(__file__).parent / 'figures' / 'mission_sensitivity.png'))

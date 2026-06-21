"""
SpacePK — Interactive Flight Surgeon Decision Support
Full ISS medical kit catalog + PBPK simulation
"""

import sys
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit as st

CODE_DIR = Path(__file__).resolve().parent / 'code'
sys.path.insert(0, str(CODE_DIR))

from drug_properties import analyze_drug  # noqa: E402
from iss_drug_catalog import (  # noqa: E402
    ISS_DRUG_CATALOG,
    get_catalog_table,
    get_categories,
    search_drugs,
)
from pbpk_model import run_pbpk_analysis  # noqa: E402

# ── Page config ───────────────────────────────────────────────────────────────
st.set_page_config(
    page_title='SpacePK | ISS Pharmacokinetics',
    page_icon='🛸',
    layout='wide',
    initial_sidebar_state='expanded',
)

st.markdown("""
<style>
    .block-container { padding-top: 1.2rem; max-width: 1400px; }
    [data-testid="stSidebar"] { background: linear-gradient(180deg, #0d1117 0%, #161b22 100%); }
    .hero {
        background: linear-gradient(135deg, #1a237e 0%, #0d47a1 50%, #1565c0 100%);
        padding: 1.5rem 2rem; border-radius: 12px; margin-bottom: 1.5rem;
        border: 1px solid #30363d;
    }
    .hero h1 { color: #fff !important; margin: 0; font-size: 1.8rem; }
    .hero p { color: #b3d4fc !important; margin: 0.4rem 0 0 0; }
    div[data-testid="stMetric"] {
        background: #161b22; border: 1px solid #30363d; border-radius: 10px;
        padding: 0.75rem 1rem;
    }
    .stTabs [data-baseweb="tab-list"] { gap: 8px; }
    .stTabs [data-baseweb="tab"] {
        background: #21262d; border-radius: 8px; padding: 8px 16px;
        border: 1px solid #30363d;
    }
</style>
""", unsafe_allow_html=True)


@st.cache_data(show_spinner=False)
def simulate(drug_name, mission_days, body_weight, dose_mg=None):
    profile = analyze_drug(
        drug_name, mission_days=mission_days,
        body_weight_kg=body_weight, dose_mg=dose_mg, verbose=False,
    )
    result = run_pbpk_analysis(
        profile, mission_days=mission_days,
        body_weight=body_weight, verbose=False,
    )
    return profile, result


def mission_badge(days):
    if days == 0:
        return '🌍 Earth baseline', 'info'
    if days <= 3:
        return '⚡ Acute phase (fluid shift)', 'warning'
    if days <= 14:
        return '🔄 Adaptation phase', 'info'
    return '🛰 Chronic phase (long-duration)', 'success'


def plot_pk_curves(result, drug, mission_days):
    t = result['t']
    fig = make_subplots(
        rows=1, cols=3,
        subplot_titles=('Systemic blood', 'Liver', 'Peripheral tissue'),
    )
    traces = [
        (result['C_earth'], result['C_space']),
        (result['C_liver_e'], result['C_liver_s']),
        (result['C_tissue_e'], result['C_tissue_s']),
    ]
    for col, (ce, cs) in enumerate(traces, 1):
        fig.add_trace(
            go.Scatter(x=t, y=ce, name='Earth', line=dict(color='#42a5f5', width=2.5),
                       showlegend=(col == 1)),
            row=1, col=col,
        )
        fig.add_trace(
            go.Scatter(x=t, y=cs, name='Space', line=dict(color='#ef5350', width=2.5),
                       showlegend=(col == 1)),
            row=1, col=col,
        )
    fig.update_layout(
        height=380,
        template='plotly_dark',
        paper_bgcolor='#0d1117',
        plot_bgcolor='#161b22',
        title=dict(text=f'{drug} — Mission day {mission_days}', x=0.5),
        legend=dict(orientation='h', y=1.12),
        margin=dict(t=60, b=40),
    )
    fig.update_xaxes(title_text='Time (h)')
    fig.update_yaxes(title_text='Conc (mg/L)')
    return fig


def plot_dose_timeline(drug, body_weight):
    days = [0, 1, 3, 7, 14, 30, 60, 90, 180]
    factors = []
    for d in days:
        _, res = simulate(drug, d, body_weight)
        dr = res.get('dose_recommendation')
        factors.append(dr['adjustment_factor'] if dr else 1.0)
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=days, y=factors, mode='lines+markers',
        line=dict(color='#66bb6a', width=3),
        marker=dict(size=8),
        name='Dose factor',
    ))
    fig.add_hline(y=1.0, line_dash='dash', line_color='#888')
    fig.update_layout(
        height=320, template='plotly_dark',
        paper_bgcolor='#0d1117', plot_bgcolor='#161b22',
        title=f'{drug} — Space dose factor vs mission day',
        xaxis_title='Mission day', yaxis_title='Recommended dose factor (×Earth)',
    )
    return fig


# ── Header ────────────────────────────────────────────────────────────────────
n_drugs = len(ISS_DRUG_CATALOG)
n_iss = sum(1 for e in ISS_DRUG_CATALOG.values() if e.get('iss_kit'))
n_space = sum(1 for e in ISS_DRUG_CATALOG.values() if e.get('space_data'))

st.markdown(f"""
<div class="hero">
  <h1>🛸 SpacePK — Spaceflight Pharmacokinetics</h1>
  <p>{n_drugs} drugs · {n_iss} ISS medical kit · {n_space} with spaceflight PK literature
  · RDKit → 7-compartment PBPK → dose recommendation</p>
</div>
""", unsafe_allow_html=True)

# ── Sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header('🔍 Find a drug')
    search_q = st.text_input('Search by name', placeholder='e.g. ibuprofen, cipro…')
    category = st.selectbox('Category', ['All'] + get_categories())
    iss_only = st.checkbox('ISS medical kit only', value=False)
    space_only = st.checkbox('Space PK literature only', value=False)

    filtered = search_drugs(
        query=search_q,
        category=None if category == 'All' else category,
        iss_only=iss_only,
        space_data_only=space_only,
    )

    if not filtered:
        st.warning('No drugs match filters.')
        st.stop()

    drug = st.selectbox(
        f'Drug ({len(filtered)} available)',
        filtered,
        index=0,
    )

    meta = ISS_DRUG_CATALOG.get(drug, {})
    st.caption(f"**{meta.get('category', '')}** · "
               f"{'ISS kit ✅' if meta.get('iss_kit') else 'General'} · "
               f"{'Space data ✅' if meta.get('space_data') else 'Estimated PK'}")

    st.divider()
    st.header('🚀 Mission')
    mission_days = st.slider('Mission day', 0, 180, 30, help='0 = Earth control')
    body_weight = st.number_input('Body weight (kg)', 50, 120, 75, step=1)

    default_dose = meta.get('dose_mg', 500)
    custom_dose = st.number_input('Dose (mg)', 0.1, 5000.0, float(default_dose), step=0.1)

    badge, level = mission_badge(mission_days)
    getattr(st, level)(badge)

# ── Run simulation ────────────────────────────────────────────────────────────
profile, result = simulate(drug, mission_days, body_weight, dose_mg=custom_dose)
pk_e, pk_s = result['pk_earth'], result['pk_space']
dose_rec = result['dose_recommendation']

# ── Top metrics ───────────────────────────────────────────────────────────────
m1, m2, m3, m4, m5, m6 = st.columns(6)
m1.metric('BCS class', profile['BCS'])
m2.metric('Cmax Earth', f"{pk_e['Cmax']:.3f}", f"Space {pk_s['Cmax']:.3f}")
m3.metric('AUC Earth', f"{pk_e['AUC']:.2f}", f"Space {pk_s['AUC']:.2f}")
m4.metric('F (Earth→Space)', f"{profile['F_earth']:.2f}→{profile['F_space']:.2f}")
if dose_rec:
    m5.metric('Space dose', f"{dose_rec['space_dose_mg']:.0f} mg", f"{dose_rec['adjustment_pct']:+.0f}%")
    m6.metric('Factor', f"×{dose_rec['adjustment_factor']:.2f}")
else:
    m5.metric('Space dose', 'N/A')
    m6.metric('Factor', '—')

# ── Tabs ──────────────────────────────────────────────────────────────────────
tab_sim, tab_dose, tab_props, tab_catalog, tab_custom = st.tabs([
    '📈 Simulation', '💊 Dose plan', '🧬 Properties', '📋 Full drug list', '⚗️ Custom SMILES',
])

with tab_sim:
    st.plotly_chart(plot_pk_curves(result, drug, mission_days), use_container_width=True)

    c1, c2 = st.columns(2)
    with c1:
        st.subheader('PK comparison')
        df_pk = pd.DataFrame({
            'Parameter': ['Cmax (mg/L)', 'Tmax (h)', 'AUC (mg·h/L)', 't½ (h)'],
            'Earth': [pk_e['Cmax'], pk_e['Tmax'], pk_e['AUC'], pk_e['t12']],
            'Space': [pk_s['Cmax'], pk_s['Tmax'], pk_s['AUC'], pk_s['t12']],
        })
        df_pk['Change %'] = [
            f"{((pk_s[k] - pk_e[k]) / pk_e[k] * 100):+.1f}" if pk_e[k] else '—'
            for k in ['Cmax', 'Tmax', 'AUC', 't12']
            if isinstance(pk_e.get('Cmax'), (int, float))
        ]
        # simpler change column
        changes = []
        for param in ['Cmax', 'Tmax', 'AUC', 't12']:
            ve, vs = pk_e[param], pk_s[param]
            if isinstance(ve, (int, float)) and isinstance(vs, (int, float)) and ve:
                changes.append(f"{((vs - ve) / ve) * 100:+.1f}%")
            else:
                changes.append('—')
        df_pk['Change %'] = changes
        st.dataframe(df_pk, hide_index=True, use_container_width=True)

    with c2:
        st.subheader('Mission timeline')
        st.plotly_chart(plot_dose_timeline(drug, body_weight), use_container_width=True)

with tab_dose:
    if dose_rec:
        st.subheader('Flight surgeon dose recommendation')
        d1, d2, d3, d4 = st.columns(4)
        d1.metric('Standard Earth dose', f"{dose_rec['earth_dose_mg']:.0f} mg")
        d2.metric('Recommended space dose', f"{dose_rec['space_dose_mg']:.0f} mg")
        d3.metric('AUC-matched dose', f"{dose_rec.get('space_dose_auc_mg', 0):.0f} mg")
        d4.metric('Adjustment', f"{dose_rec['adjustment_pct']:+.1f}%")

        if dose_rec['adjustment_factor'] > 1.15:
            st.warning('⚠️ Increase >15% — monitor peak exposure and hepatotoxicity limits.')
        elif dose_rec['adjustment_factor'] < 0.85:
            st.warning('⚠️ Reduction suggested — confirm efficacy at trough concentrations.')
        else:
            st.success('✅ Dose change within ±15% — standard Earth dosing may be acceptable with monitoring.')

        st.info(
            f"**Rationale:** {dose_rec['rationale']}. "
            f"Phase: **{profile['space_phase']}** (day {mission_days}). "
            f"Based on 7-compartment PBPK with literature space modifiers "
            f"(Gandia 2003, Kovachevich 2009, Polyakov 2021)."
        )
    else:
        st.error('Dose recommendation unavailable (IV bolus or zero exposure).')

    mods = profile['space_modifiers']
    st.subheader('Applied space modifiers')
    mc1, mc2, mc3, mc4 = st.columns(4)
    mc1.metric('Absorption (ka)', f"×{mods['ka_factor']}")
    mc2.metric('Bioavailability (F)', f"×{mods['F_factor']}")
    mc3.metric('Distribution (Vd)', f"×{mods['Vd_factor']}")
    mc4.metric('Elimination (ke)', f"×{mods['ke_factor']}")

with tab_props:
    p1, p2, p3, p4 = st.columns(4)
    p1.metric('MW', f"{profile['MW']:.1f} g/mol")
    p2.metric('LogP', f"{profile['logP']:.2f}")
    p3.metric('PSA', f"{profile['PSA']:.1f} Å²")
    p4.metric('QED', f"{profile['QED']:.3f}")

    p5, p6, p7, p8 = st.columns(4)
    p5.metric('HBD / HBA', f"{profile['HBD']} / {profile['HBA']}")
    p6.metric('Protein binding', f"{profile['protein_binding']*100:.0f}%")
    p7.metric('Vd Earth', f"{profile['Vd_earth']:.1f} L")
    p8.metric('Vd Space', f"{profile['Vd_space']:.1f} L")

    st.code(meta.get('smiles', profile.get('smiles', '')), language=None)
    st.caption(f"CYP: {meta.get('cyp', '—')} · Transporter: {meta.get('transporter', '—')}")

with tab_catalog:
    st.subheader(f'Complete catalog — {n_drugs} compounds')
    df_cat = pd.DataFrame(get_catalog_table())
    st.dataframe(
        df_cat,
        hide_index=True,
        use_container_width=True,
        height=480,
        column_config={
            'Dose (mg)': st.column_config.NumberColumn(format='%.1f'),
            't½ (h)': st.column_config.NumberColumn(format='%.1f'),
        },
    )
    st.download_button(
        'Download CSV',
        df_cat.to_csv(index=False),
        file_name='spacepk_drug_catalog.csv',
        mime='text/csv',
    )

with tab_custom:
    st.subheader('Analyze any molecule by SMILES')
    cname = st.text_input('Compound name', value='Custom drug')
    smiles = st.text_input('SMILES string', placeholder='e.g. CC(=O)Nc1ccc(O)cc1')
    cdose = st.number_input('Dose mg (custom)', 0.1, 5000.0, 500.0)
    if st.button('Run custom analysis', type='primary'):
        if not smiles.strip():
            st.error('Enter a valid SMILES string.')
        else:
            try:
                cp = analyze_drug(
                    cname, smiles=smiles.strip(), dose_mg=cdose,
                    mission_days=mission_days, body_weight_kg=body_weight, verbose=False,
                )
                cr = run_pbpk_analysis(
                    cp, mission_days=mission_days,
                    body_weight=body_weight, verbose=False,
                )
                st.plotly_chart(plot_pk_curves(cr, cname, mission_days), use_container_width=True)
                if cr.get('dose_recommendation'):
                    dr = cr['dose_recommendation']
                    st.success(
                        f"Space dose: **{dr['space_dose_mg']:.0f} mg** "
                        f"(×{dr['adjustment_factor']:.2f})"
                    )
            except Exception as exc:
                st.error(f'Analysis failed: {exc}')

st.divider()
st.caption('SpacePK · [GitHub](https://github.com/prakash023-hub/spacepk) · CPT:PSP research framework')

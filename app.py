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

from drug_properties import analyze_drug, HAS_RDKIT  # noqa: E402
from iss_drug_catalog import (  # noqa: E402
    ISS_DRUG_CATALOG,
    get_catalog_table,
    get_categories,
    search_drugs,
)
from pbpk_model import run_pbpk_analysis  # noqa: E402

st.set_page_config(
    page_title='SpacePK | ISS Pharmacokinetics',
    page_icon='🛸',
    layout='wide',
    initial_sidebar_state='expanded',
)

st.markdown("""
<style>
    .block-container { padding-top: 1rem; max-width: 1280px; }
    div[data-testid="stMetric"] {
        background: #161b22; border: 1px solid #30363d;
        border-radius: 8px; padding: 0.5rem 0.75rem;
    }
    .hero-box {
        background: linear-gradient(120deg, #1565c0 0%, #0d47a1 100%);
        border: 1px solid #30363d; border-radius: 10px;
        padding: 1.25rem 1.5rem; margin-bottom: 1rem;
    }
    .hero-box h2 { color: #fff; margin: 0; font-size: 1.5rem; }
    .hero-box p { color: #cfe8ff; margin: 0.35rem 0 0 0; font-size: 0.95rem; }
</style>
""", unsafe_allow_html=True)


@st.cache_data(show_spinner='Running PBPK simulation…')
def simulate(drug_name, mission_days, body_weight, dose_mg):
    profile = analyze_drug(
        drug_name, mission_days=mission_days,
        body_weight_kg=body_weight, dose_mg=dose_mg, verbose=False,
    )
    result = run_pbpk_analysis(
        profile, mission_days=mission_days,
        body_weight=body_weight, verbose=False,
    )
    return profile, result


@st.cache_data(show_spinner='Building mission timeline…')
def dose_timeline_data(drug_name, body_weight, dose_mg):
    days = [0, 1, 3, 7, 14, 30, 60, 90, 180]
    factors = []
    for d in days:
        _, res = simulate(drug_name, d, body_weight, dose_mg)
        dr = res.get('dose_recommendation')
        factors.append(dr['adjustment_factor'] if dr else 1.0)
    return days, factors


def pk_change_pct(earth_val, space_val):
    if isinstance(earth_val, (int, float)) and isinstance(space_val, (int, float)) and earth_val:
        return f"{((space_val - earth_val) / earth_val) * 100:+.1f}%"
    return '—'


def plot_pk_curves(result, drug, mission_days):
    t = result['t']
    fig = make_subplots(
        rows=1, cols=3,
        subplot_titles=('Systemic blood', 'Liver', 'Peripheral tissue'),
        horizontal_spacing=0.08,
    )
    pairs = [
        (result['C_earth'], result['C_space']),
        (result['C_liver_e'], result['C_liver_s']),
        (result['C_tissue_e'], result['C_tissue_s']),
    ]
    for col, (ce, cs) in enumerate(pairs, 1):
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
        height=400,
        template='plotly_dark',
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='#161b22',
        title=dict(text=f'{drug} — Mission day {mission_days}', x=0.5, font=dict(size=16)),
        legend=dict(orientation='h', yanchor='bottom', y=1.02, x=0.5, xanchor='center'),
        margin=dict(t=70, b=40, l=50, r=20),
    )
    fig.update_xaxes(title_text='Time (h)', gridcolor='#30363d')
    fig.update_yaxes(title_text='Concentration (mg/L)', gridcolor='#30363d')
    return fig


def plot_dose_timeline(drug, body_weight, dose_mg):
    days, factors = dose_timeline_data(drug, body_weight, dose_mg)
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=days, y=factors, mode='lines+markers',
        line=dict(color='#66bb6a', width=3),
        marker=dict(size=9, color='#81c784'),
        hovertemplate='Day %{x}<br>Factor ×%{y:.2f}<extra></extra>',
    ))
    fig.add_hline(y=1.0, line_dash='dash', line_color='#888', annotation_text='No change')
    fig.add_vrect(x0=0, x1=3, fillcolor='#ff9800', opacity=0.08, line_width=0)
    fig.add_vrect(x0=3, x1=14, fillcolor='#2196f3', opacity=0.08, line_width=0)
    fig.add_vrect(x0=14, x1=180, fillcolor='#4caf50', opacity=0.08, line_width=0)
    fig.update_layout(
        height=340, template='plotly_dark',
        paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='#161b22',
        title=dict(text=f'{drug} — dose factor vs mission day', x=0.5),
        xaxis_title='Mission day', yaxis_title='Space dose factor (× Earth dose)',
        margin=dict(t=50, b=40),
    )
    return fig


# ── Stats ─────────────────────────────────────────────────────────────────────
n_drugs = len(ISS_DRUG_CATALOG)
n_iss = sum(1 for e in ISS_DRUG_CATALOG.values() if e.get('iss_kit'))
n_space = sum(1 for e in ISS_DRUG_CATALOG.values() if e.get('space_data'))

st.title('SpacePK — Spaceflight Pharmacokinetics')
st.markdown(
    f'**{n_drugs} drugs** · **{n_iss} ISS medical kit** · '
    f'**{n_space} with spaceflight PK literature** · '
    'RDKit → PBPK → dose recommendation'
)

with st.expander('URLs for publication / Devpost / paper (important)'):
    st.markdown("""
**Do not submit** (only work on your Mac):
- `http://localhost:8501`
- `http://192.168.x.x:8501`
- `http://117.254.x.x:8501` (your IP — stops when laptop is off)

**Submit these:**
| What | URL |
|------|-----|
| **Code & reproducibility** | https://github.com/prakash023-hub/spacepk |
| **Live demo** | Deploy at [share.streamlit.io](https://share.streamlit.io) → `https://YOUR-NAME-spacepk.streamlit.app` |
| **Paper** | GitHub repo + figure files in `figures/` |

**Deploy in 5 min:** GitHub → share.streamlit.io → New app → repo `prakash023-hub/spacepk` → main file `app.py` → Advanced → Python version 3.12 → use `environment.yml`.
    """)

# ── Sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    if not HAS_RDKIT:
        st.caption('☁️ Cloud mode — 41 catalog drugs (precomputed descriptors)')
    st.header('Find a drug')
    search_q = st.text_input('Search name or category', placeholder='ibuprofen, antibiotic…')
    category = st.selectbox('Category', ['All categories'] + get_categories())
    col_f1, col_f2 = st.columns(2)
    with col_f1:
        iss_only = st.checkbox('ISS kit only', value=False)
    with col_f2:
        space_only = st.checkbox('Space data only', value=False)

    filtered = search_drugs(
        query=search_q,
        category=None if category == 'All categories' else category,
        iss_only=iss_only,
        space_data_only=space_only,
    )

    if not filtered:
        st.warning('No drugs match. Clear filters or search.')
        st.stop()

    drug = st.selectbox(f'Select drug ({len(filtered)})', filtered)
    meta = ISS_DRUG_CATALOG[drug]

    st.markdown(
        f"**{meta['category']}**  \n"
        f"{'ISS medical kit' if meta.get('iss_kit') else 'Supplemental'} · "
        f"{'Space PK literature' if meta.get('space_data') else 'Estimated PK'}"
    )

    st.divider()
    st.header('Mission settings')
    mission_days = st.slider('Mission day', 0, 180, 30, help='0 = Earth baseline')

    if mission_days == 0:
        st.info('Earth baseline — no space modifiers')
    elif mission_days <= 3:
        st.warning('Acute phase — fluid shift')
    elif mission_days <= 14:
        st.info('Adaptation phase')
    else:
        st.success('Chronic phase — long-duration mission')

    body_weight = st.number_input('Body weight (kg)', min_value=50, max_value=120, value=75, step=1)

    default_dose = float(meta.get('dose_mg', 500))
    dose_key = f'dose_{drug.replace(" ", "_")}'
    if dose_key not in st.session_state:
        st.session_state[dose_key] = default_dose

    dose_mg = st.number_input(
        'Dose (mg)', min_value=0.1, max_value=5000.0,
        step=0.5, key=dose_key,
    )

    route = 'IV bolus' if meta.get('lit_ka_h') == 'IV' else 'Oral'
    st.caption(f"Route: **{route}** · Literature default: **{default_dose:g} mg**")

# ── Simulation ────────────────────────────────────────────────────────────────
profile, result = simulate(drug, mission_days, body_weight, dose_mg)
pk_e, pk_s = result['pk_earth'], result['pk_space']
dose_rec = result['dose_recommendation']

# Header metrics
c1, c2, c3, c4, c5, c6 = st.columns(6)
c1.metric('Drug', drug)
c2.metric('BCS class', str(profile['BCS']))
c3.metric('Cmax (mg/L)', f"{pk_e['Cmax']:.3f}", delta=f"Space {pk_s['Cmax']:.3f}", delta_color='off')
c4.metric('AUC (mg·h/L)', f"{pk_e['AUC']:.2f}", delta=f"Space {pk_s['AUC']:.2f}", delta_color='off')
if dose_rec:
    c5.metric('Space dose', f"{dose_rec['space_dose_mg']:.1f} mg")
    c6.metric('Change', f"{dose_rec['adjustment_pct']:+.0f}%")
else:
    c5.metric('Space dose', '—')
    c6.metric('Change', '—')

tab_sim, tab_dose, tab_props, tab_catalog, tab_custom = st.tabs([
    'Simulation', 'Dose plan', 'Properties', f'All {n_drugs} drugs', 'Custom SMILES',
])

with tab_sim:
    st.plotly_chart(plot_pk_curves(result, drug, mission_days), width='stretch')

    left, right = st.columns([1, 1])
    with left:
        st.markdown('**PK parameter comparison**')
        st.dataframe(
            pd.DataFrame({
                'Parameter': ['Cmax (mg/L)', 'Tmax (h)', 'AUC (mg·h/L)', 'Half-life (h)'],
                'Earth': [pk_e['Cmax'], pk_e['Tmax'], pk_e['AUC'], pk_e['t12']],
                'Space': [pk_s['Cmax'], pk_s['Tmax'], pk_s['AUC'], pk_s['t12']],
                'Change': [
                    pk_change_pct(pk_e[k], pk_s[k])
                    for k in ['Cmax', 'Tmax', 'AUC', 't12']
                ],
            }),
            hide_index=True,
            width='stretch',
        )
    with right:
        st.markdown('**Mission-day dose factor**')
        st.caption('Shaded bands: acute (0–3d) · adaptation (3–14d) · chronic (14–180d)')
        st.plotly_chart(plot_dose_timeline(drug, body_weight, dose_mg), width='stretch')

with tab_dose:
    if dose_rec:
        d1, d2, d3, d4 = st.columns(4)
        d1.metric('Earth dose', f"{dose_rec['earth_dose_mg']:.0f} mg")
        d2.metric('Recommended space dose', f"{dose_rec['space_dose_mg']:.0f} mg")
        d3.metric('AUC-matched dose', f"{dose_rec.get('space_dose_auc_mg', 0):.0f} mg")
        d4.metric('Adjustment', f"{dose_rec['adjustment_pct']:+.1f}%")

        factor = dose_rec['adjustment_factor']
        if factor > 1.15:
            st.warning('Increase above 15% — monitor peak exposure and toxicity limits.')
        elif factor < 0.85:
            st.warning('Reduction suggested — confirm efficacy at trough concentrations.')
        else:
            st.success('Within ±15% — Earth dosing may be acceptable with monitoring.')

        st.info(
            f"**Phase:** {profile['space_phase']} (mission day {mission_days}). "
            f"**{dose_rec['rationale']}** "
            f"Refs: Gandia 2003, Kovachevich 2009, Polyakov 2021."
        )
    else:
        st.warning('Could not compute dose recommendation for this profile.')

    st.markdown('**Space physiology modifiers**')
    mods = profile['space_modifiers']
    m1, m2, m3, m4 = st.columns(4)
    m1.metric('Absorption (ka)', f"× {mods['ka_factor']}")
    m2.metric('Bioavailability (F)', f"× {mods['F_factor']}")
    m3.metric('Volume (Vd)', f"× {mods['Vd_factor']}")
    m4.metric('Elimination (ke)', f"× {mods['ke_factor']}")

with tab_props:
    r1c1, r1c2, r1c3, r1c4 = st.columns(4)
    r1c1.metric('Molecular weight', f"{profile['MW']:.1f} g/mol")
    r1c2.metric('LogP', f"{profile['logP']:.2f}")
    r1c3.metric('PSA', f"{profile['PSA']:.1f} Å²")
    r1c4.metric('QED score', f"{profile['QED']:.3f}")

    r2c1, r2c2, r2c3, r2c4 = st.columns(4)
    r2c1.metric('H-bond donors / acceptors', f"{profile['HBD']} / {profile['HBA']}")
    r2c2.metric('Protein binding', f"{profile['protein_binding'] * 100:.0f}%")
    r2c3.metric('Vd (Earth)', f"{profile['Vd_earth']:.1f} L")
    r2c4.metric('Vd (Space)', f"{profile['Vd_space']:.1f} L")

    r3c1, r3c2 = st.columns(2)
    r3c1.metric('Route', profile.get('route', route))
    r3c2.metric('Space phase', profile['space_phase'])

    st.markdown('**SMILES**')
    st.code(profile.get('smiles', ''), language=None)
    st.caption(f"CYP: {meta.get('cyp', '—')} · Transporter: {meta.get('transporter', '—')}")

with tab_catalog:
    st.markdown(f'**{n_drugs} compounds** — ISS kit and space-relevant medications')
    df_cat = pd.DataFrame(get_catalog_table())
    st.dataframe(
        df_cat,
        hide_index=True,
        width='stretch',
        height=500,
        column_config={
            'Dose (mg)': st.column_config.NumberColumn(format='%.2f'),
            't½ (h)': st.column_config.NumberColumn(format='%.1f'),
        },
    )
    st.download_button(
        'Download catalog (CSV)',
        df_cat.to_csv(index=False),
        file_name='spacepk_drug_catalog.csv',
        mime='text/csv',
    )

with tab_custom:
    st.markdown('Analyze any molecule by entering a **SMILES** string.')
    cname = st.text_input('Compound name', value='Custom compound', key='custom_name')
    smiles_in = st.text_input('SMILES', placeholder='CC(=O)Nc1ccc(O)cc1', key='custom_smiles')
    cdose = st.number_input('Dose (mg)', min_value=0.1, max_value=5000.0, value=500.0, key='custom_dose_val')

    if st.button('Run analysis', type='primary', key='custom_run'):
        if not smiles_in.strip():
            st.error('Enter a valid SMILES string.')
        else:
            try:
                cp = analyze_drug(
                    cname, smiles=smiles_in.strip(), dose_mg=cdose,
                    mission_days=mission_days, body_weight_kg=body_weight, verbose=False,
                )
                cr = run_pbpk_analysis(cp, mission_days=mission_days, body_weight=body_weight, verbose=False)
                st.plotly_chart(plot_pk_curves(cr, cname, mission_days), width='stretch')
                if cr.get('dose_recommendation'):
                    dr = cr['dose_recommendation']
                    st.success(f"Recommended space dose: **{dr['space_dose_mg']:.0f} mg** (×{dr['adjustment_factor']:.2f})")
                else:
                    st.info('No oral dose recommendation for this route/profile.')
            except Exception as exc:
                st.error(f'Analysis failed: {exc}')

st.divider()
st.caption(
    'SpacePK research framework · '
    '[GitHub](https://github.com/prakash023-hub/spacepk) · '
    'For publication use the public GitHub or Streamlit Cloud URL — not localhost.'
)

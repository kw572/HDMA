import os
import runpy
from pathlib import Path

import streamlit as st

APP_DIR = Path(__file__).resolve().parent
APPS = {
    "FASTA Editor": APP_DIR / "fasta.py",
    "Genome Browser": APP_DIR / "genome.browser.py",
}

st.set_page_config(
    page_title="ChromBPNet Tools",
    page_icon="DNA",
    layout="wide",
)

with st.sidebar:
    st.markdown("### Tools")
    selected_app = st.radio(
        "Choose a tool",
        list(APPS.keys()),
        label_visibility="collapsed",
    )

st.title(selected_app)
os.environ["CHROMBPNET_EMBEDDED_APP"] = "1"
try:
    runpy.run_path(str(APPS[selected_app]), run_name="__main__")
finally:
    os.environ.pop("CHROMBPNET_EMBEDDED_APP", None)

# FASTA_processing/streamlit_app.py
from __future__ import annotations
import io, zipfile
from collections import Counter
from typing import Iterable, Tuple

import pandas as pd
import numpy as np
import streamlit as st
from Bio import SeqIO

st.set_page_config(page_title="FASTA Processing – Stats", page_icon="🧬", layout="wide")
st.title("🧬 FASTA Processing – Statystyki z plików FASTA")

st.sidebar.header("Wejście")
uploaded = st.sidebar.file_uploader("Wgraj FASTA lub ZIP z FASTA", type=["fa", "fasta", "fna", "faa", "fas", "fa.gz", "gz", "zip"])

min_len = st.sidebar.number_input("Minimalna długość", min_value=0, value=0, step=100)
max_len = st.sidebar.number_input("Maksymalna długość (0 = bez limitu)", min_value=0, value=0, step=1000)
topn_kmers = st.sidebar.number_input("Top k-mer (k=3) – ile pozycji pokazać", min_value=0, value=10, step=1)

def read_fasta_stream(file_name: str, data: bytes) -> Iterable[SeqIO.SeqRecord]:
    # .gz wprost (Streamlit rozpakowuje gz jako raw bytes – BioPython ogarnie)
    if file_name.endswith(".gz"):
        return SeqIO.parse(io.TextIOWrapper(io.BytesIO(data), encoding="utf-8", newline=""), "fasta")

    # ZIP – wyciągamy wszystkie pliki *.fa*, scalając rekordy
    if file_name.endswith(".zip"):
        zf = zipfile.ZipFile(io.BytesIO(data))
        for name in zf.namelist():
            if any(name.lower().endswith(ext) for ext in (".fa", ".fasta", ".fna", ".faa", ".fas")):
                with zf.open(name) as fh:
                    yield from SeqIO.parse(io.TextIOWrapper(fh, encoding="utf-8", newline=""), "fasta")
        return

    # Zwykły FASTA
    return SeqIO.parse(io.StringIO(data.decode("utf-8", errors="ignore")), "fasta")

def basic_stats(seq: str) -> Tuple[int, float, float]:
    seq_upper = seq.upper()
    length = len(seq_upper)
    if length == 0:
        return 0, 0.0, 0.0
    gc = (seq_upper.count("G") + seq_upper.count("C")) / length * 100.0
    n_pct = seq_upper.count("N") / length * 100.0
    return length, gc, n_pct

def kmers(seq: str, k: int = 3) -> Counter:
    s = seq.upper()
    return Counter(s[i:i+k] for i in range(0, max(len(s)-k+1, 0)))

if not uploaded:
    st.info("Wgraj plik FASTA (lub ZIP) po lewej, aby rozpocząć.")
    st.stop()

data = uploaded.read()
records = list(read_fasta_stream(uploaded.name, data))
if isinstance(records, Iterable) and not isinstance(records, list):
    records = list(records)  # w razie generatora

if len(records) == 0:
    st.error("Nie znaleziono żadnych rekordów FASTA.")
    st.stop()

# Tabela sekwencji
rows = []
kmer_counter = Counter()
for r in records:
    seq = str(r.seq)
    length, gc, n_pct = basic_stats(seq)
    rows.append({
        "id": r.id,
        "description": r.description,
        "length": length,
        "gc_percent": gc,
        "n_percent": n_pct,
    })
    if topn_kmers > 0:
        kmer_counter.update(kmers(seq, k=3))

df = pd.DataFrame(rows)

# Filtrowanie po długości
if min_len > 0:
    df = df[df["length"] >= min_len]
if max_len > 0:
    df = df[df["length"] <= max_len]

c1, c2, c3, c4 = st.columns(4)
c1.metric("Liczba sekwencji", f"{len(df):,}".replace(",", " "))
c2.metric("Średnia długość", f"{df['length'].mean():.1f}" if len(df) else "–")
c3.metric("Mediana długości", f"{df['length'].median():.1f}" if len(df) else "–")
c4.metric("Średni GC%", f"{df['gc_percent'].mean():.2f}" if len(df) else "–")

st.subheader("Tabela sekwencji")
st.dataframe(df, use_container_width=True)

st.subheader("Histogram długości")
if len(df):
    hist, edges = np.histogram(df["length"], bins=50)
    st.bar_chart(hist)

st.subheader("Top k-mery (k=3)")
if topn_kmers > 0:
    km_df = pd.DataFrame(kmer_counter.most_common(topn_kmers), columns=["kmer", "count"])
    st.table(km_df)

# Eksport
st.subheader("Eksport")
st.download_button("Pobierz CSV", df.to_csv(index=False).encode("utf-8"), file_name="fasta_stats.csv", mime="text/csv")

st.caption("Made with ❤️ Streamlit + Biopython")

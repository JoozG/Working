import streamlit as st
import pandas as pd
from fasta_processing import read_fasta_bytes, dnaOperations, rnaConverter, proteinOperations

st.set_page_config(page_title="FASTA Processing", page_icon="🧬", layout="wide")
st.title("🧬 FASTA Processing – demo Twojej biblioteki")

upl = st.file_uploader("Wgraj plik FASTA", type=["fa","fasta","fna","faa","fas"])
mode = st.selectbox("Tryb analizy", ["DNA", "RNA→Protein", "Protein"])

if not upl:
    st.info("Wgraj plik, aby rozpocząć.")
    st.stop()

seqs = read_fasta_bytes(upl.read())  # <= Twoja funkcja z io_utils
rows = []

if mode == "DNA":
    for sid, seq in seqs.items():
        dna = dnaOperations(seq)
        rows.append({
            "id": sid,
            "length": len(seq),
            "GC%": dna.gc_content(),
            "A": dna.nuc_count()["A"],
            "C": dna.nuc_count()["C"],
            "G": dna.nuc_count()["G"],
            "T": dna.nuc_count()["T"],
            "revcomp_5to3": "".join(dna.reverse_complement()),
        })
elif mode == "RNA→Protein":
    for sid, rna in seqs.items():
        prot = rnaConverter(rna).rna_to_protein(one_letter=True)
        rows.append({
            "id": sid,
            "aa_len": len(prot),
            "protein(1-letter)": "".join(prot[:100]) + ("..." if len(prot) > 100 else "")
        })
else:  # Protein
    for sid, aa in seqs.items():
        p = proteinOperations(aa)
        rows.append({
            "id": sid,
            "length": len(aa),
            "MW (Da)": p.molecular_weight(),
            "pI~": p.isoelectric_point(),
            "Hydroph.": p.hydrophobicity_score(),
        })

st.dataframe(pd.DataFrame(rows), use_container_width=True)

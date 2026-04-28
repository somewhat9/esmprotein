"""
variant summary file
Size: 436,222,584
Released: 2024-03-31 19:30:16
Last modified: 2026-04-21 11:22:06

Pipeline:
  1. Load variant_summary.txt
  2. Filter: GRCh38 + single SNV + missense (HGVS p.) + 2-star+ review + benign/pathogenic
  3. Map ClinVar GeneID -> UniProt accession via HUMAN_9606_idmapping.dat
  4. Verify wt residue against the UniProt FASTA sequence; build mutated sequence
  5. Save outputs/clinvar_subs.csv
"""

import re
from pathlib import Path

import pandas as pd

_SCRIPT_DIR = Path(__file__).parent
_DEFAULT_CLINVAR_DIR = _SCRIPT_DIR.parent / "data" / "clinvar"
_DEFAULT_UNIPROT_DIR = _SCRIPT_DIR.parent / "data" / "UniProt"

_DTYPES = {
    "#AlleleID": "int32",
    "Type": "category",
    "Name": "string",
    "GeneID": "int32",
    "GeneSymbol": "category",
    "HGNC_ID": "string",
    "ClinicalSignificance": "category",
    "ClinSigSimple": "int8",
    "LastEvaluated": "string",
    "RS# (dbSNP)": "int64",
    "nsv/esv (dbVar)": "string",
    "RCVaccession": "string",
    "PhenotypeIDS": "string",
    "PhenotypeList": "string",
    "Origin": "category",
    "OriginSimple": "category",
    "Assembly": "category",
    "ChromosomeAccession": "string",
    "Chromosome": "category",
    "Start": "int32",
    "Stop": "int32",
    "ReferenceAllele": "string",
    "AlternateAllele": "string",
    "Cytogenetic": "string",
    "ReviewStatus": "category",
    "NumberSubmitters": "int16",
    "Guidelines": "string",
    "TestedInGTR": "category",
    "OtherIDs": "string",
    "SubmitterCategories": "int8",
    "VariationID": "int32",
    "PositionVCF": "int32",
    "ReferenceAlleleVCF": "string",
    "AlternateAlleleVCF": "string",
    "SomaticClinicalImpact": "string",
    "SomaticClinicalImpactLastEvaluated": "string",
    "ReviewStatusClinicalImpact": "string",
    "Oncogenicity": "string",
    "OncogenicityLastEvaluated": "string",
    "ReviewStatusOncogenicity": "string",
    "SCVsForAggregateGermlineClassification": "string",
    "SCVsForAggregateSomaticClinicalImpact": "string",
    "SCVsForAggregateOncogenicityClassification": "string",
}

_MISSENSE_RE = r'p\.[A-Z][a-z]{2}\d+(?!Ter|fs|del|dup|ins|ext|=)[A-Z][a-z]{2}'

# "NM_007294.4(BRCA1):c.5096G>A (p.Arg1699Gln)" -> ("Arg", "1699", "Gln")
_PROT_RE = re.compile(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)')

_TWO_STAR_PLUS = {
    "reviewed by expert panel",          # 3 stars
    "practice guideline",                # 4 stars
    "criteria provided, multiple submitters, no conflicts",  # 2 stars
}

_BENIGN = ['Benign', 'Likely benign', 'Benign/Likely benign']
_PATHOGENIC = ['Pathogenic', 'Likely pathogenic', 'Pathogenic/Likely pathogenic']

_AA_3TO1 = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
}


def _load_clinvar_filtered(clinvar_dir):
    df = pd.read_csv(
        Path(clinvar_dir) / 'variant_summary.txt',
        sep="\t",
        dtype=_DTYPES,
        na_values=["-", "na", ""],
    )
    df = df[df['Assembly'] == 'GRCh38']
    df = df[df['Type'] == 'single nucleotide variant']
    df = df[df['Name'].str.contains(_MISSENSE_RE, regex=True, na=False)]
    df = df[df['ReviewStatus'].isin(_TWO_STAR_PLUS)]
    df = df[df['ClinicalSignificance'].isin(_BENIGN + _PATHOGENIC)]
    df = df.copy()
    df['label'] = df['ClinicalSignificance'].isin(_PATHOGENIC).astype(int)
    return df


def _load_geneid_to_uniprot(idmapping_path):
    """GeneID (str of int) -> list of UniProt accessions (canonical first)."""
    mapping = {}
    with open(idmapping_path) as f:
        for line in f:
            ac, id_type, value = line.rstrip("\n").split("\t")
            if id_type != "GeneID":
                continue
            # Prefer canonical accession (no '-N' isoform suffix) at the front.
            bucket = mapping.setdefault(value, [])
            if "-" in ac:
                bucket.append(ac)
            else:
                bucket.insert(0, ac)
    return mapping


def _load_fasta(fasta_path):
    """UniProt accession -> sequence string (sp| entries only)."""
    seqs = {}
    ac = None
    chunks = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                if ac is not None:
                    seqs[ac] = "".join(chunks)
                header = line[1:].split(None, 1)[0]  # e.g. "sp|P31946|1433B_HUMAN"
                parts = header.split("|")
                ac = parts[1] if len(parts) >= 2 and parts[0] == "sp" else None
                chunks = []
            elif ac is not None:
                chunks.append(line.strip())
        if ac is not None:
            seqs[ac] = "".join(chunks)
    return seqs


def _parse_protein_change(name):
    m = _PROT_RE.search(name)
    if m is None:
        return None
    wt3, pos_str, mut3 = m.groups()
    wt = _AA_3TO1.get(wt3)
    mut = _AA_3TO1.get(mut3)
    if wt is None or mut is None:
        return None
    return wt, int(pos_str), mut


def load_clinvar(clinvar_dir=None, uniprot_dir=None):
    if clinvar_dir is None:
        clinvar_dir = _DEFAULT_CLINVAR_DIR
    if uniprot_dir is None:
        uniprot_dir = _DEFAULT_UNIPROT_DIR

    df = _load_clinvar_filtered(clinvar_dir)
    geneid_to_acs = _load_geneid_to_uniprot(Path(uniprot_dir) / "HUMAN_9606_idmapping.dat")
    ac_to_seq = _load_fasta(Path(uniprot_dir) / "uniprot_human_reviewed.fasta")

    records = []
    for row in df.itertuples(index=False):
        parsed = _parse_protein_change(row.Name)
        if parsed is None:
            continue
        wt_aa, pos_1based, mut_aa = parsed
        idx = pos_1based - 1

        candidates = geneid_to_acs.get(str(row.GeneID), [])
        chosen_ac = None
        wt_seq = None
        for ac in candidates:
            seq = ac_to_seq.get(ac)
            if seq is None:
                continue
            if 0 <= idx < len(seq) and seq[idx] == wt_aa:
                chosen_ac = ac
                wt_seq = seq
                break
        if chosen_ac is None:
            continue

        mut_seq = wt_seq[:idx] + mut_aa + wt_seq[idx + 1:]
        records.append({
            "protein_id": chosen_ac,
            "gene_symbol": str(row.GeneSymbol),
            "wt_seq": wt_seq,
            "mut_seq": mut_seq,
            "position": idx,
            "wt_aa": wt_aa,
            "mut_aa": mut_aa,
            "label": int(row.label),
        })

    return pd.DataFrame(records).drop_duplicates(
        subset=["protein_id", "position", "wt_aa", "mut_aa"]
    ).reset_index(drop=True)


if __name__ == "__main__":
    out_dir = _SCRIPT_DIR.parent / "outputs"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "clinvar_subs.csv"
    df = load_clinvar()
    df.to_csv(out_path, index=False)

    n = len(df)
    n_path = int(df["label"].sum())
    n_benign = n - n_path
    ratio = max(n_path, n_benign) / max(min(n_path, n_benign), 1)
    summary = (
        f"Total variants: {n}\n"
        f"Benign:         {n_benign} ({100 * n_benign / n:.1f}%)\n"
        f"Pathogenic:     {n_path} ({100 * n_path / n:.1f}%)\n"
        f"Majority:minority ratio: {ratio:.2f}:1\n"
    )
    (out_dir / "clinvar_class_balance.txt").write_text(summary)
    print(f"Saved {n} records to {out_path}\n")
    print(summary)
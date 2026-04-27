import pandas as pd
from pathlib import Path

_SCRIPT_DIR = Path(__file__).parent
_DEFAULT_DATA_DIR = _SCRIPT_DIR.parent / "data" / "proteingym"


def load_proteingym_subs(data_dir=None):
    if data_dir is None:
        data_dir = _DEFAULT_DATA_DIR
    ref = pd.read_csv(Path(data_dir) / "DMS_substitutions_reference.csv")

    all_assays = []
    assay_dir = Path(data_dir) / "DMS_ProteinGym_substitutions"

    for csv_path in assay_dir.glob("*.csv"):
        df = pd.read_csv(csv_path)
        dms_id = csv_path.stem

        #match reference
        ref_row = ref[ref["DMS_id"] == dms_id].iloc[0]
        wt_seq = ref_row["target_seq"]
        protein_id = ref_row["UniProt_ID"]

        for _, row in df.iterrows():
            mut = row["mutant"]

            if ":" in mut:  # skip if multiple mutations
                continue
            wt_aa = mut[0]
            mut_aa = mut[-1]
            position = int(mut[1:-1]) - 1  # convert to 0-index

            mut_seq = row["mutated_sequence"]

            all_assays.append({
                "protein_id": protein_id,
                "wt_seq": wt_seq,
                "mut_seq": mut_seq,
                "position": position,
                "wt_aa": wt_aa,
                "mut_aa": mut_aa,
                "DMS_score": row["DMS_score"]
            })
    return pd.DataFrame(all_assays)


if __name__ == "__main__":
    out_path = _SCRIPT_DIR.parent / "outputs" / "proteingym_subs.csv"
    df = load_proteingym_subs()
    df.to_csv(out_path, index=False)
    print(f"Saved {len(df)} records to {out_path}")
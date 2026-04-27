"""
variant summary file
Size: 436,222,584
Released: 2024-03-31 19:30:16
Last modified: 2026-04-21 11:22:06
"""

import pandas as pd
from pathlib import Path

_SCRIPT_DIR = Path(__file__).parent
_DEFAULT_DATA_DIR = _SCRIPT_DIR.parent / "data" / "clinvar"

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

df = pd.read_csv(
    Path(_DEFAULT_DATA_DIR) / 'variant_summary.txt',
    sep="\t",
    dtype=_DTYPES,
    na_values=["-", "na", ""],
)

df = df[df['Assembly'] == 'GRCh38']
df = df[df['Type'] == 'single nucleotide variant']

_MISSENSE_RE = r'p\.[A-Z][a-z]{2}\d+(?!Ter|fs|del|dup|ins|ext|=)[A-Z][a-z]{2}'
df = df[df['Name'].str.contains(_MISSENSE_RE, regex=True, na=False)]

TWO_STAR_PLUS = {
    "reviewed by expert panel",  # 3 stars
    "practice guideline",  # 4 stars
    "criteria provided, multiple submitters, no conflicts"  # 2 stars
}

df = df[df['ReviewStatus'].isin(TWO_STAR_PLUS)]

benign = ['Benign', 'Likely benign', 'Benign/Likely benign']
pathogenic = ['Pathogenic', 'Likely pathogenic', 'Pathogenic/Likely pathogenic']

df = df[df['ClinicalSignificance'].isin(benign + pathogenic)]
df['label'] = df['ClinicalSignificance'].isin(pathogenic).astype(int)



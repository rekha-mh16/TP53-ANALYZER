from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
from collections import Counter


# ============================================================
# TP53 COMPARATIVE PROTEIN ANALYSIS
# ============================================================

# FASTA files present in the project folder
fasta_files = [
    "human tp53.fasta",
    "mouse tp53.fasta",
    "rat tp53.fasta",
    "rabbit tp53.fasta",
    "rhesus tp53.fasta",
    "dog tp53.fasta",
    "cattle tp53.fasta"
]


# Standard amino acids
standard_amino_acids = list("ACDEFGHIKLMNPQRSTVWY")


# ============================================================
# FUNCTION TO ANALYZE ONE PROTEIN
# ============================================================

def analyze_protein(filename):

    record = SeqIO.read(filename, "fasta")

    sequence = str(record.seq).upper()

    # Remove spaces and newline characters
    sequence = sequence.replace(" ", "").replace("\n", "")

    analysis = ProteinAnalysis(sequence)

    # --------------------------------------------------------
    # BASIC INFORMATION
    # --------------------------------------------------------

    length = len(sequence)

    molecular_weight = analysis.molecular_weight()

    isoelectric_point = analysis.isoelectric_point()

    aromaticity = analysis.aromaticity()

    instability = analysis.instability_index()

    gravy = analysis.gravy()

    # --------------------------------------------------------
    # AMINO ACID COMPOSITION
    # --------------------------------------------------------

    amino_acid_count = Counter(sequence)

    amino_acid_percentage = {}

    for aa in standard_amino_acids:

        count = amino_acid_count.get(aa, 0)

        percentage = (count / length) * 100

        amino_acid_percentage[aa] = round(percentage, 2)

    # --------------------------------------------------------
    # IMPORTANT TP53 RESIDUES
    # --------------------------------------------------------

    important_residues = {

        "Arg175":
            sequence[174] if length >= 175 else "-",

        "Gly245":
            sequence[244] if length >= 245 else "-",

        "Arg248":
            sequence[247] if length >= 248 else "-",

        "Arg249":
            sequence[248] if length >= 249 else "-",

        "Arg273":
            sequence[272] if length >= 273 else "-",

        "Arg282":
            sequence[281] if length >= 282 else "-"
    }

    # --------------------------------------------------------
    # RESULT
    # --------------------------------------------------------

    result = {

        "File": filename,

        "Sequence_ID": record.id,

        "Length": length,

        "Molecular_Weight": round(
            molecular_weight, 2
        ),

        "pI": round(
            isoelectric_point, 2
        ),

        "Aromaticity": round(
            aromaticity, 4
        ),

        "Instability_Index": round(
            instability, 2
        ),

        "GRAVY": round(
            gravy, 4
        )
    }

    return (
        result,
        amino_acid_percentage,
        sequence,
        important_residues
    )


# ============================================================
# ANALYZE ALL FASTA FILES
# ============================================================

results = []

sequences = {}

amino_acid_data = {}

important_residue_data = {}


for file in fasta_files:

    try:

        (
            result,
            aa_percentage,
            sequence,
            important
        ) = analyze_protein(file)

        results.append(result)

        sequences[file] = sequence

        amino_acid_data[file] = aa_percentage

        important_residue_data[file] = important

        print("\nAnalyzed:", file)

    except FileNotFoundError:

        print("\nFile not found:", file)

    except Exception as e:

        print("\nError analyzing:", file)

        print(e)


# ============================================================
# TP53 SUMMARY TABLE
# ============================================================

df = pd.DataFrame(results)


print("\n")
print("=" * 80)
print("TP53 ANALYSIS SUMMARY")
print("=" * 80)

print(
    df.to_string(index=False)
)


# ============================================================
# SAVE SUMMARY
# ============================================================

df.to_csv(
    "TP53_analysis_summary.csv",
    index=False
)

print("\nSummary saved as:")

print(
    "TP53_analysis_summary.csv"
)


# ============================================================
# AM
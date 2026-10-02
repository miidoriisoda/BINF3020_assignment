# -----------------------------------------------------------------
# FINDING RELATED PROTEINS
# -----------------------------------------------------------------
# ------------------- Set-up --------------------------------------
from helper_functions import global_alignment, blosum62_matrix, \
    percent_identity
from Bio import Entrez, SeqIO
from Bio.SeqFeature import SimpleLocation
import sys; sys.path.append(".")
import json
from Bio import Entrez
Entrez.email = "z5476602@unsw.edu.au"

accession_codes = {
    # 6 known human coronaviruses
    "Human-SARS": "NC_004718",
    "Human-MERS": "NC_019843",
    "Human-HCoV-OC43": "NC_006213",
    "Human-HCoV-229E": "NC_002645",
    "Human-HCoV-NL63": "NC_005831",
    "Human-HCoV-HKU1": "NC_006577",

    # Bat
    "Bat-CoV MOP1": "EU420138",
    "Bat-CoV HKU8": "NC_010438",
    "Bat-CoV HKU2": "NC_009988",
    "Bat-CoV HKU5": "NC_009020",
    "Bat-CoV RaTG13": "MN996532",
    "Bat-CoV-ENT": "NC_003045",

    # Other animals
    "Hedgehog-CoV 2012-174/GER/2012": "NC_039207",
    "Pangolin-CoV MP789": "MT121216",
    "Rabbit-CoV HKU14": "NC_017083",
    "Duck-CoV isolate DK/GD/27/2014": "NC_048214",
    "Feline infectious peritonitis virus": "NC_002306",
    "Giraffe-CoV US/OH3/2003": "EF424623",
    "Murine-CoV MHV/BHKR_lab/USA/icA59_L94P/2012": "KF268338",
    "Equine-CoV Obihiro12-2": "LC061274",
}

# ------------------- Extract spike proteins ----------------------
# get sars-cov-2 spike from provided fasta file
sequence = SeqIO.read("data/sars_cov_2.fa", "fasta")
location = SimpleLocation(21562, 25384, 1)
sars_cov_2 = location.extract(sequence.seq).translate(cds=True)

# extract spike protein sequence from candidates in accession_codes
spike_sequence = {}
for name, id in accession_codes.items():
        stream = Entrez.efetch(db="nuccore", id=f"{id}", rettype="gb",\
                        retmode="text")
        record = SeqIO.read(stream, "gb")
        for feature in record.features:
            if feature.type == "CDS" and (
                feature.qualifiers.get("gene", [None])[0] == "S"
                or feature.qualifiers.get("product", [None])[0] == "spike protein"
            ):
                nucleotide_sequence = feature.location.extract(record.seq)
                translation = nucleotide_sequence.translate(cds=True)
                spike_sequence[name] = translation

# ------------------- Globally align extracted sequences ----------
# store results (alignment & % identity) in all_results
all_results = []

for name, seq in spike_sequence.items():
    alignment = global_alignment(sars_cov_2, seq, blosum62_matrix)

    all_results.append({
        "name": name,
        "score": alignment[2],
        "seq1": alignment[0],
        "seq2": alignment[1],
        "ident": percent_identity(alignment[0], alignment[1])
    })

# sort results by descending alignment score
all_results = sorted(
    all_results,
    key=lambda x: x["score"],
    reverse=True
)

# print top 3 results
for result in all_results[:3]:
    print(
        result["name"],
        "score:", result["score"],
        "identity:", result["ident"]
    )

# save all alignments to JSON
with open("all_alignments.json", "w") as file:
    json.dump(all_results, file, indent=4)

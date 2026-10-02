# ------------------- Set-up --------------------------------------
import json
from Bio import Entrez, SeqIO
from Bio.SeqFeature import SimpleLocation
from helper_functions import local_alignment, blosum90_matrix, alignment_stats
Entrez.email = "z5476602@unsw.edu.au"

og_sequence = SeqIO.read("data/sars_cov_2.fa", "fasta")

orf_candidates = {
    "ORF-1": (1, 11995, 13483),
    "ORF-2": (1, 26792, 27191),
    "ORF-3": (1, 23650, 25384),
    "ORF-4": (-1, 421, 667),
    "ORF-5": (1, 9133, 13483),
}

top_3 = {
    "Human-SARS": "NC_004718",
    "Bat-CoV RaTG13": "MN996532",
    "Pangolin-CoV MP789": "MT121216",
}

ref_database = []

# ------------------- Extract CDS for top 3 genomes ---------------
# Extract coding regions for each genome
ref_database = []
for name, acc in top_3.items():
    with Entrez.efetch(db="nuccore", id=acc, rettype="gb", retmode="text") as stream:
        record = SeqIO.read(stream, "gb")
    for feature in record.features:
        if feature.type != "CDS":
            continue
        try:
            prot = str(feature.location.extract(record.seq).translate(cds=True))
        except Exception:
            # fall back just in case
            prot = feature.qualifiers["translation"][0]
        ref_database.append({
            "virus_name": name,
            "gene_name": feature.qualifiers.get("gene", [None])[0],
            "function": feature.qualifiers.get("product", [None])[0],
            "prot_sequence": prot,
        })

# ------------------- Compute local alignment ---------------------
# Find protein sequence for ORF candidates
orf_database = []
for orf, pos in orf_candidates.items():
    location = SimpleLocation(pos[1], pos[2], pos[0])
    cds = location.extract(og_sequence.seq).translate(cds=True)
    orf_database.append({
        "name": orf,
        "position": pos,
        "sequence": cds
    })

# Local alignment with BLOSUM90
# As looking for 'perform same function' use high fidelity BLOSUM
# matrix (i.e. BLOSUM90)
results = []
for orf in orf_database:
    for ref in ref_database:
        a1, a2, score = local_alignment(ref["prot_sequence"], orf["sequence"],
                                        blosum90_matrix)
        ident, cov = alignment_stats(a1, a2, len(orf["sequence"]))
        results.append({
            "orf": orf["name"], "virus": ref["virus_name"],
            "gene": ref["gene_name"],
            "function": ref["function"],
            "score": score, "identity": ident, "coverage": cov,
            "seq1": a1, "seq2": a2,
        })

best = {}
for r in results:
    if r["orf"] not in best or r["score"] > best[r["orf"]]["score"]:
        best[r["orf"]] = r

for orf, r in best.items():
    print(f'{orf}: {r["virus"]} | {r["gene"]} | score={r["score"]:.0f} '
          f'identity={r["identity"]:.2f} coverage={r["coverage"]:.2f}')

# Write it to a file in case vs-code doesn't load again QAQ
with open("best_local_alignments.json", "w") as file:
    json.dump(best, file, indent=4)

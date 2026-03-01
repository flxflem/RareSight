[README.md](https://github.com/user-attachments/files/25665271/README.md)
# RareSight — Gene Therapy Target Analysis Tool

A local web application for physicians and researchers to analyze therapeutic targets
for monogenic diseases. Queries OMIM, UniProt, ClinVar, Ensembl, gnomAD, Open Targets,
and PubMed to produce a structured clinical report.

## Quick Start

```bash
# 1. Install dependencies
pip install flask requests reportlab

# 2. Start the server
python app.py

# 3. Open in browser
http://localhost:5000
```

Or simply run:
```bash
bash start.sh
```

## Files

```
RareSight/
├── app.py              ← Flask backend (API endpoints + PDF generation)
├── pipeline.py         ← Analysis pipeline (all database queries)
├── templates/
│   └── index.html      ← Frontend (single-file HTML/CSS/JS)
├── start.sh            ← Startup script
└── README.md
```

## Usage

1. Enter a **gene symbol** (e.g. `SCN1A`, `MECP2`, `BRAF`)
2. Optionally enter a **patient variant** in HGVS cDNA format (e.g. `c.4849C>T`)
   - Activates variant-level ClinVar lookup
   - Intersects variant position against UniProt protein features
   - Queries gnomAD for population frequency
   - Designs CRISPR guide RNAs from Ensembl genomic sequence
   - Scores literature evidence for gene + variant
3. Click **Run Analysis** (or press Enter)
4. Review the report on the right panel
5. Click **Export PDF Report** to download a formatted clinical report

## Output Sections

| Section | Source |
|---------|--------|
| Gene Identity & Sequence | OMIM, MyGene.info |
| Disease Mechanism | OMIM text mining |
| Literature Evidence | PubMed / NCBI Entrez |
| Patient Variant Analysis | ClinVar, gnomAD, UniProt |
| Therapy Recommendation | Rule-based engine |
| CRISPR Guide RNA Design | Ensembl REST API |
| Protein Features | UniProt |
| Key ClinVar Variants | ClinVar |
| Pathways & Downstream Targets | KEGG, Reactome via MyGene |
| Known Drugs | Open Targets GraphQL |
| Clinical Notes & Warnings | Synthesized from all sources |
| OMIM Summary | OMIM |

## OMIM API Key

The OMIM API key is hardcoded in `pipeline.py` at the top of the file:
```python
OMIM_API_KEY = "your_key_here"
```
Replace with your institution's key if needed.

## Notes

- All data stays local — no patient data is sent to external AI services
- Network calls go only to public biomedical databases (OMIM, NCBI, UniProt, etc.)
- Not intended as a standalone clinical decision tool — always validate with a
  qualified clinical geneticist or bioinformatician

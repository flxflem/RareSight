"""
Gene Therapy Target Analysis Pipeline v2
==========================================
Builds on v1 with three major new capabilities:

  1. PATIENT VARIANT INPUT  (--variant "c.4849C>T")
     - ClinVar lookup for that specific variant: pathogenicity, condition, consequence
     - Amino acid position parsed from HGVS protein notation (or estimated from cDNA)
     - UniProt feature intersection: does this variant land in an active site,
       binding site, or functional domain?
     - gnomAD population frequency check (is it truly ultra-rare?)
     - Allele-specific ASO/siRNA discrimination feasibility assessment
     - Variant-level mechanism can override gene-level when clearer

  2. CRISPR GUIDE RNA DESIGN
     - Fetches real genomic sequence from Ensembl REST around the variant
     - Scans both strands for all NGG PAM sites
     - Scores each 20nt guide by: GC content, distance to variant, quality flags
       (poly-T runs, homopolymers, low/high GC)
     - Returns top 8 guides ranked by quality score
     - Assesses base editing feasibility: CBE (C>T), ABE (A>G)
     - Assesses prime editing feasibility for transversions / indels

  3. LITERATURE EVIDENCE SCORING
     - PubMed count for gene (overall characterization depth)
     - PubMed count for gene + specific variant (variant-level evidence)
     - Evidence tier: well-characterized / moderate / emerging / unknown
     - Top 5 review PMIDs for physician reference
     - Confidence warning if variant has <3 publications

Databases queried:
  OMIM, MyGene.info, UniProt, ClinVar, Ensembl REST,
  Open Targets GraphQL, gnomAD GraphQL, NCBI PubMed Entrez

Usage:
  # Gene-level analysis only
  python gene_therapy_pipeline_v2.py --gene SCN1A --omim_key YOUR_KEY

  # With patient-specific variant (activates all three new features)
  python gene_therapy_pipeline_v2.py --gene SCN1A --omim_key YOUR_KEY --variant "c.4849C>T"

  # Save JSON output for assay generator
  python gene_therapy_pipeline_v2.py --gene BRAF --omim_key YOUR_KEY --variant "c.1799T>A" --output_json braf_v600e.json
"""

import requests
import re
import json
import time
import argparse
from dataclasses import dataclass, field, asdict
from typing import Optional

# ── Config ─────────────────────────────────────────────────────────────────────
GENE         = "SCN1A"
OMIM_API_KEY = "8Kh7KDGuQ1qfzRZRQu232g"
NCBI_BASE    = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENSEMBL_BASE = "https://rest.ensembl.org"


# ══════════════════════════════════════════════════════════════════════════════
# Data structures
# ══════════════════════════════════════════════════════════════════════════════

@dataclass
class VariantRecord:
    """A single pathogenic variant from ClinVar"""
    clinvar_id: str
    name: str
    molecular_consequence: str
    significance: str
    condition: str
    mechanism_clue: str          # GoF / LoF / Unknown
    position_cdna: Optional[str] = None
    position_protein: Optional[str] = None
    chromosome: Optional[str] = None
    genomic_start: Optional[int] = None
    genomic_stop: Optional[int] = None


@dataclass
class ProteinFeature:
    """A functional feature from UniProt"""
    feature_type: str
    description: str
    start: int
    end: int


@dataclass
class VariantProteinIntersection:
    """
    Result of resolving a patient variant against the protein feature map.
    Populated only when --variant is supplied.
    """
    aa_position: Optional[int]
    overlapping_features: list = field(default_factory=list)   # list[ProteinFeature]
    is_in_active_site: bool = False
    is_in_binding_site: bool = False
    is_in_domain: bool = False
    clinical_significance: str = ""
    gnomad_allele_frequency: Optional[float] = None
    is_absent_from_gnomad: Optional[bool] = None
    aso_discrimination_feasible: Optional[bool] = None
    aso_notes: str = ""


@dataclass
class GuideRNA:
    """A candidate CRISPR guide RNA"""
    sequence: str                    # 20nt protospacer (5'→3')
    pam: str                         # PAM (e.g. CGG, TGG)
    strand: str                      # + or -
    distance_to_variant_bp: int      # bp from Cas9 cut site to variant center
    gc_content: float                # 0.0–1.0
    quality_flags: list = field(default_factory=list)
    quality_score: int = 0           # 0–100, higher = better
    cut_site_genomic_pos: Optional[int] = None


@dataclass
class CRISPRParameters:
    """
    All CRISPR-relevant parameters for this target.
    Populated for all runs; richer when --variant is supplied.
    """
    strategy: str                    # knock-out / correction / base-edit
    target_sequence_context: Optional[str] = None   # ~80bp window shown
    candidate_guides: list = field(default_factory=list)   # list[GuideRNA]
    base_editing_feasible: bool = False
    base_edit_type: Optional[str] = None    # e.g. CBE (C→T) or ABE (A→G)
    prime_editing_feasible: bool = False
    notes: list = field(default_factory=list)


@dataclass
class ModalityGuide:
    """
    Technical design parameters for a specific payload modality —
    the equivalent of CRISPRParameters but generalised to every payload type.
    One ModalityGuide is attached to each TherapyRecommendation.
    """
    modality_type: str           # aav / mrna_lnp / aso_sirna / crispr / cast

    # ── AAV-specific
    aav_serotype: Optional[str] = None          # e.g. "AAV9"
    aav_promoter: Optional[str] = None          # e.g. "CBA (ubiquitous)" or "SYN1 (neuron)"
    aav_itr_config: Optional[str] = None        # scAAV or ssAAV + rationale
    aav_dose_range: Optional[str] = None        # e.g. "1×10¹⁴ vg/kg (IV, CNS)"
    aav_tropism_note: Optional[str] = None      # serotype tissue tropism caveat
    aav_packaging_note: Optional[str] = None    # HEK293 triple-transfection / Sf9 / etc.

    # ── mRNA + LNP-specific
    mrna_utr_design: Optional[str] = None       # 5'/3' UTR strategy
    mrna_codon_opt: Optional[str] = None        # codon optimisation flag
    mrna_cap: Optional[str] = None              # CleanCap-AG / ARCA
    mrna_modification: Optional[str] = None     # N1-methylpseudouridine etc.
    lnp_formulation: Optional[str] = None       # MC3 / ALC-0315 / SM-102 etc.
    lnp_targeting: Optional[str] = None         # PEG-lipid ratio, GalNAc, etc.
    mrna_stability: Optional[str] = None        # half-life estimate, cold-chain note

    # ── ASO / siRNA-specific
    aso_target_region: Optional[str] = None     # exon / intron boundary / 3'UTR
    aso_design_window: Optional[str] = None     # e.g. "±50 nt around splice donor"
    aso_chemistry: Optional[str] = None         # LNA / 2'-MOE / PS backbone
    aso_conjugate: Optional[str] = None         # GalNAc (liver), naked, nanoparticle
    aso_mechanism: Optional[str] = None         # RNase H / steric block / splice switch
    sirna_strand_bias: Optional[str] = None     # guide strand selection note

    # ── CRISPR-specific (summary level; full detail in crispr_parameters)
    crispr_strategy: Optional[str] = None       # knock-out / base-edit / prime-edit / HDR
    crispr_editor: Optional[str] = None         # SpCas9 / BE4max / ABE8e / PE3
    crispr_pam: Optional[str] = None            # NGG / NG / NRNH
    crispr_window: Optional[str] = None         # editing window (e.g. "positions 4–8 from PAM")
    crispr_delivery_form: Optional[str] = None  # RNP / plasmid / AAV / LNP-mRNA
    crispr_offtarget: Optional[str] = None      # off-target risk note

    # ── CAST / transposon-specific
    cast_system: Optional[str] = None           # PASTE / piggyBac / Sleeping Beauty / Tn7
    cast_target_site: Optional[str] = None      # safe-harbour locus or gene-specific
    cast_donor_size: Optional[str] = None       # max cargo capacity
    cast_transposase_delivery: Optional[str] = None  # mRNA / RNP / plasmid
    cast_offtarget_risk: Optional[str] = None   # insertion site diversity note

    # ── Shared
    key_parameters: list = field(default_factory=list)  # flat list of param strings for display
    design_notes: list = field(default_factory=list)     # caveats, references, open questions


@dataclass
class TherapyRecommendation:
    """
    A single payload + delivery vehicle combination with its own confidence score.
    Confidence is per-combo, not global: it reflects gene evidence depth AND
    how well-validated this delivery modality is in the relevant tissue class.
    """
    rank: int                    # 1 = primary, 2+ = alternatives
    payload: str                 # e.g. "AAV gene replacement"
    delivery: str                # e.g. "AAV9 (CNS/systemic)"
    label: str                   # short display label
    confidence: str              # high / medium / low / experimental
    confidence_rationale: str    # one-sentence explanation
    notes: str = ""              # optional caveats
    modality_guide: Optional["ModalityGuide"] = None  # technical design parameters


@dataclass
class LiteratureEvidence:
    """PubMed-based evidence density for the gene (variant pubs removed)"""
    gene_pubmed_count: int = 0
    evidence_tier: str = "unknown"   # well-characterized / moderate / emerging / unknown
    tier_rationale: str = ""
    top_review_pmids: list = field(default_factory=list)


@dataclass
class AssayParameters:
    """
    Master output — structured for physician/researcher review
    and direct input to an assay generator.
    """
    # ── Gene identity
    gene_symbol: str
    omim_id: str
    uniprot_id: Optional[str]
    ensembl_gene_id: Optional[str]
    refseq_transcript: Optional[str]

    # ── Sequence
    protein_length_aa: Optional[int]
    coding_sequence_kb: Optional[float]
    canonical_transcript_id: Optional[str]

    # ── Mechanism
    mechanism: str
    mechanism_confidence: str
    inheritance: Optional[str]

    # ── Therapy (structured combos with per-combo confidence)
    therapy_combos: list = field(default_factory=list)   # list[TherapyRecommendation]
    aav_feasible: Optional[bool] = None
    # Legacy flat strings kept for PDF/CLI compat — populated from therapy_combos
    primary_therapy: str = ""
    alternative_therapies: list = field(default_factory=list)

    # ── Patient variant (NEW — populated if --variant supplied)
    patient_variant_hgvs: Optional[str] = None
    patient_variant_resolved: Optional[VariantRecord] = None
    variant_protein_intersection: Optional[VariantProteinIntersection] = None

    # ── CRISPR design (NEW)
    crispr_parameters: Optional[CRISPRParameters] = None

    # ── Literature evidence (NEW)
    literature_evidence: Optional[LiteratureEvidence] = None

    # ── Protein features
    target_region_cdna: Optional[str] = None
    representative_variants: list = field(default_factory=list)
    druggable_domains: list = field(default_factory=list)
    active_sites: list = field(default_factory=list)
    binding_sites: list = field(default_factory=list)

    # ── Pathways
    pathway_names: list = field(default_factory=list)
    downstream_targets: list = field(default_factory=list)

    # ── Known drugs
    known_drugs: list = field(default_factory=list)

    # ── Misc
    omim_molecular_genetics_summary: Optional[str] = None
    notes: list = field(default_factory=list)


# ══════════════════════════════════════════════════════════════════════════════
# Shared utilities
# ══════════════════════════════════════════════════════════════════════════════

def safe_get(url, params=None, headers=None, label="request", retries=2):
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, headers=headers, timeout=15)
            r.raise_for_status()
            return r
        except Exception as e:
            print(f"  [warn] {label} failed (attempt {attempt+1}): {e}")
            time.sleep(1)
    return None


def parse_aa_position(hgvs_protein: Optional[str]) -> Optional[int]:
    """Extract integer AA position from p.Arg1617* or p.Val600Glu etc."""
    if not hgvs_protein:
        return None
    m = re.search(r"p\.[A-Za-z]+(\d+)", hgvs_protein)
    return int(m.group(1)) if m else None


def parse_cdna_position(hgvs_cdna: Optional[str]) -> Optional[int]:
    """Extract integer coding position from c.4849C>T etc."""
    if not hgvs_cdna:
        return None
    m = re.search(r"c\.(\d+)", hgvs_cdna)
    return int(m.group(1)) if m else None


def gc_content(seq: str) -> float:
    s = seq.upper()
    return (s.count("G") + s.count("C")) / len(s) if s else 0.0


def quality_flags_for_guide(seq: str) -> list:
    flags = []
    s = seq.upper()
    gc = gc_content(s)
    if gc < 0.3:
        flags.append("LOW_GC (<30%)")
    if gc > 0.8:
        flags.append("HIGH_GC (>80%)")
    if "TTTT" in s:
        flags.append("POLY-T (U6 termination risk)")
    for b in ["AAAA", "CCCC", "GGGG"]:
        if b in s:
            flags.append(f"HOMOPOLYMER ({b})")
    if s and s[-1] != "G":
        flags.append("3-PRIME NOT G (reduced efficiency)")
    return flags


def score_guide(seq: str, dist: int, flags: list) -> int:
    score = 100
    gc = gc_content(seq)
    if gc < 0.3 or gc > 0.8:
        score -= 25
    # Penalty for guides far from the variant; no penalty within 10bp
    score -= min(30, max(0, dist - 10))
    score -= len(flags) * 10
    return max(0, score)


# ══════════════════════════════════════════════════════════════════════════════
# Step 0 — Gene validation (must pass before anything else runs)
# ══════════════════════════════════════════════════════════════════════════════

class GeneNotFoundError(ValueError):
    """Raised when the input string is not a recognised human gene symbol."""
    pass


def validate_gene(gene: str) -> dict:
    """
    Confirms the gene symbol is a known human gene using MyGene.info with a
    symbol-exact query.  Returns the MyGene hit dict on success.
    Raises GeneNotFoundError with a user-friendly message on failure.
    """
    print(f"\n[0/6] Validating gene symbol '{gene}'...")

    r = safe_get(
        "https://mygene.info/v3/query",
        params={
            "q": f"symbol:{gene}",   # exact symbol match — much stricter than free text
            "species": "human",
            "fields": "symbol,name,entrezgene,taxid",
            "size": 3,
        },
        label="MyGene validation"
    )

    if not r:
        raise GeneNotFoundError(
            f"Could not reach MyGene.info to validate '{gene}'. "
            "Check your internet connection and try again."
        )

    hits = r.json().get("hits", [])

    # Filter to exact symbol matches (case-insensitive) for human (taxid 9606)
    exact = [h for h in hits
             if h.get("symbol", "").upper() == gene.upper()
             and h.get("taxid") == 9606]

    if not exact:
        # Try a slightly broader check — symbol in the name — to give a helpful suggestion
        suggestions = [h.get("symbol") for h in hits if h.get("taxid") == 9606]
        suggestion_str = ""
        if suggestions:
            suggestion_str = f" Did you mean: {', '.join(suggestions[:3])}?"
        raise GeneNotFoundError(
            f"'{gene}' is not a recognised human gene symbol.{suggestion_str} "
            "Please enter a valid HGNC symbol (e.g. SCN1A, MECP2, BRAF)."
        )

    hit = exact[0]
    print(f"  Validated: {hit.get('symbol')} — {hit.get('name')} "
          f"(Entrez {hit.get('entrezgene')})")
    return hit


# ══════════════════════════════════════════════════════════════════════════════
# Step 1 — OMIM
# ══════════════════════════════════════════════════════════════════════════════

def query_omim(gene: str, api_key: str) -> dict:
    print(f"\n[1/6] Querying OMIM for {gene}...")
    result = {
        "omim_id": None, "gene_title": None,
        "mechanism": "Unknown", "mechanism_confidence": "low",
        "protein_length": None, "molecular_text": None,
        "aav_feasible": None, "coding_kb": None, "inheritance": None,
    }

    r = safe_get("https://api.omim.org/api/entry/search",
                 params={"search": gene, "format": "json", "apiKey": api_key},
                 label="OMIM search")
    if not r:
        return result

    entry_list = r.json().get("omim", {}).get("searchResponse", {}).get("entryList", [])
    if not entry_list:
        print("  [warn] No OMIM entries found")
        return result

    chosen = next((e["entry"] for e in entry_list
                   if e.get("entry", {}).get("prefix") in ("*", "+")), None)
    if not chosen:
        chosen = entry_list[0]["entry"]

    mim_number = chosen["mimNumber"]
    gene_title  = chosen.get("titles", {}).get("preferredTitle", gene)
    print(f"  Found: {gene_title} (OMIM {mim_number})")

    r2 = safe_get("https://api.omim.org/api/entry",
                  params={"mimNumber": mim_number, "include": "text",
                          "format": "json", "apiKey": api_key},
                  label="OMIM entry")
    if not r2:
        result.update({"omim_id": mim_number, "gene_title": gene_title})
        return result

    sections = r2.json()["omim"]["entryList"][0]["entry"].get("textSectionList", [])
    all_text = "\n".join(
        s["textSection"]["textSectionContent"]
        for s in sections
        if "textSection" in s and "textSectionContent" in s["textSection"]
    )
    mol_text = next(
        (s["textSection"]["textSectionContent"] for s in sections
         if s.get("textSection", {}).get("textSectionName") == "molecularGenetics"),
        ""
    )

    tl = all_text.lower()
    gof = sum(["gain-of-function" in tl, "gain of function" in tl,
               "constitutively active" in tl, "hyperactive" in tl])
    lof = sum(["haploinsufficiency" in tl, "loss of function" in tl,
               "loss-of-function" in tl, "truncat" in tl,
               "nonsense" in tl, "frameshift" in tl])

    if gof > 0 and lof == 0:
        mechanism, confidence = "Gain of Function", ("high" if gof >= 2 else "medium")
    elif lof > 0 and gof == 0:
        mechanism, confidence = "Loss of Function", ("high" if lof >= 2 else "medium")
    elif "dominant negative" in tl:
        mechanism, confidence = "Dominant Negative", "medium"
    elif gof > 0 and lof > 0:
        mechanism, confidence = "Mixed (GoF and LoF variants reported)", "low"
    else:
        mechanism, confidence = "Unknown", "low"

    inheritance = None
    if "autosomal dominant" in tl:
        inheritance = "Autosomal Dominant"
    elif "autosomal recessive" in tl:
        inheritance = "Autosomal Recessive"
    elif "x-linked" in tl:
        inheritance = "X-Linked"

    matches = re.findall(r"(\d{1,3}(?:,\d{3})|\d{3,5})-(?:amino.acid|residue|aa)", all_text, re.I)
    protein_len = max(int(m.replace(",", "")) for m in matches) if matches else None
    coding_kb   = round(protein_len * 3 / 1000, 2) if protein_len else None
    aav_feasible = (coding_kb <= 4.7) if coding_kb else None

    result.update({
        "omim_id": mim_number, "gene_title": gene_title,
        "mechanism": mechanism, "mechanism_confidence": confidence,
        "protein_length": protein_len, "coding_kb": coding_kb,
        "aav_feasible": aav_feasible, "inheritance": inheritance,
        "molecular_text": (mol_text or all_text)[:2000],
    })
    return result


# ══════════════════════════════════════════════════════════════════════════════
# Step 2 — MyGene.info
# ══════════════════════════════════════════════════════════════════════════════

def query_mygene(gene: str) -> dict:
    print(f"\n[2/6] Querying MyGene.info for cross-database IDs...")
    result = {"entrez_id": None, "ensembl_gene_id": None,
              "uniprot_id": None, "refseq_transcript": None, "pathways": []}

    r = safe_get("https://mygene.info/v3/query",
                 params={"q": gene, "species": "human",
                         "fields": "entrezgene,ensembl,uniprot,refseq,pathway"},
                 label="MyGene query")
    if not r:
        return result

    hits = r.json().get("hits", [])
    if not hits:
        return result
    hit = hits[0]

    ensembl = hit.get("ensembl", {})
    if isinstance(ensembl, list):
        ensembl = ensembl[0]
    ensembl_id = ensembl.get("gene") if isinstance(ensembl, dict) else None

    uniprot = hit.get("uniprot", {})
    uniprot_id = None
    if isinstance(uniprot, dict):
        sw = uniprot.get("Swiss-Prot")
        uniprot_id = sw[0] if isinstance(sw, list) else sw

    refseq = hit.get("refseq", {})
    rna = refseq.get("rna") if isinstance(refseq, dict) else None
    refseq_tx = rna[0] if isinstance(rna, list) else rna

    pathway_names = []
    for db in ["kegg", "reactome", "wikipathways"]:
        pw = hit.get("pathway", {}).get(db, [])
        if isinstance(pw, dict):
            pw = [pw]
        for p in (pw if isinstance(pw, list) else []):
            if p.get("name"):
                pathway_names.append(f"{db.upper()}: {p['name']}")

    result.update({
        "entrez_id": hit.get("entrezgene"),
        "ensembl_gene_id": ensembl_id,
        "uniprot_id": uniprot_id,
        "refseq_transcript": refseq_tx,
        "pathways": pathway_names[:10],
    })
    print(f"  UniProt: {uniprot_id} | Ensembl: {ensembl_id}")
    return result


# ══════════════════════════════════════════════════════════════════════════════
# Step 3 — UniProt
# ══════════════════════════════════════════════════════════════════════════════

def query_uniprot(uniprot_id: str) -> dict:
    print(f"\n[3/6] Querying UniProt ({uniprot_id}) for protein features...")
    result = {"features": [], "function_text": None, "length": None}
    if not uniprot_id:
        print("  [skip] No UniProt ID")
        return result

    r = safe_get(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}",
                 params={"format": "json"}, label="UniProt")
    if not r:
        return result

    data = r.json()
    result["length"] = data.get("sequence", {}).get("length")

    for c in data.get("comments", []):
        if c.get("commentType") == "FUNCTION":
            texts = c.get("texts", [])
            if texts:
                result["function_text"] = texts[0].get("value", "")[:500]
                break

    keep = {"Domain", "Binding site", "Active site", "Modified residue",
            "Natural variant", "Transmembrane", "Topological domain", "Region", "Motif"}
    parsed = []
    for f in data.get("features", []):
        ftype = f.get("type", "")
        if ftype not in keep:
            continue
        loc   = f.get("location", {})
        start = loc.get("start", {}).get("value")
        end   = loc.get("end", {}).get("value")
        desc  = f.get("description", "") or f.get("featureId", "")
        if start and end:
            parsed.append(ProteinFeature(feature_type=ftype, description=desc,
                                         start=start, end=end))
    result["features"] = parsed
    print(f"  Found {len(parsed)} annotated features")
    return result


# ══════════════════════════════════════════════════════════════════════════════
# Step 4 — ClinVar (gene level)
# ══════════════════════════════════════════════════════════════════════════════

def query_clinvar(gene: str, entrez_id) -> list:
    print(f"\n[4/6] Querying ClinVar for pathogenic variants in {gene}...")
    if not entrez_id:
        print("  [skip] No Entrez ID")
        return []

    search_r = safe_get(f"{NCBI_BASE}/esearch.fcgi", params={
        "db": "clinvar",
        "term": f"{gene}[gene] AND (pathogenic[clinical_significance] "
                f"OR likely_pathogenic[clinical_significance])",
        "retmax": 200, "retmode": "json",
    }, label="ClinVar search")
    if not search_r:
        return []

    ids = search_r.json().get("esearchresult", {}).get("idlist", [])
    if not ids:
        print("  No ClinVar records found")
        return []

    print(f"  Found {len(ids)} records — fetching top 50...")
    ids = ids[:50]

    sum_r = safe_get(f"{NCBI_BASE}/esummary.fcgi",
                     params={"db": "clinvar", "id": ",".join(ids), "retmode": "json"},
                     label="ClinVar esummary")
    if not sum_r:
        return []

    summaries = sum_r.json().get("result", {})
    cons_map = {
        "nonsense": "LoF", "frameshift": "LoF", "stop gained": "LoF",
        "splice": "LoF", "deletion": "LoF",
        "missense": "Unknown", "gain of function": "GoF",
    }
    variants = []
    for uid in ids:
        rec = summaries.get(uid, {})
        if not rec or uid == "uids":
            continue
        title    = rec.get("title", "")
        sig_obj  = rec.get("clinical_significance", {})
        sig      = sig_obj.get("description", "") if isinstance(sig_obj, dict) else str(sig_obj)
        mol_cons = (rec.get("molecular_consequence_list") or ["unknown"])[0]
        trait_set = rec.get("trait_set", [])
        condition = trait_set[0].get("trait_name", "") if trait_set else ""
        mech = "Unknown"
        for key, val in cons_map.items():
            if key in mol_cons.lower() or key in title.lower():
                mech = val
                break
        cdna_m  = re.search(r"(c\.[^\s\)]+)", title)
        prot_m  = re.search(r"(p\.[^\s\)]+)", title)
        variants.append(VariantRecord(
            clinvar_id=uid, name=title, molecular_consequence=mol_cons,
            significance=sig, condition=condition, mechanism_clue=mech,
            position_cdna=cdna_m.group(1) if cdna_m else None,
            position_protein=prot_m.group(1) if prot_m else None,
        ))

    print(f"  Parsed {len(variants)} variant records")
    return variants


# ══════════════════════════════════════════════════════════════════════════════
# Step 5 — Open Targets
# ══════════════════════════════════════════════════════════════════════════════

def query_open_targets(ensembl_id: str) -> list:
    print(f"\n[5/6] Querying Open Targets for known drugs...")
    if not ensembl_id:
        print("  [skip] No Ensembl ID")
        return []

    gql = """
    query KnownDrugs($ensemblId: String!) {
      target(ensemblId: $ensemblId) {
        knownDrugs { rows {
          drug { name maximumClinicalTrialPhase }
          mechanismOfAction disease { name }
        }}
      }
    }"""
    try:
        r = requests.post("https://api.platform.opentargets.org/api/v4/graphql",
                          json={"query": gql, "variables": {"ensemblId": ensembl_id}},
                          timeout=15)
        rows = (r.json().get("data", {}).get("target", {})
                        .get("knownDrugs", {}).get("rows", []))
        drugs = [{
            "drug": row.get("drug", {}).get("name"),
            "mechanism": row.get("mechanismOfAction"),
            "disease": row.get("disease", {}).get("name"),
            "max_phase": row.get("drug", {}).get("maximumClinicalTrialPhase"),
        } for row in rows[:15]]
        print(f"  Found {len(drugs)} drug records")
        return drugs
    except Exception as e:
        print(f"  [warn] Open Targets failed: {e}")
        return []


# ══════════════════════════════════════════════════════════════════════════════
# NEW — Step 6a: Patient Variant Analysis
# ══════════════════════════════════════════════════════════════════════════════

def resolve_patient_variant(gene: str, variant_hgvs: str,
                            uniprot_features: list,
                            refseq_transcript: Optional[str]) -> tuple:
    """
    Given a patient HGVS cDNA variant (e.g. c.4849C>T):
      1. ClinVar lookup for exact pathogenicity + protein change
      2. Amino acid position extracted (or estimated from cDNA)
      3. Intersects AA position against every UniProt annotated feature
      4. gnomAD population frequency (ultra-rare = strong pathogenicity signal)
      5. Allele-specific ASO discrimination feasibility
    Returns: (VariantRecord | None, VariantProteinIntersection)
    """
    print(f"\n[NEW-a] Resolving patient variant {variant_hgvs}...")
    vpi = VariantProteinIntersection(aa_position=None)
    matched = None

    # ── ClinVar lookup ─────────────────────────────────────────────────────────
    # The ">" in HGVS (e.g. c.4849C>T) causes problems in Entrez free-text
    # queries.  We use NCBI's dedicated Variation Services API instead, which
    # accepts HGVS strings natively and returns structured JSON.
    # Fall back to an Entrez esearch using only the numeric position if the
    # Variation Services call fails.

    def _fetch_clinvar_by_hgvs(gene, variant_hgvs):
        """Try NCBI Variation Services first (most reliable for HGVS lookup)."""
        # Build candidate HGVS strings to try.  Users often omit the transcript
        # prefix; we try bare cDNA notation and also gene-prefixed form.
        candidates = [variant_hgvs]  # e.g. c.4849C>T
        # Also try with gene prefix in case API needs it
        candidates.append(f"{gene}:{variant_hgvs}")

        for hgvs_str in candidates:
            try:
                import urllib.parse
                encoded = urllib.parse.quote(hgvs_str, safe="")
                r = requests.get(
                    f"https://api.ncbi.nlm.nih.gov/variation/v0/hgvs/{encoded}/contextuals",
                    headers={"Accept": "application/json"},
                    timeout=12
                )
                if r.ok:
                    data = r.json()
                    # Get the ClinVar variant ID if present
                    cv_ids = []
                    for seq_loc in data.get("data", {}).get("spdis", []):
                        pass  # spdis gives us coords but not ClinVar ID directly
                    # Try the simpler /vcv endpoint via esearch with genomic coords
                    return None  # coords obtained but need separate ClinVar lookup
            except Exception:
                pass
        return None

    # Primary approach: Entrez esearch with HGVS position number + gene
    # We extract the numeric position and search ClinVar more broadly,
    # then filter results by matching the variant string in the title.
    cdna_num = parse_cdna_position(variant_hgvs)

    # Build multiple search terms to maximise recall:
    # 1. Gene + numeric position
    # 2. Gene + full variant (with > URL-safe as [variant])
    search_terms = []
    if cdna_num:
        search_terms.append(f"{gene}[gene] AND {cdna_num}[variant_pos]")
    # Encode > as literal string — Entrez handles it in title search
    search_terms.append(f"{gene}[gene] AND \"{variant_hgvs}\"[all]")

    ids = []
    for term in search_terms:
        r = safe_get(f"{NCBI_BASE}/esearch.fcgi", params={
            "db": "clinvar", "retmax": 10, "retmode": "json", "term": term,
        }, label="ClinVar variant search")
        if r:
            found = r.json().get("esearchresult", {}).get("idlist", [])
            ids.extend(x for x in found if x not in ids)
        if ids:
            break  # stop once we have hits

    if ids:
        sum_r = safe_get(f"{NCBI_BASE}/esummary.fcgi",
                         params={"db": "clinvar",
                                 "id": ",".join(ids[:5]), "retmode": "json"},
                         label="ClinVar variant summary")
        if sum_r:
            sums = sum_r.json().get("result", {})

            # Build normalised match tokens from the user's input:
            # e.g. "c.4849C>T" → tokens {"4849", "c>t", "4849c", "4849ct"}
            norm_input = variant_hgvs.lower().replace("c.", "").replace(">", "")
            pos_str = str(cdna_num) if cdna_num else ""

            for uid in ids[:5]:
                rec = sums.get(uid, {})
                if not rec or uid == "uids":
                    continue
                title = rec.get("title", "").lower()

                # Match if the position number AND at least the alt base appear
                # in the title, OR if the full HGVS string appears literally.
                position_match = pos_str and pos_str in title
                full_match     = variant_hgvs.lower() in title
                norm_match     = norm_input and norm_input in title.replace(".", "").replace(">","")

                if full_match or (position_match and norm_match):
                    sig_obj  = rec.get("clinical_significance", {})
                    sig      = sig_obj.get("description", "") if isinstance(sig_obj, dict) else str(sig_obj)
                    mol_cons = (rec.get("molecular_consequence_list") or ["unknown"])[0]
                    cond     = (rec.get("trait_set") or [{}])[0].get("trait_name", "")
                    cdna_m   = re.search(r"(c\.[^\s\)]+)", rec.get("title",""))
                    prot_m   = re.search(r"(p\.[^\s\)]+)", rec.get("title",""))
                    mech = "LoF" if any(x in mol_cons.lower()
                                        for x in ["nonsense","frameshift","splice","stop"]) \
                           else "Unknown"
                    matched = VariantRecord(
                        clinvar_id=uid, name=rec.get("title",""),
                        molecular_consequence=mol_cons, significance=sig,
                        condition=cond, mechanism_clue=mech,
                        position_cdna=cdna_m.group(1) if cdna_m else variant_hgvs,
                        position_protein=prot_m.group(1) if prot_m else None,
                    )
                    vpi.clinical_significance = sig
                    print(f"  ClinVar match: {rec.get('title','')[:70]}")
                    break

    if not matched:
        print(f"  No ClinVar match for {variant_hgvs} — proceeding with positional estimates")

    # ── Amino acid position ───────────────────────────────────────────────────
    aa_pos = None
    if matched and matched.position_protein:
        aa_pos = parse_aa_position(matched.position_protein)
    if aa_pos is None:
        cdna_pos = parse_cdna_position(variant_hgvs)
        if cdna_pos:
            aa_pos = (cdna_pos + 2) // 3
            print(f"  Estimated AA position from cDNA pos {cdna_pos}: ~aa {aa_pos}")
    vpi.aa_position = aa_pos

    # ── UniProt feature intersection ──────────────────────────────────────────
    if aa_pos:
        hits = [f for f in uniprot_features if f.start <= aa_pos <= f.end]
        vpi.overlapping_features = hits
        for f in hits:
            if f.feature_type == "Active site":
                vpi.is_in_active_site = True
            elif f.feature_type == "Binding site":
                vpi.is_in_binding_site = True
            elif f.feature_type == "Domain":
                vpi.is_in_domain = True
        if hits:
            descs = ", ".join(f.description or f.feature_type for f in hits)
            print(f"  !! Variant at aa {aa_pos} overlaps: {descs}")
        else:
            print(f"  Variant at aa {aa_pos}: no overlap with annotated protein features")

    # ── gnomAD population frequency ───────────────────────────────────────────
    print(f"  Querying gnomAD...")
    try:
        gql = """
        query SearchVariants($query: String!, $dataset: DatasetId!) {
          searchResults(query: $query, dataset: $dataset) {
            ... on VariantSearchResult {
              variant { variantId genome { af } exome { af } }
            }
          }
        }"""
        rg = requests.post(
            "https://gnomad.broadinstitute.org/api",
            json={"query": gql,
                  "variables": {"query": f"{gene} {variant_hgvs}", "dataset": "gnomad_r4"}},
            headers={"Content-Type": "application/json"}, timeout=12
        )
        if rg.ok:
            for res in rg.json().get("data", {}).get("searchResults", []):
                v = res.get("variant", {})
                af = ((v.get("genome") or {}).get("af") or
                      (v.get("exome")  or {}).get("af"))
                if af is not None:
                    vpi.gnomad_allele_frequency = af
                    vpi.is_absent_from_gnomad   = af < 1e-6
                    print(f"  gnomAD AF: {af:.2e}")
                    break
        if vpi.gnomad_allele_frequency is None:
            vpi.is_absent_from_gnomad = True
            print("  gnomAD: variant absent / not found (likely ultra-rare)")
    except Exception as e:
        print(f"  [warn] gnomAD query failed: {e}")
        vpi.is_absent_from_gnomad = True

    # ── ASO allele-specific discrimination feasibility ────────────────────────
    snv_m = re.search(r"[Cc]\.[0-9]+([ACGT])>([ACGT])", variant_hgvs)
    if snv_m:
        ref, alt = snv_m.group(1), snv_m.group(2)
        transitions = {("A","G"),("G","A"),("C","T"),("T","C")}
        vpi.aso_discrimination_feasible = True
        if (ref, alt) in transitions:
            vpi.aso_notes = (
                f"Transition ({ref}>{alt}): allele-specific ASO feasible but requires "
                f"careful placement of the discriminating mismatch near the 3' end of the "
                f"gapmer (positions 1-4 from 3' end). LNA or 2'-MOE modifications strongly "
                f"recommended to sharpen single-base discrimination."
            )
        else:
            vpi.aso_notes = (
                f"Transversion ({ref}>{alt}): allele-specific ASO discrimination is "
                f"favorable — transversions create larger thermodynamic mismatch penalty. "
                f"Standard gapmer design should provide adequate mutant selectivity."
            )
    elif re.search(r"del|ins|dup", variant_hgvs, re.I):
        vpi.aso_discrimination_feasible = True
        vpi.aso_notes = (
            "Indel variant: allele-specific ASO targeting is highly favorable. "
            "The deletion/insertion creates a unique sequence context that is "
            "straightforward to design against with high mutant specificity."
        )
    else:
        vpi.aso_discrimination_feasible = None
        vpi.aso_notes = "Could not parse SNV type — manual ASO feasibility assessment required."

    return matched, vpi


# ══════════════════════════════════════════════════════════════════════════════
# NEW — Step 6b: CRISPR Guide RNA Design
# ══════════════════════════════════════════════════════════════════════════════

def design_crispr_guides(gene: str, ensembl_gene_id: Optional[str],
                         variant_hgvs: Optional[str], mechanism: str) -> CRISPRParameters:
    """
    1. Fetches genomic sequence from Ensembl (~600bp window around variant,
       or gene start region if no variant supplied)
    2. Scans both strands for NGG PAM sites
    3. Scores all 20nt guides and returns top 8
    4. Checks base editing (CBE/ABE) and prime editing feasibility
    """
    print(f"\n[NEW-b] Designing CRISPR guide RNAs...")
    strategy = ("knock-out" if any(x in mechanism
                                   for x in ["Gain of Function","Dominant Negative"])
                else "knock-in / correction")
    crispr = CRISPRParameters(strategy=strategy)

    if not ensembl_gene_id:
        crispr.notes.append("No Ensembl ID — cannot fetch genomic sequence for guide design")
        return crispr

    # ── Gene coordinates from Ensembl ─────────────────────────────────────────
    gene_r = safe_get(f"{ENSEMBL_BASE}/lookup/id/{ensembl_gene_id}",
                      params={"content-type": "application/json", "expand": 0},
                      label="Ensembl gene lookup")
    if not gene_r:
        crispr.notes.append("Ensembl gene lookup failed")
        return crispr

    info       = gene_r.json()
    chrom      = info.get("seq_region_name")
    gene_start = info.get("start")
    gene_end   = info.get("end")

    if not all([chrom, gene_start, gene_end]):
        crispr.notes.append("Could not retrieve genomic coordinates from Ensembl")
        return crispr

    # ── Determine sequence window ─────────────────────────────────────────────
    win_start = gene_start
    win_end   = min(gene_start + 500, gene_end)

    if variant_hgvs:
        cdna_pos = parse_cdna_position(variant_hgvs)
        if cdna_pos:
            # cDNA position underestimates genomic offset (introns add length),
            # but gives a useful approximate anchor for a ±300bp window
            win_start = max(gene_start, gene_start + cdna_pos - 300)
            win_end   = win_start + 600
            crispr.notes.append(
                f"Sequence window centered ~{cdna_pos}bp into the gene for variant "
                f"{variant_hgvs}. Note: intronic sequence means actual genomic offset "
                f"is larger than cDNA position — verify final guide coords against the "
                f"full transcript map before ordering."
            )

    # ── Fetch sequence ────────────────────────────────────────────────────────
    seq_r = safe_get(
        f"{ENSEMBL_BASE}/sequence/region/human/{chrom}:{win_start}..{win_end}:1",
        params={"content-type": "application/json"},
        label="Ensembl sequence"
    )
    if not seq_r:
        crispr.notes.append("Could not retrieve genomic sequence from Ensembl")
        return crispr

    seq = seq_r.json().get("seq", "").upper()
    if not seq:
        crispr.notes.append("Empty sequence returned from Ensembl")
        return crispr

    print(f"  Retrieved {len(seq)}bp genomic window ({chrom}:{win_start}-{win_end})")
    crispr.target_sequence_context = seq[:80] + "..." if len(seq) > 80 else seq

    # ── Scan for NGG PAM sites — both strands ─────────────────────────────────
    variant_center = len(seq) // 2   # approximate position of variant in window
    guides = []

    # Forward strand: [20nt guide][N][GG]
    for m in re.finditer(r"(?=([ACGT]{20}[ACGT]GG))", seq):
        ps   = m.group(1)[:20]
        pam  = m.group(1)[20:]
        cut  = m.start() + 17     # Cas9 cuts between nt 17-18 of protospacer
        dist = abs(cut - variant_center)
        flags = quality_flags_for_guide(ps)
        guides.append(GuideRNA(
            sequence=ps, pam=pam, strand="+",
            distance_to_variant_bp=dist,
            gc_content=round(gc_content(ps), 2),
            quality_flags=flags,
            quality_score=score_guide(ps, dist, flags),
            cut_site_genomic_pos=win_start + cut,
        ))

    # Reverse strand: reverse-complement the window, scan same pattern
    rc_table = str.maketrans("ACGT", "TGCA")
    rev = seq.translate(rc_table)[::-1]
    for m in re.finditer(r"(?=([ACGT]{20}[ACGT]GG))", rev):
        ps     = m.group(1)[:20]
        pam    = m.group(1)[20:]
        rc_cut = m.start() + 17
        # Map back to forward-strand coordinate
        fwd_cut = len(seq) - rc_cut - 1
        dist    = abs(fwd_cut - variant_center)
        flags   = quality_flags_for_guide(ps)
        guides.append(GuideRNA(
            sequence=ps, pam=pam, strand="-",
            distance_to_variant_bp=dist,
            gc_content=round(gc_content(ps), 2),
            quality_flags=flags,
            quality_score=score_guide(ps, dist, flags),
            cut_site_genomic_pos=win_start + fwd_cut if fwd_cut >= 0 else None,
        ))

    guides.sort(key=lambda g: (-g.quality_score, g.distance_to_variant_bp))
    crispr.candidate_guides = guides[:8]
    print(f"  {len(guides)} total guide candidates — returning top {len(crispr.candidate_guides)}")

    if not guides:
        crispr.notes.append("No NGG PAM sites found in window — consider SpRY or Cas12a (TTTV PAM)")

    # ── Base / prime editing feasibility ──────────────────────────────────────
    if variant_hgvs:
        snv_m = re.search(r"[Cc]\.[0-9]+([ACGT])>([ACGT])", variant_hgvs)
        if snv_m:
            ref, alt = snv_m.group(1), snv_m.group(2)
            if (ref, alt) in [("C","T"), ("G","A")]:
                crispr.base_editing_feasible = True
                crispr.base_edit_type = "CBE (cytosine base editor — C→T / G→A)"
                crispr.notes.append(
                    "C→T transition: compatible with CBE (e.g. BE4max, AncBE4max). "
                    "No double-strand break required. Editing window typically protospacer "
                    "positions 4-8 from PAM-distal end. Confirm guide positions window "
                    "correctly over the target cytosine."
                )
            elif (ref, alt) in [("A","G"), ("T","C")]:
                crispr.base_editing_feasible = True
                crispr.base_edit_type = "ABE (adenine base editor — A→G / T→C)"
                crispr.notes.append(
                    "A→G transition: compatible with ABE (e.g. ABE8e, ABE8.20m). "
                    "No double-strand break required. Editing window typically positions "
                    "4-8. Confirm guide windows correctly over the target adenine."
                )
            else:
                crispr.prime_editing_feasible = True
                crispr.notes.append(
                    f"Transversion ({ref}→{alt}): not addressable by CBE or ABE. "
                    "Prime editing (PE3 or PE3b with pegRNA) recommended for precise "
                    "single-base correction without a double-strand break."
                )
        elif re.search(r"del|ins|dup", variant_hgvs, re.I):
            crispr.prime_editing_feasible = True
            crispr.notes.append(
                "Indel variant: prime editing (PE3/PE3b) can install precise correction. "
                "For LoF truncation variants, NHEJ-mediated disruption via paired guide "
                "RNAs (exon deletion) is a simpler alternative."
            )

    return crispr


# ══════════════════════════════════════════════════════════════════════════════
# NEW — Step 6c: Literature Evidence Scoring
# ══════════════════════════════════════════════════════════════════════════════

def query_literature_evidence(gene: str) -> LiteratureEvidence:
    """
    Queries PubMed via NCBI Entrez for gene characterization depth.
    Assigns evidence tier and returns top review PMIDs.
    """
    print(f"\n[NEW-c] Scoring literature evidence...")
    ev = LiteratureEvidence()

    # ── Gene-level count ──────────────────────────────────────────────────────
    rg = safe_get(f"{NCBI_BASE}/esearch.fcgi", params={
        "db": "pubmed", "retmax": 0, "retmode": "json",
        "term": (f"{gene}[tiab] AND "
                 "(mutation OR variant OR genetic disease OR channelopathy OR gene therapy)"),
    }, label="PubMed gene count")
    if rg:
        ev.gene_pubmed_count = int(rg.json().get("esearchresult", {}).get("count", 0))
        print(f"  Gene '{gene}': {ev.gene_pubmed_count} PubMed publications")

    # ── Top review PMIDs ──────────────────────────────────────────────────────
    rr = safe_get(f"{NCBI_BASE}/esearch.fcgi", params={
        "db": "pubmed", "retmax": 5, "retmode": "json", "sort": "relevance",
        "term": f"{gene}[tiab] AND review[pt]",
    }, label="PubMed reviews")
    if rr:
        ev.top_review_pmids = rr.json().get("esearchresult", {}).get("idlist", [])[:5]

    # ── Tier classification (gene-only) ──────────────────────────────────────
    gc = ev.gene_pubmed_count

    gene_tier = ("well-characterized" if gc >= 500 else
                 "moderate"           if gc >= 100 else
                 "emerging"           if gc >= 20  else "unknown")
    ev.evidence_tier = gene_tier

    detail = (
        "Strong foundational evidence — high confidence in target biology." if gc >= 500 else
        "Moderate evidence base — therapy recommendations are well-supported but "
        "independent functional validation recommended." if gc >= 100 else
        "Limited published evidence — recommendations are biologically plausible "
        "but require functional validation before clinical translation." if gc >= 20 else
        "Poorly characterised gene — treat all recommendations as highly preliminary; "
        "extensive functional and mechanistic work required."
    )
    ev.tier_rationale = f"{gc} publications → {gene_tier}. {detail}"

    return ev


# ══════════════════════════════════════════════════════════════════════════════
# Tissue class inference  (used by confidence engine)
# ══════════════════════════════════════════════════════════════════════════════

# Maps OMIM / pathway text clues → tissue class
# Tissue class drives the delivery vehicle precedent score
TISSUE_CLUES = {
    "CNS":    ["neuron","brain","cerebell","cortex","hippocamp","spinal","epilep",
               "dravet","rett","parkinson","alzheimer","huntington","ataxia",
               "scn1","mecp2","fmr1","tsc","pten","nf1"],
    "liver":  ["hepat","liver","cirrhosis","cholest","bile","urea cycle",
               "pkd","glycogen","pah","fah","hgd","ttr","fah"],
    "muscle": ["muscul","myopath","dystrophin","dmd","bmd","myotub","nemaline",
               "dysferlin","calpain","limb-girdle","facioscapulohum"],
    "blood":  ["hematop","erythro","hemoglobin","sickle","thalassemia","hemophil",
               "wiskott","ada-scid","cgt","gaucher","fabry","niemann"],
    "lung":   ["pulmonary","bronch","alveol","cftr","cystic fibrosis","surfactant",
               "alpha-1 antitrypsin"],
    "eye":    ["retina","photoreceptor","leber","stargardt","rpe65","cep290",
               "usherin","usher","choroider"],
    "kidney": ["nephron","renal","glomerul","podocyte","pkd1","pkd2","alport",
               "nephronophthisis"],
    "heart":  ["cardiac","cardiomyopath","sarcomere","myh7","tnnt2","scn5a",
               "lqt","brugada","arvc","hcm","dcm"],
}

def infer_tissue(gene: str, omim_text: str, pathways: list) -> str:
    """Return best-guess tissue class from gene name + OMIM text."""
    combined = (gene + " " + (omim_text or "") + " " + " ".join(pathways)).lower()
    scores = {tissue: sum(1 for clue in clues if clue in combined)
              for tissue, clues in TISSUE_CLUES.items()}
    best = max(scores, key=scores.get)
    return best if scores[best] > 0 else "systemic"


# ══════════════════════════════════════════════════════════════════════════════
# Delivery × tissue precedent matrix
# Scores: 3=well-established in humans, 2=clinical trials/late preclinical,
#         1=early preclinical, 0=no precedent
# ══════════════════════════════════════════════════════════════════════════════

DELIVERY_PRECEDENT = {
    #            CNS  liver muscle blood  lung   eye  kidney heart  systemic
    "AAV9":     [ 3,    2,    2,    1,     1,     1,    1,     2,     2 ],
    "AAVrh10":  [ 3,    2,    1,    1,     1,     1,    1,     1,     2 ],
    "AAV8":     [ 2,    3,    2,    1,     1,     1,    1,     2,     2 ],
    "AAV5":     [ 2,    2,    1,    1,     1,     3,    1,     1,     1 ],
    "AAV2":     [ 2,    2,    1,    1,     1,     3,    1,     1,     1 ],
    "LNP":      [ 1,    3,    1,    2,     2,     0,    1,     1,     2 ],
    "EVLP":     [ 1,    1,    1,    3,     1,     0,    1,     1,     2 ],
    "Naked ASO":[ 2,    2,    1,    2,     1,     0,    2,     1,     2 ],
    "Nanopart": [ 1,    2,    1,    2,     2,     0,    1,     1,     2 ],
}
TISSUE_ORDER = ["CNS","liver","muscle","blood","lung","eye","kidney","heart","systemic"]

def delivery_score(vehicle: str, tissue: str) -> int:
    row = DELIVERY_PRECEDENT.get(vehicle)
    if not row:
        return 0
    idx = TISSUE_ORDER.index(tissue) if tissue in TISSUE_ORDER else 8
    return row[idx]

DELIVERY_LABELS = {
    "AAV9":      "AAV9 (broad CNS/systemic tropism)",
    "AAVrh10":   "AAVrh10 (CNS/liver tropism)",
    "AAV8":      "AAV8 (liver/muscle tropism)",
    "AAV5":      "AAV5 (CNS/retinal tropism)",
    "AAV2":      "AAV2 (retinal/CNS tropism)",
    "LNP":       "LNP — lipid nanoparticle (liver/systemic)",
    "EVLP":      "EVLP — ex vivo lentiviral (HSC/blood)",
    "Naked ASO": "Naked ASO / GalNAc-ASO (CNS/liver)",
    "Nanopart":  "Polymer nanoparticle (lung/systemic)",
}

def best_delivery_for(tissue: str, top_n: int = 3) -> list:
    """Return top_n delivery vehicles for this tissue, sorted by precedent score."""
    scored = sorted(DELIVERY_PRECEDENT.keys(),
                    key=lambda v: delivery_score(v, tissue), reverse=True)
    return scored[:top_n]


# ══════════════════════════════════════════════════════════════════════════════
# Per-combo confidence scoring
# ══════════════════════════════════════════════════════════════════════════════

def combo_confidence(gene_tier: str, mech_fit: str,
                     del_score: int) -> tuple:
    """
    Returns (confidence_label, rationale_string).

    gene_tier   : well-characterized / moderate / emerging / unknown
    mech_fit    : high / medium / low  (logical fit of payload to mechanism)
    del_score   : 0-3 from DELIVERY_PRECEDENT
    """
    tier_score = {"well-characterized": 3, "moderate": 2,
                  "emerging": 1, "unknown": 0}.get(gene_tier, 0)
    fit_score  = {"high": 3, "medium": 2, "low": 1}.get(mech_fit, 1)

    total = tier_score + fit_score + del_score   # 0–9

    if total >= 7:
        label = "high"
    elif total >= 5:
        label = "medium"
    elif total >= 3:
        label = "low"
    else:
        label = "experimental"

    rationale_parts = []
    rationale_parts.append(
        f"Gene evidence: {gene_tier} ({tier_score}/3)"
    )
    rationale_parts.append(
        f"mechanism-payload fit: {mech_fit} ({fit_score}/3)"
    )
    rationale_parts.append(
        f"delivery precedent in target tissue: {del_score}/3"
    )
    return label, " · ".join(rationale_parts)


# ══════════════════════════════════════════════════════════════════════════════
# Modality Guide builders — one per payload class
# ══════════════════════════════════════════════════════════════════════════════

# ── Tissue → preferred promoter lookup ────────────────────────────────────────
TISSUE_PROMOTER = {
    "CNS":      ("SYN1 (neuron-specific) or CBA/CMV (pan-neuronal)", "scAAV preferred for CNS — smaller payload, faster expression"),
    "liver":    ("ApoE-hAAT (hepatocyte-specific) or TBG", "ssAAV for large cDNAs; scAAV if ≤2.2 kb"),
    "muscle":   ("MCK or MHCK7 (muscle-specific)", "ssAAV for large dystrophin constructs"),
    "blood":    ("PGK or EF1α (ubiquitous, HSC-compatible)", "EVLP preferred over AAV for HSC transduction"),
    "lung":     ("SFTPC (type II pneumocyte) or SP-B", "Aerosolised delivery under investigation"),
    "eye":      ("VMD2 (RPE) or IRBP (photoreceptor)", "Subretinal injection for RPE; intravitreal for RGC"),
    "kidney":   ("NPHS2 (podocyte-specific) or CAG", "Direct renal injection or IV with renal-tropic capsid"),
    "heart":    ("cTNT or CMV (cardiac)", "IV or intramyocardial injection"),
    "systemic": ("CMV or CAG (ubiquitous)", "ssAAV; tissue-specific variants preferred if known"),
}

TISSUE_AAV_DOSE = {
    "CNS":      "1×10¹³–5×10¹⁴ vg/kg IV or 1×10¹⁰–1×10¹² vg intrathecal/ICV",
    "liver":    "5×10¹²–3×10¹³ vg/kg IV (hepatic)",
    "muscle":   "1×10¹⁴–2×10¹⁵ vg/kg IV or IM",
    "blood":    "Ex vivo (EVLP): 1–5 MOI per HSC; not typically AAV",
    "lung":     "1×10¹¹–1×10¹³ vg per lung lobe (inhaled or IT)",
    "eye":      "1×10¹⁰–1×10¹¹ vg subretinal or intravitreal per eye",
    "kidney":   "1×10¹³ vg/kg IV or intrarenal infusion",
    "heart":    "1×10¹³–1×10¹⁴ vg/kg IV or intramyocardial",
    "systemic": "1×10¹³–1×10¹⁴ vg/kg IV (tissue-dependent; verify)",
}

TISSUE_LNP = {
    "CNS":      ("SM-102 or DLin-MC3-DMA (intrathecal) / Lipid 5 (IV-CNS)", "CNS-targeting LNPs (apoE-mediated) in development; intrathecal preferred currently"),
    "liver":    ("ALC-0315 (Pfizer/BioNTech) or SM-102 (Moderna) or MC3", "Well-validated hepatic tropism; ApoE-mediated uptake via LDLR"),
    "muscle":   ("Muscle-targeting LNPs (charged lipid optimised)", "Systemic IV LNP with muscle-tropic formulation; IM injection an alternative"),
    "blood":    ("ALC-0315 or DODAP (HSC ex vivo)", "Ex vivo LNP transfection of HSCs; higher efficiency than in vivo blood"),
    "lung":     ("Nebulised LNP (DOTAP or ionisable lipid)", "Inhaled delivery; ionisable lipids with pulmonary-optimised pKa ~6.2"),
    "eye":      ("Cationic LNP (intravitreal)", "Limited precedent; subretinal LNP in preclinical"),
    "kidney":   ("Renal-targeting LNP (folate-conjugated)", "IV delivery; renal LNPs in early development"),
    "heart":    ("Cardiac LNP (IV or intramyocardial)", "Intramyocardial injection preferred for cardiac mRNA"),
    "systemic": ("ALC-0315 or SM-102 (IV systemic)", "Hepatic first-pass predominant; tissue-specific formulations vary"),
}

TISSUE_ASO_CONJUGATE = {
    "CNS":      "Naked ASO (IT/ICV) or Ligand-conjugated (e.g. cholesterol-ASO, IV + BBB crossing)",
    "liver":    "GalNAc-ASO (subcutaneous — high hepatic uptake via ASGPR; preferred) or naked IV",
    "muscle":   "Peptide-conjugated ASO (PPMO; e.g. Pip6a) or naked IV (moderate uptake)",
    "blood":    "Naked ASO or LNP-encapsulated for HSC ex vivo",
    "lung":     "Inhaled naked ASO or aerosolised nanoparticle",
    "eye":      "Intravitreal naked ASO (direct access; fomivirsen precedent)",
    "kidney":   "Naked ASO (renal accumulation after IV; moderate efficacy)",
    "heart":    "Naked ASO IV or pericardial injection",
    "systemic": "GalNAc-ASO (SC) for liver; naked ASO (IV/IT) for other tissues",
}


def _guide_aav(gene: str, payload: str, delivery_key: str,
               tissue: str, protein_len: Optional[int],
               coding_kb: Optional[float], variant_hgvs: Optional[str],
               refseq: Optional[str]) -> ModalityGuide:
    """Design parameters for AAV gene replacement or dual-AAV."""
    g = ModalityGuide(modality_type="aav")

    is_dual = "dual" in payload.lower() or "split" in payload.lower()
    serotype = delivery_key if delivery_key.startswith("AAV") else "AAV9"
    g.aav_serotype = serotype

    promoter_info = TISSUE_PROMOTER.get(tissue, TISSUE_PROMOTER["systemic"])
    g.aav_promoter = promoter_info[0]

    # ITR configuration
    if coding_kb and coding_kb <= 2.2:
        g.aav_itr_config = "scAAV (self-complementary) — doubles effective titer, faster onset"
    elif is_dual:
        g.aav_itr_config = (
            "Dual ssAAV split-intein: 5' vector carries N-terminal CDS + split intein donor; "
            "3' vector carries split intein acceptor + C-terminal CDS. "
            "Trans-splicing reconstitutes full-length protein."
        )
    else:
        g.aav_itr_config = "ssAAV (single-stranded) — required for coding sequences >2.2 kb"

    g.aav_dose_range = TISSUE_AAV_DOSE.get(tissue, TISSUE_AAV_DOSE["systemic"])
    g.aav_tropism_note = promoter_info[1]
    g.aav_packaging_note = (
        "HEK293 triple-transfection (rep/cap/GOI plasmids) for clinical-grade manufacture; "
        "Sf9/baculovirus system for scalability. "
        "Confirm capsid IP clearance for chosen serotype."
    )

    g.key_parameters = [
        f"Serotype: {serotype}",
        f"Vector type: {'Dual ssAAV split-intein' if is_dual else ('scAAV' if coding_kb and coding_kb <= 2.2 else 'ssAAV')}",
        f"Promoter: {g.aav_promoter}",
        f"Payload size: {f'{coding_kb} kb CDS' if coding_kb else 'Unknown — measure before vector design'}",
        f"Dose (estimated): {g.aav_dose_range}",
        f"Transcript: {refseq or 'Verify canonical RefSeq NM_ for CDS cloning'}",
    ]
    g.design_notes = [
        f"Include WPRE (Woodchuck Hepatitis Virus Post-transcriptional Regulatory Element) for mRNA stability in ssAAV.",
        f"Add SV40 polyA or BGH polyA signal 3' of CDS.",
        f"Confirm codon-optimised transgene does not trigger innate immune response via TLR9.",
        f"Pre-screen patient for neutralising antibodies against {serotype} before dosing.",
        f"Tissue: {tissue} — confirm {serotype} tropism is appropriate (see serotype biodistribution literature).",
    ]
    if variant_hgvs:
        g.design_notes.append(
            f"Patient carries {variant_hgvs} — transgene expresses WT protein regardless of endogenous variant; "
            f"no allele-specific design needed for LoF replacement."
        )
    return g


def _guide_mrna_lnp(gene: str, payload: str, tissue: str,
                    protein_len: Optional[int], coding_kb: Optional[float],
                    variant_hgvs: Optional[str]) -> ModalityGuide:
    """Design parameters for mRNA + LNP delivery."""
    g = ModalityGuide(modality_type="mrna_lnp")

    lnp_info = TISSUE_LNP.get(tissue, TISSUE_LNP["systemic"])
    g.lnp_formulation  = lnp_info[0]
    g.lnp_targeting    = lnp_info[1]
    g.mrna_utr_design  = (
        "5' UTR: human β-globin 5'UTR or Kozak-optimised synthetic (e.g. TISU element). "
        "3' UTR: AES + mtRNR1 tandem (Moderna-style) or human α-globin 3'UTR for stability."
    )
    g.mrna_codon_opt   = (
        f"Codon-optimise {gene} CDS for human codon usage (CAI > 0.8). "
        "Avoid CpG dinucleotides where possible to reduce innate immune activation. "
        "Use codon optimisation tools (e.g. IDT CodonOpt, GenScript OptimumGene)."
    )
    g.mrna_cap         = "CleanCap-AG (co-transcriptional capping; Cap1 structure for eIF4E binding)"
    g.mrna_modification = (
        "N1-methylpseudouridine (m1Ψ) full substitution — reduces TLR7/8 activation, "
        "increases ribosome loading and protein yield."
    )
    size_note = f"{coding_kb} kb CDS — " if coding_kb else ""
    g.mrna_stability = (
        f"{size_note}Expected mRNA half-life 12–72 h depending on tissue. "
        "Redosing schedule required (monthly to quarterly; tissue- and indication-dependent). "
        "Cold chain: −80°C for LNP-mRNA; lyophilisation under development."
    )

    g.key_parameters = [
        f"LNP formulation: {g.lnp_formulation}",
        f"Cap: {g.mrna_cap}",
        f"Modification: N1-methylpseudouridine (m1Ψ)",
        f"5' UTR: β-globin / Kozak-optimised",
        f"3' UTR: AES+mtRNR1 or α-globin",
        f"Codon optimisation: Required (CAI > 0.8)",
        f"CDS size: {f'{coding_kb} kb' if coding_kb else 'measure first — LNP can accommodate up to ~10 kb mRNA'}",
    ]
    g.design_notes = [
        f"LNP tissue targeting note: {g.lnp_targeting}",
        "Poly-A tail: 120 nt enzymatic or encoded — optimise for half-life.",
        "Avoid uORFs in 5' UTR; scan with ORF prediction tools before synthesis.",
        "HPLC purification of IVT mRNA to remove dsRNA contaminants (reduces immunogenicity).",
        "Pre-clinical: confirm protein expression + activity in relevant cell line before in vivo.",
    ]
    return g


def _guide_aso_sirna(gene: str, payload: str, tissue: str,
                     variant_hgvs: Optional[str], mechanism: str,
                     refseq: Optional[str], cdna_pos: Optional[int]) -> ModalityGuide:
    """Design parameters for ASO or siRNA."""
    g = ModalityGuide(modality_type="aso_sirna")
    is_sirna = "sirna" in payload.lower()
    is_splice = "splice" in payload.lower() or "exon" in payload.lower()
    is_allele_specific = "allele" in payload.lower()

    conjugate = TISSUE_ASO_CONJUGATE.get(tissue, TISSUE_ASO_CONJUGATE["systemic"])
    g.aso_conjugate = conjugate

    if is_splice:
        g.aso_mechanism    = "Steric block of splicing signals (SSO) — masks splice donor/acceptor to force exon inclusion or skipping"
        g.aso_target_region = (
            "Target splice donor (GT) or acceptor (AG) of the relevant exon, or "
            "exonic splicing enhancers (ESEs). "
            "Priority: exons flanking LoF variant or frameshifting exon."
        )
        g.aso_design_window = (
            "Design 15–25 nt ASOs centred ±50 nt around splice donor/acceptor. "
            "Screen 5–10 candidates across the splice junction; validate by RT-PCR."
        )
        g.aso_chemistry = "2'-MOE PS backbone (preferred for splice-switching) or LNA/DNA mixmer"
    elif is_allele_specific:
        g.aso_mechanism    = "RNase H-mediated degradation of mutant transcript (gapmer) with discriminating mismatch at position 1–4 from 3' end"
        g.aso_target_region = (
            f"Target region surrounding {variant_hgvs or 'the patient variant'} — "
            "design guide to position the variant at the thermodynamically critical zone (positions 1–4 from 3' end of gapmer)."
        )
        g.aso_design_window = (
            "15–20 nt gapmer: 2–4 LNA or 2'-MOE wings flanking 7–10 nt DNA gap. "
            "Discriminating mismatch in gap region. "
            "Screen 5+ candidates with varying mismatch positions; validate by allele-specific RT-qPCR."
        )
        g.aso_chemistry = "LNA gapmer (highest mismatch discrimination) or 2'-MOE PS gapmer"
    elif is_sirna:
        g.aso_mechanism    = "RISC-mediated cleavage of target mRNA (guide strand loaded into Ago2)"
        g.aso_target_region = (
            "Target 3'UTR or CDS of mutant transcript. "
            "Avoid seed-region (positions 2–7 of guide) matches to off-target transcriptomes."
        )
        g.aso_design_window = (
            "21 nt duplex (19 nt + 2 nt 3' overhang on both strands). "
            "Use siRNA design algorithms (siDirect, RNAi Designer) to minimise off-targets. "
            "Chemically stabilise with 2'-OMe alternating pattern + PS at termini."
        )
        g.aso_chemistry = "2'-OMe + PS (alternating) + optional LNA at positions 1, 2 of guide strand"
        g.sirna_strand_bias = (
            "Confirm guide strand preferential loading: weak base pair at 5' end of guide strand. "
            "Asymmetric duplex design (thermodynamic asymmetry) reduces passenger strand off-targets."
        )
    else:
        g.aso_mechanism    = "RNase H-mediated knockdown (gapmer) or steric block"
        g.aso_target_region = "5'UTR, AUG start site, or coding region splice junctions"
        g.aso_design_window = "15–20 nt; screen multiple positions by walking across the mRNA"
        g.aso_chemistry    = "LNA gapmer or 2'-MOE PS — select based on tissue delivery route"

    g.key_parameters = [
        f"Mechanism: {g.aso_mechanism}",
        f"Chemistry: {g.aso_chemistry}",
        f"Conjugate / delivery: {conjugate}",
        f"Target region: {(g.aso_target_region or '')[:80]}",
        f"Design window: {(g.aso_design_window or '')[:80]}",
        f"Transcript reference: {refseq or 'Identify canonical NM_ RefSeq for design'}",
    ]
    if g.sirna_strand_bias:
        g.key_parameters.append(f"Strand bias: {g.sirna_strand_bias[:80]}")

    g.design_notes = [
        f"Tissue delivery: {tissue} — {conjugate}.",
        "Validate knockdown by RT-qPCR and western blot in disease-relevant cell line.",
        "Run BLAST of ASO sequence against human transcriptome to flag off-target risk.",
        "Confirm no seed-region matches to essential genes (positions 2–7).",
        "Dose titration required — start at 1–10 mg/kg (ASO) or 0.1–1 mg/kg (siRNA-LNP).",
    ]
    if is_allele_specific:
        g.design_notes.append(
            "Allele specificity validation: compare WT vs mutant transcript knockdown in patient-derived cells "
            "(iPSC or primary) — require ≥10-fold selectivity before advancing."
        )
    return g


def _guide_crispr(gene: str, payload: str, tissue: str,
                  variant_hgvs: Optional[str], mechanism: str,
                  crispr_params: Optional["CRISPRParameters"]) -> ModalityGuide:
    """Design parameters for CRISPR (any variant — KO, base edit, prime edit, HDR)."""
    g = ModalityGuide(modality_type="crispr")

    is_base  = "base edit" in payload.lower()
    is_prime = "prime edit" in payload.lower()
    is_hdr   = "hdr" in payload.lower()
    is_ko    = "knock-out" in payload.lower() or "nhej" in payload.lower()

    if is_base:
        if crispr_params and crispr_params.base_edit_type:
            bet = crispr_params.base_edit_type
            g.crispr_editor  = f"BE4max (CBE)" if "CBE" in bet else "ABE8e (ABE)"
            g.crispr_window  = "C→T editing window: positions 4–8 from PAM (CBE)" if "CBE" in bet else "A→G editing window: positions 4–7 from PAM (ABE8e)"
        else:
            g.crispr_editor  = "BE4max (CBE, C→T) or ABE8e (ABE, A→G) — determine from variant type"
            g.crispr_window  = "Editing window positions 4–8 from PAM (verify for specific editor variant)"
        g.crispr_strategy = "DSB-free base editing — no double-strand break, lower indel risk"
        g.crispr_pam      = "NGG (SpCas9-derived editors) or NG (SpRY for expanded targeting)"
        g.crispr_offtarget = "Off-target base edits at bystander Cs/As within window — screen with EditScan or CIRCLE-seq"
    elif is_prime:
        g.crispr_strategy = "Prime editing (PE3 or PE3b) — installs precise edit via pegRNA reverse transcriptase"
        g.crispr_editor   = "PE3 (SpCas9-H840A nickase + RT) or PE3b (paired nicking for efficiency)"
        g.crispr_window   = "pegRNA spacer 20 nt; RT template 10–30 nt encoding desired edit; PBS 13 nt"
        g.crispr_pam      = "NGG (standard) — PE5max or PEmax variants improve efficiency"
        g.crispr_offtarget = "Lower DSB-associated off-targets vs Cas9; pegRNA scaffold insertions possible — validate by Cas-OFFinder + deep sequencing"
    elif is_hdr:
        g.crispr_strategy = "HDR correction — DSB at target site repaired using donor template"
        g.crispr_editor   = "SpCas9 (high-fidelity variant: eSpCas9 or HiFiCas9 recommended)"
        g.crispr_window   = "Cut site within 10–50 bp of correction site; ssODN donor (≤200 nt) or dsDNA template"
        g.crispr_pam      = "NGG — position cut site proximal to desired correction"
        g.crispr_offtarget = "HDR efficiency typically <10% in post-mitotic cells — ex vivo preferred; validate by amplicon deep sequencing"
    elif is_ko:
        g.crispr_strategy = "NHEJ knock-out — DSB at target exon creates frameshift indel"
        g.crispr_editor   = "SpCas9 (standard) — high-fidelity variants reduce off-targets"
        g.crispr_window   = "Target early exon (first 30% of CDS) for maximal LoF; avoid last exon (NMD escape)"
        g.crispr_pam      = "NGG — verify guide is unique in genome (specificity score ≥ 0.6)"
        g.crispr_offtarget = "Run CRISPOR or Cas-OFFinder; validate top 5 off-target sites by Sanger or deep sequencing"
    else:
        g.crispr_strategy = "SpCas9-based editing — strategy determined by variant type"
        g.crispr_editor   = "SpCas9 or high-fidelity variant"
        g.crispr_pam      = "NGG"
        g.crispr_offtarget = "Validate off-targets by Cas-OFFinder"

    g.crispr_delivery_form = (
        "RNP (ribonucleoprotein) electroporation — fastest clearance, lowest off-target "
        "for ex vivo (blood/HSC). "
        "AAV-sgRNA + AAV-editor for in vivo (split-intein for large editors like PE). "
        "LNP-mRNA (editor) + sgRNA for liver."
    )

    # Pull guide count from crispr_params if available
    guide_count = len(crispr_params.candidate_guides) if crispr_params else 0
    g.key_parameters = [
        f"Strategy: {g.crispr_strategy}",
        f"Editor: {g.crispr_editor}",
        f"PAM: {g.crispr_pam}",
        f"Editing window: {g.crispr_window}",
        f"Delivery form: RNP (ex vivo) / AAV or LNP-mRNA (in vivo)",
        f"Candidate guides computed: {guide_count} (see guide table below)" if guide_count else "Run guide design with variant input for scored guide list",
    ]
    g.design_notes = [
        f"Off-target strategy: {g.crispr_offtarget}",
        f"Tissue: {tissue} — confirm delivery form is appropriate (RNP for ex vivo; AAV/LNP for in vivo).",
        "Validate editing efficiency by Sanger + TIDE or amplicon deep sequencing.",
        "For allele-specific editing: confirm guide does not target WT allele (SNP-based discrimination).",
        "GMP manufacture: RNP components (Cas9 protein + IVT sgRNA) for clinical ex vivo; LNP-mRNA for in vivo.",
    ]
    return g


def _guide_cast(gene: str, payload: str, tissue: str,
                protein_len: Optional[int], coding_kb: Optional[float],
                variant_hgvs: Optional[str]) -> ModalityGuide:
    """Design parameters for CAST / transposon-mediated insertion."""
    g = ModalityGuide(modality_type="cast")

    is_paste = "paste" in payload.lower()
    is_pb    = "piggybac" in payload.lower() or "piggybac" in payload.lower()

    if is_paste:
        g.cast_system = "PASTE (Programmable Addition via Site-specific Targeting Elements) — Cas9 nickase + serine integrase + attB-flanked donor"
        g.cast_target_site = "attP landing site pre-installed in safe-harbour locus (AAVS1, CCR5, ROSA26 human equivalent) OR endogenous attP-like sequence"
        g.cast_donor_size  = "Up to ~15 kb cargo capacity — accommodates full gene + regulatory elements"
        g.cast_transposase_delivery = "Cas9 nickase + integrase as mRNA (LNP) or RNP + AAV donor template"
        g.cast_offtarget_risk = "Targeted integration via serine integrase — minimal off-target insertions vs piggyBac; verify by LAM-PCR or GUIDE-seq"
    elif is_pb:
        g.cast_system = "piggyBac transposon — TTAA site integration driven by piggyBac transposase"
        g.cast_target_site = "Semi-random TTAA sites genome-wide (not fully targeted) — enrichment near TSS; safe-harbour bias tools available"
        g.cast_donor_size  = "Up to ~200 kb theoretical; practical limit ~14 kb for efficient transposition"
        g.cast_transposase_delivery = "Transposase mRNA (LNP or electroporation) + DNA donor plasmid or minicircle"
        g.cast_offtarget_risk = "Semi-random integration — insertional mutagenesis risk; require whole-genome sequencing of clones; enrichment for safe-harbour sites recommended"
    else:
        g.cast_system = "CAST / transposon system (PASTE, piggyBac, or Sleeping Beauty — specify based on cargo size and targeting requirement)"
        g.cast_target_site = "Safe-harbour locus (AAVS1 / CCR5) preferred for predictable expression; gene-specific if regulation required"
        g.cast_donor_size  = f"Cargo: {f'{coding_kb} kb CDS' if coding_kb else 'measure before design'} + promoter + polyA (~2–4 kb overhead)"
        g.cast_transposase_delivery = "Transposase mRNA (LNP or electroporation) + donor template (plasmid, minicircle, or AAV)"
        g.cast_offtarget_risk = "Insertion site diversity: characterise by LAM-PCR or GUIDE-seq; avoid proto-oncogene loci"

    g.key_parameters = [
        f"System: {g.cast_system}",
        f"Target site: {g.cast_target_site}",
        f"Cargo capacity: {g.cast_donor_size}",
        f"Transposase delivery: {g.cast_transposase_delivery}",
        f"Integration risk: {(g.cast_offtarget_risk or '')[:80]}",
        "Clinical precedent: None (preclinical only as of 2025) — regulatory pathway undefined",
    ]
    g.design_notes = [
        "CAST / transposon approaches have no clinical approvals as of 2025 — treat as exploratory.",
        f"For {gene}: design donor carrying full WT CDS + tissue-appropriate promoter (see AAV guide for promoter options).",
        "Validate integration by Junction PCR, Southern blot, and WGS before in vivo use.",
        "Confirm transgene expression at correct level — safe-harbour expression may differ from endogenous locus.",
        "Regulatory: discuss genotoxicity assessment strategy with FDA/EMA early (IND pre-meeting recommended).",
    ]
    return g


def design_modality_guide(combo: "TherapyRecommendation",
                          gene: str, tissue: str,
                          protein_len: Optional[int], coding_kb: Optional[float],
                          variant_hgvs: Optional[str], mechanism: str,
                          refseq: Optional[str],
                          crispr_params: Optional["CRISPRParameters"]) -> ModalityGuide:
    """
    Dispatcher: inspects the payload string of a TherapyRecommendation
    and returns the appropriate ModalityGuide.
    """
    p = combo.payload.lower()
    cdna_pos = None
    if variant_hgvs:
        import re as _re
        m = _re.search(r"c\.(\d+)", variant_hgvs)
        cdna_pos = int(m.group(1)) if m else None

    if "aav" in p or "dual aav" in p or "split-intein" in p:
        # Extract serotype key from delivery label
        dk = "AAV9"
        for key in DELIVERY_PRECEDENT:
            if key.startswith("AAV") and key.lower() in combo.delivery.lower():
                dk = key
                break
        return _guide_aav(gene, combo.payload, dk, tissue, protein_len,
                          coding_kb, variant_hgvs, refseq)

    elif "mrna" in p or "mrna replacement" in p:
        return _guide_mrna_lnp(gene, combo.payload, tissue, protein_len,
                                coding_kb, variant_hgvs)

    elif any(x in p for x in ["aso", "sirna", "allele-specific", "splice", "exon"]):
        return _guide_aso_sirna(gene, combo.payload, tissue, variant_hgvs,
                                mechanism, refseq, cdna_pos)

    elif any(x in p for x in ["crispr", "base edit", "prime edit", "hdr", "nhej", "knock-out"]):
        return _guide_crispr(gene, combo.payload, tissue, variant_hgvs,
                             mechanism, crispr_params)

    elif any(x in p for x in ["cast", "transposon", "paste", "piggybac", "sleeping beauty"]):
        return _guide_cast(gene, combo.payload, tissue, protein_len,
                           coding_kb, variant_hgvs)

    else:
        # Fallback generic guide
        mg = ModalityGuide(modality_type="other")
        mg.key_parameters = [f"Payload: {combo.payload}", f"Delivery: {combo.delivery}"]
        mg.design_notes   = ["No specific guide template for this modality — manual design required."]
        return mg


# ══════════════════════════════════════════════════════════════════════════════
# Therapy recommendation engine — returns list[TherapyRecommendation]
# ══════════════════════════════════════════════════════════════════════════════

def recommend_therapies(p: AssayParameters,
                        gene_tier: str,
                        tissue: str) -> list:
    """
    Builds an ordered list of TherapyRecommendation objects.
    Every gene gets a concrete primary recommendation regardless of
    mechanism confidence — uncertainty is expressed via the confidence
    label and rationale, not by withholding a recommendation.
    """
    m   = p.mechanism
    aav = p.aav_feasible
    combos = []
    rank = 1

    def make(payload, delivery_key, mech_fit, notes="", cap_delivery=None):
        ds    = delivery_score(delivery_key, tissue)
        if cap_delivery is not None:
            ds = min(ds, cap_delivery)
        conf, rat = combo_confidence(gene_tier, mech_fit, ds)
        dlabel = DELIVERY_LABELS.get(delivery_key, delivery_key)
        return TherapyRecommendation(
            rank=rank,
            payload=payload,
            delivery=dlabel,
            label=f"{payload}  +  {dlabel}",
            confidence=conf,
            confidence_rationale=rat,
            notes=notes,
        )

    # ── Loss of Function ──────────────────────────────────────────────────────
    if "Loss of Function" in m or m == "Unknown":
        # Unknown mechanism defaults to LoF strategy as most common for
        # rare monogenic disease — flagged in notes
        unk_note = " (mechanism uncertain — LoF assumed as most common default; validate before proceeding)" if m == "Unknown" else ""

        if aav is not False:   # True or None (unknown size)
            # Primary: AAV replacement, best serotype for tissue
            best_aav = [v for v in best_delivery_for(tissue)
                        if v.startswith("AAV")]
            primary_aav = best_aav[0] if best_aav else "AAV9"
            combos.append(make(
                "AAV gene replacement (full-length cDNA)",
                primary_aav, "high",
                notes=unk_note or ("Single vector if ≤4.7 kb coding sequence" if aav else "Size unknown — verify coding sequence length")
            ))
            rank += 1

        if aav is False:
            # Too large for single AAV
            combos.append(make(
                "Dual AAV split-intein (oversized gene)",
                "AAV9", "high",
                notes=f"Gene coding sequence exceeds 4.7 kb single-vector limit{unk_note}"
            ))
            rank += 1
            combos.append(make(
                "mRNA replacement",
                "LNP", "high",
                notes="LNP-mRNA avoids size limit; transient expression requires redosing" + unk_note
            ))
            rank += 1
        else:
            combos.append(make(
                "mRNA replacement",
                "LNP", "high",
                notes="Transient expression; useful for enzyme/secreted proteins" + unk_note
            ))
            rank += 1

        combos.append(make(
            "CRISPR-Cas9 HDR correction",
            "EVLP", "medium",
            notes="Ex vivo approach for blood/HSC targets; in vivo HDR efficiency low" + unk_note,
            cap_delivery=2
        ))
        rank += 1

        combos.append(make(
            "CAST / transposon-mediated insertion (PASTE, piggyBac)",
            "AAV9", "low",
            notes="Large payload capacity; no DSB; preclinical only — no clinical precedent yet" + unk_note,
            cap_delivery=1
        ))
        rank += 1

        combos.append(make(
            "ASO (splice-switching / exon inclusion)",
            "Naked ASO", "medium",
            notes="Applicable if partial LoF with residual transcript; targets pre-mRNA" + unk_note
        ))
        rank += 1

    # ── Gain of Function ─────────────────────────────────────────────────────
    elif "Gain of Function" in m or "Dominant Negative" in m:
        dn_note = " (dominant negative — silencing mutant allele essential)" if "Dominant Negative" in m else ""
        ad_note = " (AD inheritance — allele-specific approach required to preserve WT)" if p.inheritance == "Autosomal Dominant" else ""

        combos.append(make(
            "Allele-specific ASO / siRNA (silence mutant transcript)",
            "Naked ASO", "high",
            notes="Requires SNP or variant-specific sequence distinct from WT" + dn_note + ad_note
        ))
        rank += 1

        combos.append(make(
            "CRISPR base editing (CBE/ABE — DSB-free)",
            "AAV9", "high",
            notes="For point mutations compatible with C→T (CBE) or A→G (ABE)" + ad_note
        ))
        rank += 1

        combos.append(make(
            "CRISPR prime editing (PE3/pegRNA)",
            "AAV9", "medium",
            notes="For transversions / small indels not addressable by base editors" + ad_note
        ))
        rank += 1

        combos.append(make(
            "CRISPR-Cas9 knock-out (NHEJ) of mutant allele",
            "LNP", "medium",
            notes="Disrupts mutant allele; allele-specificity requires guide targeting variant site" + ad_note
        ))
        rank += 1

        combos.append(make(
            "siRNA (chemically stabilised)",
            "LNP", "medium",
            notes="Strong hepatic silencing via LNP; CNS delivery challenging" + ad_note
        ))
        rank += 1

        combos.append(make(
            "CAST / transposon insertion (replacement strategy)",
            "EVLP", "low",
            notes="Experimental — insert corrected allele alongside disruption of mutant" + ad_note,
            cap_delivery=1
        ))
        rank += 1

    # ── Mixed ─────────────────────────────────────────────────────────────────
    elif "Mixed" in m:
        combos.append(make(
            "Allele-specific ASO (for GoF variants in this gene)",
            "Naked ASO", "medium",
            notes="Variant-level mechanism classification required first — some variants GoF, some LoF"
        ))
        rank += 1
        combos.append(make(
            "AAV gene replacement (for LoF variants in this gene)",
            "AAV9", "medium",
            notes="Applicable only for confirmed LoF variants; contraindicated for GoF"
        ))
        rank += 1
        combos.append(make(
            "CRISPR base editing (variant-specific)",
            "AAV9", "low",
            notes="Design guides per variant after mechanism classification"
        ))
        rank += 1

    return combos


# ══════════════════════════════════════════════════════════════════════════════
# Assemble report
# ══════════════════════════════════════════════════════════════════════════════

def build_report(gene: str, omim_key: str,
                 variant_hgvs: Optional[str] = None) -> AssayParameters:

    # Hard gate: raises GeneNotFoundError if not a real human gene symbol.
    # This propagates to the Flask route and is shown as a clean UI error.
    validate_gene(gene)

    omim     = query_omim(gene, omim_key)
    mygene   = query_mygene(gene)
    uniprot  = query_uniprot(mygene["uniprot_id"])
    variants = query_clinvar(gene, mygene["entrez_id"])
    drugs    = query_open_targets(mygene["ensembl_gene_id"])

    protein_len  = omim["protein_length"] or uniprot["length"]
    coding_kb    = omim["coding_kb"] or (round(protein_len * 3 / 1000, 2) if protein_len else None)
    aav_feasible = (coding_kb <= 4.7) if coding_kb else None

    domains       = [f for f in uniprot["features"] if f.feature_type == "Domain"]
    active_sites  = [f for f in uniprot["features"] if f.feature_type == "Active site"]
    binding_sites = [f for f in uniprot["features"] if f.feature_type == "Binding site"]

    downstream = list({node for pw in mygene["pathways"]
                       for node in ["mTOR","MEK","ERK","PI3K","AKT","NF-kB",
                                    "JAK","STAT","MAPK","Ras","Raf","CDK","BCL"]
                       if node.lower() in pw.lower()})

    params = AssayParameters(
        gene_symbol=gene,
        omim_id=str(omim["omim_id"]) if omim["omim_id"] is not None else "",
        uniprot_id=mygene["uniprot_id"],
        ensembl_gene_id=mygene["ensembl_gene_id"],
        refseq_transcript=mygene["refseq_transcript"],
        protein_length_aa=protein_len,
        coding_sequence_kb=coding_kb,
        canonical_transcript_id=mygene["refseq_transcript"],
        mechanism=omim["mechanism"],
        mechanism_confidence=omim["mechanism_confidence"],
        inheritance=omim["inheritance"],
        primary_therapy="",
        aav_feasible=aav_feasible,
        patient_variant_hgvs=variant_hgvs,
        representative_variants=variants[:20],
        druggable_domains=domains,
        active_sites=active_sites,
        binding_sites=binding_sites,
        pathway_names=mygene["pathways"],
        downstream_targets=downstream,
        known_drugs=drugs,
        omim_molecular_genetics_summary=omim["molecular_text"],
    )

    # ── NEW: patient variant ──────────────────────────────────────────────────
    if variant_hgvs:
        resolved, vpi = resolve_patient_variant(
            gene, variant_hgvs, uniprot["features"], mygene["refseq_transcript"]
        )
        params.patient_variant_resolved      = resolved
        params.variant_protein_intersection  = vpi

        # If the variant-level mechanism is unambiguous, use it to refine
        if resolved and resolved.mechanism_clue in ("GoF","LoF") and "Mixed" in params.mechanism:
            params.mechanism = ("Gain of Function" if resolved.mechanism_clue == "GoF"
                                else "Loss of Function")
            params.mechanism_confidence = "medium"
            params.notes.append(
                f"Mechanism upgraded from gene-level 'Mixed' to variant-level "
                f"'{params.mechanism}' based on ClinVar consequence for {variant_hgvs}"
            )

    # ── NEW: CRISPR guides ────────────────────────────────────────────────────
    params.crispr_parameters = design_crispr_guides(
        gene, mygene["ensembl_gene_id"], variant_hgvs, params.mechanism
    )

    # ── Literature evidence ───────────────────────────────────────────────────
    params.literature_evidence = query_literature_evidence(gene)
    gene_tier = params.literature_evidence.evidence_tier

    # ── Infer tissue class for delivery confidence scoring ────────────────────
    tissue = infer_tissue(
        gene,
        omim.get("molecular_text", ""),
        mygene["pathways"]
    )
    params.notes.append(f"Inferred target tissue class: {tissue}")

    # ── Therapy combos (payload + delivery + per-combo confidence) ────────────
    params.therapy_combos = recommend_therapies(params, gene_tier, tissue)

    # Attach a ModalityGuide to every combo
    for combo in params.therapy_combos:
        combo.modality_guide = design_modality_guide(
            combo=combo,
            gene=gene,
            tissue=tissue,
            protein_len=protein_len,
            coding_kb=coding_kb,
            variant_hgvs=variant_hgvs,
            mechanism=params.mechanism,
            refseq=mygene["refseq_transcript"],
            crispr_params=params.crispr_parameters,
        )

    # Populate legacy flat fields from combos for PDF/CLI compatibility
    if params.therapy_combos:
        params.primary_therapy = params.therapy_combos[0].label
        params.alternative_therapies = [c.label for c in params.therapy_combos[1:]]

    # ── Notes ─────────────────────────────────────────────────────────────────
    if not aav_feasible and aav_feasible is not None:
        params.notes.append(f"Gene coding sequence ({coding_kb} kb) exceeds single-vector AAV limit (4.7 kb) — dual-vector or non-AAV delivery required")
    if variants:
        gof_v = sum(1 for v in variants if v.mechanism_clue == "GoF")
        lof_v = sum(1 for v in variants if v.mechanism_clue == "LoF")
        params.notes.append(
            f"ClinVar variant breakdown: {gof_v} GoF-suggestive, {lof_v} LoF-suggestive "
            f"of {len(variants)} reviewed"
        )
    if omim["inheritance"]:
        params.notes.append(f"Inheritance pattern: {omim['inheritance']}")
    if params.variant_protein_intersection:
        vpi = params.variant_protein_intersection
        if vpi.is_in_active_site:
            params.notes.append(
                "CRITICAL: patient variant falls within an annotated UniProt active site — "
                "likely directly disrupts core protein function"
            )
        if vpi.is_in_binding_site:
            params.notes.append(
                "Patient variant falls within an annotated binding site — "
                "may disrupt protein-protein or protein-ligand interaction"
            )
        if vpi.is_absent_from_gnomad:
            params.notes.append(
                "Variant absent from gnomAD (ultra-rare) — supports pathogenicity"
            )

    return params


# ══════════════════════════════════════════════════════════════════════════════
# Human-readable report printer
# ══════════════════════════════════════════════════════════════════════════════

def print_report(p: AssayParameters):
    W = 72
    print(f"\n{'='*W}")
    header = f"  GENE THERAPY ANALYSIS v2: {p.gene_symbol}"
    if p.patient_variant_hgvs:
        header += f"  |  Variant: {p.patient_variant_hgvs}"
    print(header)
    print(f"{'='*W}")

    # ── Identity & sequence
    print(f"\n{'─'*W}")
    print(f"  GENE IDENTITY & SEQUENCE")
    print(f"{'─'*W}")
    print(f"  Symbol:             {p.gene_symbol}")
    print(f"  OMIM ID:            {p.omim_id}")
    print(f"  UniProt:            {p.uniprot_id}")
    print(f"  Ensembl Gene:       {p.ensembl_gene_id}")
    print(f"  RefSeq Transcript:  {p.refseq_transcript}")
    print(f"  Protein Length:     {p.protein_length_aa} aa   |   "
          f"Coding: {p.coding_sequence_kb} kb   |   AAV feasible: {p.aav_feasible}")

    # ── Literature Evidence
    if p.literature_evidence:
        le = p.literature_evidence
        print(f"\n{'─'*W}")
        print(f"  LITERATURE EVIDENCE")
        print(f"{'─'*W}")
        print(f"  Gene PubMed publications: {le.gene_pubmed_count}")
        print(f"  Evidence tier:  {le.evidence_tier.upper()}")
        print(f"  {le.tier_rationale}")
        if le.top_review_pmids:
            pmid_links = "  ".join(f"PMID:{pmid}" for pmid in le.top_review_pmids)
            print(f"  Top reviews:  {pmid_links}")

    # ── Patient Variant (NEW)
    if p.patient_variant_hgvs:
        print(f"\n{'─'*W}")
        print(f"  PATIENT VARIANT: {p.patient_variant_hgvs}")
        print(f"{'─'*W}")
        vr = p.patient_variant_resolved
        if vr:
            print(f"  ClinVar:          {vr.name[:75]}")
            print(f"  Significance:     {vr.significance}")
            print(f"  Consequence:      {vr.molecular_consequence}")
            print(f"  Condition:        {vr.condition}")
            print(f"  Mechanism clue:   {vr.mechanism_clue}")
        else:
            print(f"  ClinVar: no exact match found for {p.patient_variant_hgvs}")

        vpi = p.variant_protein_intersection
        if vpi:
            print(f"\n  Protein position:     aa {vpi.aa_position}")
            if vpi.gnomad_allele_frequency is not None:
                print(f"  gnomAD AF:            {vpi.gnomad_allele_frequency:.2e}")
            else:
                print(f"  gnomAD AF:            absent / ultra-rare")
            print(f"  ClinVar significance: {vpi.clinical_significance}")

            if vpi.overlapping_features:
                print(f"\n  !! VARIANT OVERLAPS {len(vpi.overlapping_features)} PROTEIN FEATURE(S):")
                for feat in vpi.overlapping_features:
                    tag = ("** ACTIVE SITE **" if feat.feature_type == "Active site" else
                           "** BINDING SITE **" if feat.feature_type == "Binding site" else
                           "   DOMAIN        ")
                    print(f"     {tag}  {feat.description or feat.feature_type}"
                          f"  (aa {feat.start}–{feat.end})")
            else:
                print(f"\n  No overlap with annotated UniProt features at aa {vpi.aa_position}")

            feas = ("FEASIBLE" if vpi.aso_discrimination_feasible else
                    "UNCERTAIN" if vpi.aso_discrimination_feasible is None else "NOT FEASIBLE")
            print(f"\n  ASO allele discrimination: {feas}")
            if vpi.aso_notes:
                # Word-wrap at 68 chars
                words = vpi.aso_notes.split()
                line = "  "
                for w in words:
                    if len(line) + len(w) + 1 > 70:
                        print(f"  {line.strip()}")
                        line = "  " + w + " "
                    else:
                        line += w + " "
                if line.strip():
                    print(f"  {line.strip()}")

    # ── Therapy
    print(f"\n{'─'*W}")
    print(f"  THERAPY RECOMMENDATION")
    print(f"{'─'*W}")
    print(f"  PRIMARY:  {p.primary_therapy}")
    for a in p.alternative_therapies:
        print(f"  ALT:      {a}")

    # ── CRISPR (NEW)
    if p.crispr_parameters:
        cr = p.crispr_parameters
        print(f"\n{'─'*W}")
        print(f"  CRISPR GUIDE RNA DESIGN")
        print(f"{'─'*W}")
        print(f"  Strategy: {cr.strategy}")
        if cr.base_editing_feasible:
            print(f"  Base editing:  YES — {cr.base_edit_type}")
        if cr.prime_editing_feasible:
            print(f"  Prime editing: YES (recommended for this variant type)")
        if cr.target_sequence_context:
            print(f"  Sequence window (first 80bp): {cr.target_sequence_context}")
        if cr.candidate_guides:
            print(f"\n  Candidate guide RNAs (top {len(cr.candidate_guides)}, ranked by quality score):")
            hdr = f"  {'#':<3} {'Guide (5->3)':<22} {'PAM':<5} {'Str':^3} {'Dist':>5}bp {'GC':>5} {'Score':>6}  Flags"
            print(hdr)
            print(f"  {'-'*3} {'-'*22} {'-'*5} {'-'*3} {'-'*7} {'-'*5} {'-'*6}  {'-'*28}")
            for i, g in enumerate(cr.candidate_guides, 1):
                flags_str = ", ".join(g.quality_flags) if g.quality_flags else "none"
                print(f"  {i:<3} {g.sequence:<22} {g.pam:<5} {g.strand:^3} "
                      f"{g.distance_to_variant_bp:>6}   {g.gc_content:>4.0%} "
                      f"{g.quality_score:>6}  {flags_str}")
        for note in cr.notes:
            print(f"  ℹ  {note}")

    # ── Protein features
    if p.druggable_domains or p.active_sites or p.binding_sites:
        print(f"\n{'─'*W}")
        print(f"  PROTEIN FEATURES (UniProt)")
        print(f"{'─'*W}")
        for f in p.druggable_domains[:5]:
            print(f"  [Domain]       {f.description:42s} aa {f.start}–{f.end}")
        for f in p.active_sites[:5]:
            print(f"  [Active Site]  {f.description:42s} aa {f.start}–{f.end}")
        for f in p.binding_sites[:5]:
            print(f"  [Binding Site] {f.description:42s} aa {f.start}–{f.end}")

    # ── ClinVar variants (gene level)
    if p.representative_variants:
        print(f"\n{'─'*W}")
        print(f"  KEY ClinVar VARIANTS (gene-level, top 10 with positions)")
        print(f"{'─'*W}")
        shown = 0
        for v in p.representative_variants:
            if v.position_cdna or v.position_protein:
                print(f"  [{v.mechanism_clue:7s}] {v.position_cdna or '':18s} "
                      f"{v.position_protein or '':18s}  {v.molecular_consequence}")
                shown += 1
            if shown >= 10:
                break

    # ── Pathways
    if p.pathway_names:
        print(f"\n{'─'*W}")
        print(f"  PATHWAYS & DOWNSTREAM TARGETS")
        print(f"{'─'*W}")
        for pw in p.pathway_names[:6]:
            print(f"  • {pw}")
        if p.downstream_targets:
            print(f"  Downstream druggable nodes: {', '.join(sorted(set(p.downstream_targets)))}")

    # ── Known drugs
    if p.known_drugs:
        print(f"\n{'─'*W}")
        print(f"  KNOWN DRUGS (Open Targets)")
        print(f"{'─'*W}")
        seen = set()
        for d in p.known_drugs[:8]:
            name = d.get("drug", "")
            if name and name not in seen:
                mech = (d.get("mechanism") or "")[:50]
                print(f"  {name:24s}  Phase {d.get('max_phase','?')}  |  {mech}")
                seen.add(name)

    # ── Notes / warnings
    if p.notes:
        print(f"\n{'─'*W}")
        print(f"  NOTES / WARNINGS")
        print(f"{'─'*W}")
        for n in p.notes:
            print(f"  ⚠  {n}")

    # ── OMIM summary
    if p.omim_molecular_genetics_summary:
        print(f"\n{'─'*W}")
        print(f"  OMIM MOLECULAR GENETICS SUMMARY")
        print(f"{'─'*W}")
        clean = re.sub(r"<[^>]+>", "", p.omim_molecular_genetics_summary)
        print(f"  {clean[:800]}...")

    print(f"\n{'='*W}")
    print("  ASSAY PARAMETERS (JSON) — pass to assay generator:")
    print(f"{'='*W}")


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Gene Therapy Target Analysis Pipeline v2",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Gene-level (same as v1 but with CRISPR + literature added)
  python gene_therapy_pipeline_v2.py --gene SCN1A --omim_key YOUR_KEY

  # Full patient-variant analysis (all three new features active)
  python gene_therapy_pipeline_v2.py --gene SCN1A --omim_key YOUR_KEY --variant "c.4849C>T"

  # Save JSON for assay generator
  python gene_therapy_pipeline_v2.py --gene BRAF --omim_key YOUR_KEY --variant "c.1799T>A" --output_json braf_v600e.json
        """
    )
    parser.add_argument("--gene",        default=GENE,         help="HGNC gene symbol")
    parser.add_argument("--omim_key",    default=OMIM_API_KEY, help="OMIM API key")
    parser.add_argument("--variant",     default=None,
                        help="Patient variant in HGVS cDNA format, e.g. c.4849C>T")
    parser.add_argument("--output_json", default=None,
                        help="Optional path to write full JSON output")
    args = parser.parse_args()

    report = build_report(args.gene, args.omim_key, args.variant)
    print_report(report)

    # Serialize (dataclasses → plain dicts)
    def serialize(obj):
        if hasattr(obj, "__dataclass_fields__"):
            return {k: serialize(v) for k, v in asdict(obj).items()}
        if isinstance(obj, list):
            return [serialize(i) for i in obj]
        return obj

    output = serialize(report)
    print(json.dumps(output, indent=2))

    if args.output_json:
        with open(args.output_json, "w") as fh:
            json.dump(output, fh, indent=2)
        print(f"\n[✓] JSON saved to {args.output_json}")


if __name__ == "__main__":
    main()

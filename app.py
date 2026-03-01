"""
TTPicker — Gene Therapy Target Analysis Web App
Flask backend: wraps gene_therapy_pipeline_v2.py
Run: python app.py  →  http://localhost:5000
"""

import json
import io
import re
import traceback
from dataclasses import asdict
from flask import Flask, request, jsonify, render_template, send_file

# Load .env file if present (for ANTHROPIC_API_KEY)
try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass

# Import pipeline (same directory)
from pipeline import build_report, AssayParameters, GeneNotFoundError

# Import assay generator
from assay_generator import generate_protocol, select_templates

app = Flask(__name__)


def serialize(obj):
    """Recursively convert dataclasses → plain dicts for JSON."""
    if hasattr(obj, "__dataclass_fields__"):
        return {k: serialize(v) for k, v in asdict(obj).items()}
    if isinstance(obj, list):
        return [serialize(i) for i in obj]
    return obj


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/analyze", methods=["POST"])
def analyze():
    data    = request.get_json()
    gene    = (data.get("gene") or "").strip().upper()
    variant = (data.get("variant") or "").strip() or None

    if not gene:
        return jsonify({"error": "Gene symbol is required"}), 400

    try:
        from pipeline import OMIM_API_KEY
        report = build_report(gene, OMIM_API_KEY, variant)
        return jsonify({"ok": True, "data": serialize(report)})
    except GeneNotFoundError as e:
        # User-facing error — no stack trace needed
        return jsonify({"error": str(e)}), 400
    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": str(e)}), 500


@app.route("/pdf", methods=["POST"])
def export_pdf():
    """Generate a PDF report from the analysis JSON and return it for download."""
    from reportlab.lib.pagesizes import letter
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import inch
    from reportlab.lib import colors
    from reportlab.platypus import (
        SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
        HRFlowable, KeepTogether
    )

    payload = request.get_json()
    d = payload.get("data", {})

    buf = io.BytesIO()
    doc = SimpleDocTemplate(
        buf, pagesize=letter,
        leftMargin=0.75*inch, rightMargin=0.75*inch,
        topMargin=0.75*inch, bottomMargin=0.75*inch,
        title=f"TTPicker Report — {d.get('gene_symbol','')}"
    )

    SS = getSampleStyleSheet()
    accent = colors.HexColor("#1a7a4a")

    title_style = ParagraphStyle("TTitle", parent=SS["Title"],
                                 fontSize=20, textColor=accent, spaceAfter=4)
    sub_style   = ParagraphStyle("TSub", parent=SS["Normal"],
                                 fontSize=10, textColor=colors.grey, spaceAfter=16)
    h1_style    = ParagraphStyle("TH1", parent=SS["Heading1"],
                                 fontSize=13, textColor=accent,
                                 spaceBefore=14, spaceAfter=4)
    body_style  = ParagraphStyle("TBody", parent=SS["Normal"],
                                 fontSize=9, leading=14, spaceAfter=4)
    warn_style  = ParagraphStyle("TWarn", parent=SS["Normal"],
                                 fontSize=9, leading=13,
                                 textColor=colors.HexColor("#b45309"),
                                 spaceAfter=3)
    code_style  = ParagraphStyle("TCode", parent=SS["Normal"],
                                 fontSize=8, fontName="Courier",
                                 backColor=colors.HexColor("#f3f4f6"),
                                 leading=12, spaceAfter=2)

    def h1(text):
        return [HRFlowable(width="100%", thickness=1,
                           color=accent, spaceAfter=2),
                Paragraph(text, h1_style)]

    def body(text):
        clean = re.sub(r"<[^>]+>", "", str(text))
        return Paragraph(clean, body_style)

    def kv(label, value):
        if not value:
            return []
        return [body(f"<b>{label}:</b>  {value}")]

    story = []

    # ── Title
    gene  = d.get("gene_symbol", "")
    var   = d.get("patient_variant_hgvs") or ""
    story.append(Paragraph(f"TTPicker Analysis Report", title_style))
    story.append(Paragraph(
        f"Gene: {gene}" + (f"  |  Variant: {var}" if var else ""),
        sub_style
    ))
    story.append(Spacer(1, 6))

    # ── Gene identity
    story += h1("Gene Identity & Sequence")
    for label, key in [
        ("OMIM ID",           "omim_id"),
        ("UniProt",           "uniprot_id"),
        ("Ensembl Gene ID",   "ensembl_gene_id"),
        ("RefSeq Transcript", "refseq_transcript"),
        ("Protein Length",    "protein_length_aa"),
        ("Coding Sequence",   "coding_sequence_kb"),
        ("AAV Feasible",      "aav_feasible"),
    ]:
        val = d.get(key)
        if val is not None:
            suffix = " aa" if key == "protein_length_aa" else " kb" if key == "coding_sequence_kb" else ""
            story += kv(label, f"{val}{suffix}")

    # ── Literature evidence
    le = d.get("literature_evidence") or {}
    if le:
        story += h1("Literature Evidence")
        story += kv("Gene PubMed publications", le.get("gene_pubmed_count"))
        story += kv("Evidence Tier", (le.get("evidence_tier") or "").upper())
        if le.get("tier_rationale"):
            story.append(body(le["tier_rationale"]))
        pmids = le.get("top_review_pmids") or []
        if pmids:
            story += kv("Top Review PMIDs",
                        "  ".join(f"PMID:{p}" for p in pmids))

    # ── Patient variant
    vpi_data = d.get("variant_protein_intersection") or {}
    vr_data  = d.get("patient_variant_resolved") or {}
    if var:
        story += h1(f"Patient Variant Analysis: {var}")
        if vr_data.get("name"):
            story += kv("ClinVar Record", vr_data["name"][:90])
            story += kv("Significance",   vr_data.get("significance"))
            story += kv("Consequence",    vr_data.get("molecular_consequence"))
            story += kv("Condition",      vr_data.get("condition"))
        else:
            story.append(body("No exact ClinVar match found — positional analysis used."))

        if vpi_data:
            story += kv("Protein Position (AA)", vpi_data.get("aa_position"))
            af = vpi_data.get("gnomad_allele_frequency")
            story += kv("gnomAD Allele Freq",
                        f"{af:.2e}" if af else "Absent / ultra-rare")
            story += kv("ClinVar Significance", vpi_data.get("clinical_significance"))

            feats = vpi_data.get("overlapping_features") or []
            if feats:
                story.append(body(f"<b>!! Variant overlaps {len(feats)} protein feature(s):</b>"))
                for f in feats:
                    tag = ("ACTIVE SITE" if f.get("feature_type") == "Active site"
                           else "BINDING SITE" if f.get("feature_type") == "Binding site"
                           else "DOMAIN")
                    story.append(body(f"  [{tag}]  {f.get('description') or f.get('feature_type')}  "
                                      f"(aa {f.get('start')}–{f.get('end')})"))
            else:
                story.append(body("No overlap with annotated UniProt protein features."))

            aso_feas = vpi_data.get("aso_discrimination_feasible")
            story += kv("ASO Discrimination",
                        "Feasible" if aso_feas else "Uncertain" if aso_feas is None else "Not Feasible")
            if vpi_data.get("aso_notes"):
                story.append(body(vpi_data["aso_notes"]))

    # ── Therapy combos
    combos = d.get("therapy_combos") or []
    if combos:
        story += h1("Therapy Recommendations")
        tdata = [["#", "Payload", "Delivery Vehicle", "Confidence", "Notes"]]
        for c in combos:
            rank_label = "★ PRIMARY" if c.get("rank") == 1 else f"#{c.get('rank','')}"
            tdata.append([
                rank_label,
                c.get("payload", ""),
                c.get("delivery", ""),
                (c.get("confidence") or "").upper(),
                (c.get("notes") or "")[:80],
            ])
        conf_colors = {
            "HIGH":         colors.HexColor("#dcfce7"),
            "MEDIUM":       colors.HexColor("#fef9c3"),
            "LOW":          colors.HexColor("#fee2e2"),
            "EXPERIMENTAL": colors.HexColor("#f3f4f6"),
        }
        tbl = Table(tdata, colWidths=[
            0.7*inch, 1.8*inch, 1.7*inch, 0.75*inch, 2.55*inch
        ])
        style_cmds = [
            ("BACKGROUND",   (0,0), (-1,0), accent),
            ("TEXTCOLOR",    (0,0), (-1,0), colors.white),
            ("FONTNAME",     (0,0), (-1,0), "Helvetica-Bold"),
            ("FONTSIZE",     (0,0), (-1,-1), 8),
            ("GRID",         (0,0), (-1,-1), 0.4, colors.HexColor("#e5e7eb")),
            ("TOPPADDING",   (0,0), (-1,-1), 4),
            ("BOTTOMPADDING",(0,0), (-1,-1), 4),
            # Primary row highlight
            ("BACKGROUND",   (0,1), (-1,1), colors.HexColor("#e8f5ee")),
            ("FONTNAME",     (0,1), (1,1),  "Helvetica-Bold"),
        ]
        # Per-confidence-cell colouring
        for i, c in enumerate(combos, 1):
            conf_key = (c.get("confidence") or "").upper()
            bg = conf_colors.get(conf_key, colors.white)
            style_cmds.append(("BACKGROUND", (3, i), (3, i), bg))
        tbl.setStyle(TableStyle(style_cmds))
        story.append(tbl)
        story.append(Spacer(1, 6))
        # Rationale footnotes
        for c in combos:
            if c.get("confidence_rationale"):
                story.append(body(f"#{c.get('rank')} rationale: {c['confidence_rationale']}"))
        story.append(Spacer(1, 4))

    # ── Modality guides (one per therapy combo)
    modality_titles = {
        "aav":       "AAV Design Parameters",
        "mrna_lnp":  "mRNA + LNP Design Parameters",
        "aso_sirna": "ASO / siRNA Design Parameters",
        "crispr":    "CRISPR Design Parameters",
        "cast":      "CAST / Transposon Design Parameters",
    }
    for c in combos:
        mg = c.get("modality_guide") or {}
        if not mg:
            continue
        mod_type  = mg.get("modality_type", "other")
        sec_title = modality_titles.get(mod_type, "Modality Design Parameters")
        rank_label = "★ PRIMARY — " if c.get("rank") == 1 else f"Alt #{c.get('rank','')}: "
        story += h1(f"{rank_label}{sec_title}")
        story.append(body(f"<b>Payload:</b> {c.get('payload','')}  |  <b>Delivery:</b> {c.get('delivery','')}"))
        story.append(Spacer(1, 4))

        # Key parameters as a 2-col table
        kps = mg.get("key_parameters") or []
        if kps:
            kp_data = []
            for kp in kps:
                colon = kp.find(":")
                k = kp[:colon].strip() if colon > -1 else "Parameter"
                v = kp[colon+1:].strip() if colon > -1 else kp
                kp_data.append([k, v])
            kp_tbl = Table(kp_data, colWidths=[1.8*inch, 5.7*inch])
            kp_tbl.setStyle(TableStyle([
                ("FONTNAME",    (0,0), (0,-1), "Helvetica-Bold"),
                ("FONTSIZE",    (0,0), (-1,-1), 8),
                ("TEXTCOLOR",   (0,0), (0,-1), accent),
                ("ROWBACKGROUNDS",(0,0),(-1,-1),
                 [colors.white, colors.HexColor("#f9fafb")]),
                ("GRID",        (0,0), (-1,-1), 0.3, colors.HexColor("#e5e7eb")),
                ("TOPPADDING",  (0,0), (-1,-1), 4),
                ("BOTTOMPADDING",(0,0),(-1,-1), 4),
                ("LEFTPADDING", (0,0), (-1,-1), 6),
            ]))
            story.append(kp_tbl)
            story.append(Spacer(1, 6))

        # Modality-specific narrative fields
        field_labels = {
            "aav_itr_config":           "Vector Configuration",
            "aav_tropism_note":         "Tropism / Promoter Note",
            "aav_packaging_note":       "Manufacture",
            "mrna_utr_design":          "UTR Design",
            "mrna_codon_opt":           "Codon Optimisation",
            "lnp_targeting":            "LNP Targeting",
            "mrna_stability":           "Stability & Redosing",
            "aso_mechanism":            "Mechanism",
            "aso_target_region":        "Target Region",
            "aso_design_window":        "Design Window",
            "sirna_strand_bias":        "Strand Bias",
            "crispr_strategy":          "Strategy",
            "crispr_window":            "Editing Window",
            "crispr_delivery_form":     "Delivery Form",
            "crispr_offtarget":         "Off-target Risk",
            "cast_system":              "System",
            "cast_target_site":         "Target Site",
            "cast_transposase_delivery":"Transposase Delivery",
            "cast_offtarget_risk":      "Integration Risk",
        }
        for field, label in field_labels.items():
            val = mg.get(field)
            if val:
                story.append(body(f"<b>{label}:</b>  {val}"))

        # CRISPR guide table (if CRISPR modality)
        if mod_type == "crispr":
            cr = d.get("crispr_parameters") or {}
            if cr.get("target_sequence_context"):
                story.append(body("<b>Sequence window:</b>"))
                story.append(Paragraph(cr["target_sequence_context"], code_style))
            guides = cr.get("candidate_guides") or []
            if guides:
                gdata = [["#", "Guide (5'→3')", "PAM", "Str", "Dist(bp)", "GC", "Score", "Flags"]]
                for i, g in enumerate(guides, 1):
                    flags = ", ".join(g.get("quality_flags") or []) or "none"
                    gdata.append([str(i), g.get("sequence",""), g.get("pam",""),
                                  g.get("strand",""), str(g.get("distance_to_variant_bp","")),
                                  f"{g.get('gc_content',0):.0%}", str(g.get("quality_score","")), flags])
                gtbl = Table(gdata, colWidths=[0.25*inch,1.7*inch,0.45*inch,0.3*inch,0.55*inch,0.35*inch,0.45*inch,2.5*inch])
                gtbl.setStyle(TableStyle([
                    ("BACKGROUND",(0,0),(-1,0),accent),("TEXTCOLOR",(0,0),(-1,0),colors.white),
                    ("FONTNAME",(0,0),(-1,0),"Helvetica-Bold"),("FONTSIZE",(0,0),(-1,-1),7),
                    ("FONTNAME",(0,1),(1,-1),"Courier"),
                    ("ROWBACKGROUNDS",(0,1),(-1,-1),[colors.white,colors.HexColor("#f9fafb")]),
                    ("GRID",(0,0),(-1,-1),0.4,colors.HexColor("#e5e7eb")),
                    ("TOPPADDING",(0,0),(-1,-1),3),("BOTTOMPADDING",(0,0),(-1,-1),3),
                ]))
                story.append(Spacer(1,4))
                story.append(gtbl)

        # Design notes
        for note in (mg.get("design_notes") or []):
            story.append(body(f"ℹ  {note}"))
        story.append(Spacer(1, 8))

    # ── Protein features
    domains  = d.get("druggable_domains") or []
    asites   = d.get("active_sites") or []
    bsites   = d.get("binding_sites") or []
    if domains or asites or bsites:
        story += h1("Protein Features (UniProt)")
        for ftype, items in [("Domain", domains[:5]),
                              ("Active Site", asites[:5]),
                              ("Binding Site", bsites[:5])]:
            for f in items:
                story.append(body(
                    f"<b>[{ftype}]</b>  {f.get('description') or '—'}  "
                    f"aa {f.get('start')}–{f.get('end')}"
                ))

    # ── Key ClinVar variants
    cvars = [v for v in (d.get("representative_variants") or [])
             if v.get("position_cdna") or v.get("position_protein")][:10]
    if cvars:
        story += h1("Key ClinVar Variants")
        vdata = [["Mechanism", "cDNA", "Protein", "Consequence"]]
        for v in cvars:
            vdata.append([
                v.get("mechanism_clue",""),
                v.get("position_cdna",""),
                v.get("position_protein",""),
                v.get("molecular_consequence",""),
            ])
        vtbl = Table(vdata, colWidths=[0.8*inch, 1.3*inch, 1.3*inch, 3.1*inch])
        vtbl.setStyle(TableStyle([
            ("BACKGROUND",  (0,0), (-1,0), accent),
            ("TEXTCOLOR",   (0,0), (-1,0), colors.white),
            ("FONTNAME",    (0,0), (-1,0), "Helvetica-Bold"),
            ("FONTSIZE",    (0,0), (-1,-1), 8),
            ("ROWBACKGROUNDS",(0,1),(-1,-1),
             [colors.white, colors.HexColor("#f9fafb")]),
            ("GRID", (0,0), (-1,-1), 0.4, colors.HexColor("#e5e7eb")),
            ("TOPPADDING",  (0,0), (-1,-1), 3),
            ("BOTTOMPADDING",(0,0),(-1,-1), 3),
        ]))
        story.append(vtbl)

    # ── Pathways
    pathways = d.get("pathway_names") or []
    if pathways:
        story += h1("Pathways & Downstream Targets")
        for pw in pathways[:6]:
            story.append(body(f"• {pw}"))
        dnodes = d.get("downstream_targets") or []
        if dnodes:
            story += kv("Downstream druggable nodes", ", ".join(sorted(set(dnodes))))

    # ── Known drugs
    drugs = d.get("known_drugs") or []
    if drugs:
        story += h1("Known Drugs (Open Targets)")
        ddata = [["Drug", "Phase", "Mechanism"]]
        seen = set()
        for dr in drugs[:10]:
            name = dr.get("drug","")
            if name and name not in seen:
                ddata.append([
                    name,
                    str(dr.get("max_phase","?")),
                    (dr.get("mechanism") or "")[:60],
                ])
                seen.add(name)
        dtbl = Table(ddata, colWidths=[1.5*inch, 0.5*inch, 4.5*inch])
        dtbl.setStyle(TableStyle([
            ("BACKGROUND",  (0,0), (-1,0), accent),
            ("TEXTCOLOR",   (0,0), (-1,0), colors.white),
            ("FONTNAME",    (0,0), (-1,0), "Helvetica-Bold"),
            ("FONTSIZE",    (0,0), (-1,-1), 8),
            ("ROWBACKGROUNDS",(0,1),(-1,-1),
             [colors.white, colors.HexColor("#f9fafb")]),
            ("GRID", (0,0), (-1,-1), 0.4, colors.HexColor("#e5e7eb")),
            ("TOPPADDING",  (0,0), (-1,-1), 3),
            ("BOTTOMPADDING",(0,0),(-1,-1), 3),
        ]))
        story.append(dtbl)

    # ── Notes
    notes = d.get("notes") or []
    if notes:
        story += h1("Clinical Notes & Warnings")
        for n in notes:
            story.append(Paragraph(f"⚠  {n}", warn_style))

    # ── OMIM summary
    omim_text = d.get("omim_molecular_genetics_summary") or ""
    if omim_text:
        story += h1("OMIM Molecular Genetics Summary")
        clean = re.sub(r"<[^>]+>", "", omim_text)[:1200]
        story.append(body(clean + "..."))

    doc.build(story)
    buf.seek(0)
    fname = f"TTPicker_{gene}" + (f"_{var.replace('>','_')}" if var else "") + ".pdf"
    return send_file(buf, mimetype="application/pdf",
                     as_attachment=True, download_name=fname)


# ── Assay Generator routes ─────────────────────────────────────────────────────

@app.route("/assay", methods=["POST"])
def generate_assay():
    """
    Accepts pipeline data + user answers from the assay question form.
    Returns JSON: { protocol_html, protocol_text, ot2_script, templates_used }
    """
    try:
        data = request.get_json()
        if not data:
            return jsonify({"error": "No JSON body provided"}), 400
        pipeline_data = data.get("pipeline_data")
        user_answers  = data.get("user_answers")
        if not pipeline_data:
            return jsonify({"error": "Missing pipeline_data"}), 400
        if not user_answers:
            return jsonify({"error": "Missing user_answers"}), 400
        result = generate_protocol(pipeline_data, user_answers)
        return jsonify(result)
    except ValueError as e:
        return jsonify({"error": str(e)}), 400
    except Exception as e:
        traceback.print_exc()
        return jsonify({"error": f"Assay generation failed: {str(e)}"}), 500


@app.route("/ot2", methods=["POST"])
def get_ot2_script():
    """
    Convenience pass-through. In practice /assay already returns ot2_script;
    this endpoint exists for direct OT-2 re-fetches if needed.
    """
    try:
        data = request.get_json()
        ot2_script = (data or {}).get("ot2_script", "")
        if not ot2_script:
            return jsonify({"error": "No OT-2 script provided"}), 400
        return jsonify({"ot2_script": ot2_script})
    except Exception as e:
        return jsonify({"error": str(e)}), 500


if __name__ == "__main__":
    print("\n  TTPicker running at http://localhost:5000\n")
    app.run(debug=False, port=5000)

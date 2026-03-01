"""
TTPicker — Assay Generator Module
==================================
Step 2 of the TTPicker pipeline: takes pipeline JSON output + user experimental
parameters and produces a human-readable validation protocol document plus an
Opentrons OT-2 liquid-handling script.

Usage:
    from assay_generator import generate_protocol, select_templates, build_claude_context

    templates = select_templates(pipeline_data)
    context = build_claude_context(pipeline_data, user_answers)
    result = generate_protocol(pipeline_data, user_answers)
    # result = { "protocol_html": str, "protocol_text": str, "ot2_script": str }
"""

import json
import os
import re
from datetime import datetime
from typing import Any

# ---------------------------------------------------------------------------
# Protocol Templates
# ---------------------------------------------------------------------------

WESTERN_BLOT_TEMPLATE = {
    "assay_type": "Western Blot - Protein Expression Validation",
    "objective": "{{OBJECTIVE}}",
    "rationale": "{{RATIONALE}}",
    "sections": [
        {
            "title": "Materials",
            "items": [
                "Primary antibody: {{PRIMARY_ANTIBODY}}",
                "Secondary antibody: HRP-conjugated anti-[species] IgG",
                "Loading control antibody: Anti-beta-actin (42 kDa) or Anti-GAPDH",
                "Lysis buffer: RIPA buffer + protease inhibitor cocktail",
                "BCA protein quantification kit",
                "SDS-PAGE gel: {{GEL_PERCENTAGE}}%",
                "PVDF membrane (0.45 um pore for >100 kDa proteins)",
                "ECL chemiluminescent substrate",
                "{{TREATMENT_REAGENT}}",
                "Vehicle control: {{VEHICLE_CONTROL}}",
                "Positive control: {{POSITIVE_CONTROL}}",
            ],
        },
        {
            "title": "Cell Treatment",
            "steps": [
                "Seed {{CELL_MODEL}} at {{SEEDING_DENSITY}} per well in {{PLATE_FORMAT}}.",
                "Allow cells to adhere/differentiate for {{PRE_TREATMENT_PERIOD}}.",
                "Treat with {{TREATMENT_REAGENT}} at doses: {{DOSES}}.",
                "Include: untreated negative control, vehicle control ({{VEHICLE_CONTROL}}), positive control ({{POSITIVE_CONTROL}}).",
                "Maintain {{BIO_REPLICATES}} biological replicates per condition.",
                "Harvest at timepoints: {{TIMEPOINTS}}.",
            ],
        },
        {
            "title": "Protein Extraction and Quantification",
            "steps": [
                "Aspirate media. Wash cells 2x with ice-cold PBS.",
                "Add 100-200 uL RIPA buffer + 1x protease inhibitor per well.",
                "Incubate on ice for 20 min with occasional vortexing.",
                "Centrifuge at 14,000 x g for 15 min at 4C. Collect supernatant.",
                "Quantify protein concentration by BCA assay (duplicate reads).",
                "Normalise all samples to {{PROTEIN_LOAD}} ug total protein.",
            ],
        },
        {
            "title": "SDS-PAGE and Transfer",
            "steps": [
                "Prepare samples in 4x Laemmli buffer. Boil at 95C for 5 min.",
                "Load {{PROTEIN_LOAD}} ug per lane on {{GEL_PERCENTAGE}}% SDS-PAGE gel.",
                "Run at 100V for ~90 min (until dye front reaches gel bottom).",
                "Transfer to PVDF membrane at 100V for 60-90 min (wet transfer) or 25V overnight (semi-dry).",
                "Block membrane in 5% non-fat milk in TBST for 1 hour at RT.",
            ],
        },
        {
            "title": "Immunodetection",
            "steps": [
                "Incubate with primary antibody ({{PRIMARY_ANTIBODY}}, {{PRIMARY_DILUTION}}) overnight at 4C.",
                "Wash 3x in TBST (5 min each).",
                "Incubate with HRP-conjugated secondary antibody (1:5000-1:10000) for 1 hour at RT.",
                "Wash 3x in TBST (5 min each).",
                "Detect with ECL substrate. Image using chemiluminescence imager.",
                "Re-probe with anti-beta-actin or anti-GAPDH for loading control.",
            ],
        },
        {
            "title": "Expected Band Sizes",
            "content": "{{EXPECTED_BANDS}}",
        },
        {
            "title": "Quantification and Statistics",
            "steps": [
                "Densitometry: normalise target band to loading control band per lane.",
                "Express as fold-change relative to untreated negative control.",
                "Statistics: one-way ANOVA with Tukey post-hoc test (n={{BIO_REPLICATES}} biological replicates x {{TECH_REPLICATES}} technical).",
                "Significance threshold: p < 0.05.",
                "{{INTERPRETATION_CRITERIA}}",
            ],
        },
        {
            "title": "Controls Layout",
            "content": "{{CONTROLS_TABLE}}",
        },
    ],
}

RT_PCR_TEMPLATE = {
    "assay_type": "RT-qPCR / RT-PCR - Transcript Level Validation",
    "objective": "{{OBJECTIVE}}",
    "rationale": "{{RATIONALE}}",
    "sections": [
        {
            "title": "Materials",
            "items": [
                "RNA extraction kit (e.g. RNeasy Mini Kit, Qiagen)",
                "DNase I treatment (to remove genomic DNA contamination)",
                "Reverse transcription kit (e.g. SuperScript IV or iScript)",
                "qPCR master mix (e.g. SYBR Green or TaqMan)",
                "Primers: {{FORWARD_PRIMER_DESIGN}}",
                "Primers: {{REVERSE_PRIMER_DESIGN}}",
                "Housekeeping gene primers: GAPDH or ACTB (provide sequences)",
                "{{SPLICE_JUNCTION_PRIMERS}}",
                "Nanodrop or Bioanalyser for RNA QC",
            ],
        },
        {
            "title": "RNA Extraction",
            "steps": [
                "Harvest cells at timepoints: {{TIMEPOINTS}}.",
                "Extract total RNA per manufacturer protocol.",
                "Assess RNA integrity: A260/A280 >= 1.8; RIN >= 7 preferred.",
                "DNase I treat to remove genomic DNA: {{DNASE_PROTOCOL}}.",
                "Normalise RNA to {{RNA_INPUT}} ng per RT reaction.",
            ],
        },
        {
            "title": "Reverse Transcription",
            "steps": [
                "Set up RT reaction with {{RNA_INPUT}} ng total RNA.",
                "Use oligo-dT or random hexamer priming (oligo-dT preferred for mRNA quantification).",
                "RT conditions: 50C x 10 min, 85C x 5 min (SuperScript IV protocol).",
                "Dilute cDNA 1:5 in nuclease-free water before qPCR.",
            ],
        },
        {
            "title": "qPCR Conditions",
            "steps": [
                "95C 3 min initial denaturation.",
                "40 cycles: 95C 15 sec, 60C 30 sec, 72C 30 sec.",
                "Melt curve: 65C to 95C (0.5C increments) to confirm single amplicon.",
                "Run {{TECH_REPLICATES}} technical replicates per biological replicate.",
                "Include no-template control (NTC) and no-RT control per plate.",
            ],
        },
        {
            "title": "Primer Design Guidance",
            "content": "{{PRIMER_DESIGN_DETAIL}}",
        },
        {
            "title": "Quantification and Statistics",
            "steps": [
                "Relative quantification: delta-delta-Ct method normalised to housekeeping gene.",
                "Accept Ct values between 15-35 (outside = flag for dilution or re-extraction).",
                "Statistics: unpaired t-test (2 groups) or one-way ANOVA (3+ groups).",
                "{{INTERPRETATION_CRITERIA}}",
            ],
        },
    ],
}

AMPLICON_SEQ_TEMPLATE = {
    "assay_type": "Amplicon Sequencing - CRISPR Editing Efficiency",
    "objective": "{{OBJECTIVE}}",
    "rationale": "{{RATIONALE}}",
    "sections": [
        {
            "title": "Materials",
            "items": [
                "Genomic DNA extraction kit (e.g. DNeasy Blood and Tissue Kit)",
                "PCR primers flanking cut site: {{AMPLICON_PRIMERS}}",
                "High-fidelity PCR polymerase (e.g. Q5 or Phusion)",
                "Gel electrophoresis: 2% agarose",
                "Amplicon sequencing service (Sanger for initial screen; NGS for quantification)",
                "TIDE or ICE (Inference of CRISPR Edits) analysis software (free online)",
                "{{CRISPR_REAGENTS}}",
            ],
        },
        {
            "title": "Guide RNA and Editor Delivery",
            "content": "{{CRISPR_DELIVERY_DETAIL}}",
        },
        {
            "title": "Amplicon PCR",
            "steps": [
                "Extract gDNA 48-72 hours post-treatment (for NHEJ) or 5-7 days (for HDR/base edit).",
                "Design amplicon primers {{AMPLICON_SIZE}} bp, flanking cut site by >=100 bp each side.",
                "PCR: 95C 30s | 35 cycles [95C 10s, 60C 30s, 72C {{EXTENSION_TIME}}s] | 72C 2min.",
                "Verify single band by 2% agarose gel.",
                "Purify PCR product (column cleanup or gel extraction).",
            ],
        },
        {
            "title": "Sequencing and Analysis",
            "steps": [
                "Sanger sequencing (initial screen): submit 50-100 ng purified amplicon per reaction.",
                "Use TIDE or ICE software to decompose Sanger trace and estimate indel frequency.",
                "NGS (quantitative): MiSeq 2x150 or 2x250 paired-end; >=10,000 reads per sample.",
                "Analysis: CRISPResso2 (free, command-line) or online CRISPResso.",
                "{{BASE_EDIT_ANALYSIS}}",
            ],
        },
        {
            "title": "Expected Outcomes",
            "content": "{{EXPECTED_OUTCOMES}}",
        },
        {
            "title": "Quantification and Statistics",
            "steps": [
                "Report: % alleles edited (indel frequency or base conversion frequency).",
                "Target editing efficiency for therapeutic relevance: {{EFFICIENCY_TARGET}}%.",
                "Statistics: mean +/- SD across biological replicates; n={{BIO_REPLICATES}}.",
            ],
        },
    ],
}

ELISA_TEMPLATE = {
    "assay_type": "ELISA - Secreted / Circulating Protein Quantification",
    "objective": "{{OBJECTIVE}}",
    "rationale": "{{RATIONALE}}",
    "sections": [
        {
            "title": "Materials",
            "items": [
                "ELISA kit: {{ELISA_KIT_SUGGESTION}}",
                "96-well flat-bottom plate (if manual ELISA)",
                "Plate reader (450 nm absorbance)",
                "Cell culture supernatant or plasma/serum (as appropriate)",
                "Standard curve: recombinant {{GENE}} protein or kit standards",
                "Assay diluent and wash buffer per kit",
            ],
        },
        {
            "title": "Sample Collection",
            "steps": [
                "Collect conditioned media at timepoints: {{TIMEPOINTS}}.",
                "For secreted proteins: spin at 300 x g to remove cells. Aliquot supernatant.",
                "For plasma/serum (in vivo validation): collect per institutional protocol.",
                "Store at -80C until assay. Avoid repeated freeze-thaw.",
                "Dilute samples to fall within linear range of standard curve.",
            ],
        },
        {
            "title": "ELISA Procedure",
            "steps": [
                "Follow manufacturer protocol for selected kit.",
                "Run {{TECH_REPLICATES}} technical replicates per sample.",
                "Include: blank (assay diluent only), standard curve (8-point, 2-fold dilution series), positive control, negative control.",
                "Read absorbance at 450 nm (reference: 570 nm) within 30 min of stop solution.",
            ],
        },
        {
            "title": "Expected Protein Levels",
            "content": "{{EXPECTED_LEVELS}}",
        },
        {
            "title": "Quantification and Statistics",
            "steps": [
                "Interpolate concentrations from standard curve (4-parameter logistic fit).",
                "Normalise to cell number or total protein if measuring intracellular.",
                "Statistics: one-way ANOVA with Tukey post-hoc; n={{BIO_REPLICATES}}.",
                "{{INTERPRETATION_CRITERIA}}",
            ],
        },
    ],
}

TEMPLATE_MAP = {
    "western_blot": WESTERN_BLOT_TEMPLATE,
    "rt_pcr": RT_PCR_TEMPLATE,
    "amplicon_seq": AMPLICON_SEQ_TEMPLATE,
    "elisa": ELISA_TEMPLATE,
}


# ---------------------------------------------------------------------------
# Context Builder
# ---------------------------------------------------------------------------

def build_claude_context(pipeline_data: dict, user_answers: dict) -> dict:
    """
    Build the minimal context object sent to Claude for protocol generation.
    Keeps total tokens under ~3000 for the context payload.
    """
    primary = (pipeline_data.get("therapy_combos") or [{}])[0]
    mg = primary.get("modality_guide") or {}
    mod_type = mg.get("modality_type", "")
    cr = pipeline_data.get("crispr_parameters") or {}
    vpi = pipeline_data.get("variant_protein_intersection") or {}
    vr = pipeline_data.get("patient_variant_resolved") or {}

    modality_fields: dict = {}
    if mod_type == "crispr":
        modality_fields = {
            "crispr_strategy": mg.get("crispr_strategy"),
            "crispr_editor": mg.get("crispr_editor"),
            "crispr_window": mg.get("crispr_window"),
            "best_guide": (cr.get("candidate_guides") or [{}])[0],
            "base_edit_type": cr.get("base_edit_type"),
            "sequence_context_snippet": (cr.get("target_sequence_context") or "")[:80],
        }
    elif mod_type == "aav":
        modality_fields = {
            "aav_serotype": mg.get("aav_serotype"),
            "aav_promoter": mg.get("aav_promoter"),
            "aav_itr_config": mg.get("aav_itr_config"),
        }
    elif mod_type == "mrna_lnp":
        modality_fields = {
            "lnp_formulation": mg.get("lnp_formulation"),
            "mrna_modification": mg.get("mrna_modification"),
        }
    elif mod_type == "aso_sirna":
        modality_fields = {
            "aso_mechanism": mg.get("aso_mechanism"),
            "aso_target_region": mg.get("aso_target_region"),
            "aso_chemistry": mg.get("aso_chemistry"),
            "aso_conjugate": mg.get("aso_conjugate"),
        }
    elif mod_type == "cast":
        modality_fields = {
            "cast_system": mg.get("cast_system"),
            "cast_target_site": mg.get("cast_target_site"),
        }

    return {
        "gene": pipeline_data.get("gene_symbol"),
        "omim_id": pipeline_data.get("omim_id"),
        "protein_length_aa": pipeline_data.get("protein_length_aa"),
        "mechanism": pipeline_data.get("mechanism"),
        "inheritance": pipeline_data.get("inheritance"),
        "variant_hgvs": pipeline_data.get("patient_variant_hgvs"),
        "variant_consequence": vr.get("molecular_consequence"),
        "variant_protein": vr.get("position_protein"),
        "variant_aa_position": vpi.get("aa_position"),
        "variant_in_active_site": vpi.get("is_in_active_site"),
        "variant_in_domain": vpi.get("is_in_domain"),
        "overlapping_domain": (vpi.get("overlapping_features") or [{}])[0].get("description"),
        "primary_payload": primary.get("payload"),
        "primary_delivery": primary.get("delivery"),
        "modality_type": mod_type,
        "modality_fields": modality_fields,
        "evidence_tier": (pipeline_data.get("literature_evidence") or {}).get("evidence_tier"),
        "cell_model": user_answers.get("cell_model"),
        "treatment_dose": user_answers.get("treatment_dose"),
        "treatment_timepoints": user_answers.get("treatment_timepoints"),
        "vehicle_control": user_answers.get("vehicle_control"),
        "readout_equipment": user_answers.get("readout_equipment"),
        "positive_control": user_answers.get("positive_control"),
        "bio_replicates": user_answers.get("bio_replicates", 3),
        "tech_replicates": user_answers.get("tech_replicates", 3),
        "cell_passage": user_answers.get("cell_passage"),
    }


# ---------------------------------------------------------------------------
# Template Selection
# ---------------------------------------------------------------------------

def _is_secreted_protein(pipeline_data: dict) -> bool:
    """Heuristic check for secreted/circulating protein."""
    keywords = {"signal", "secreted", "extracellular", "serum", "plasma", "circulating"}
    for d in pipeline_data.get("druggable_domains") or []:
        desc = (d.get("description") or "").lower()
        if any(k in desc for k in keywords):
            return True
    omim_text = (pipeline_data.get("omim_molecular_genetics_summary") or "").lower()
    if any(k in omim_text for k in keywords):
        return True
    return False


def select_templates(pipeline_data: dict) -> list:
    """
    Deterministic template dispatch based on primary payload and mechanism.
    Returns list of template keys from TEMPLATE_MAP.
    """
    primary = (pipeline_data.get("therapy_combos") or [{}])[0]
    payload = (primary.get("payload") or "").lower()
    mechanism = (pipeline_data.get("mechanism") or "").lower()
    templates = []

    if "crispr" in payload:
        templates.append("amplicon_seq")
        templates.append("western_blot")
    elif "aso" in payload and "splice" in payload:
        templates.append("rt_pcr")
    elif "aso" in payload and "knockdown" in payload:
        templates.append("western_blot")
        templates.append("rt_pcr")
    elif "sirna" in payload:
        templates.append("western_blot")
        templates.append("rt_pcr")
    elif "aav" in payload or "mrna" in payload:
        templates.append("western_blot")
    else:
        templates.append("western_blot")

    if "loss of function" in mechanism and _is_secreted_protein(pipeline_data):
        if "elisa" not in templates:
            templates.append("elisa")

    return templates


def select_templates_with_equipment(pipeline_data: dict, equipment: list = None) -> list:
    """
    Extended dispatch that also considers user equipment to add optional templates.
    """
    templates = select_templates(pipeline_data)
    equipment = equipment or []
    eq_lower = [e.lower() for e in equipment]

    if "rt_pcr" not in templates:
        if any("qpcr" in e or "pcr" in e for e in eq_lower):
            primary = (pipeline_data.get("therapy_combos") or [{}])[0]
            payload = (primary.get("payload") or "").lower()
            if "aav" in payload or "mrna" in payload or "aso" in payload:
                templates.append("rt_pcr")

    if "elisa" not in templates:
        if any("plate reader" in e for e in eq_lower):
            if _is_secreted_protein(pipeline_data):
                templates.append("elisa")

    return templates


# ---------------------------------------------------------------------------
# Claude API Prompts
# ---------------------------------------------------------------------------

SYSTEM_PROMPT = """You are an expert molecular biologist specialising in gene therapy
validation assays for rare monogenic diseases. You will receive a structured context
object describing a target gene, therapy modality, and experimental parameters.
Your job is to fill in the {{PLACEHOLDER}} fields in a protocol template to produce
a complete, accurate, human-readable laboratory protocol.

Be specific and practical. Do not hedge excessively. If you suggest an antibody,
name a real commercial one (e.g. "Abcam ab65176 anti-Nav1.1, 1:500"). If you suggest
primer design guidance, give actual rules (e.g. "design 150-200 bp amplicon flanking
the cut site at position c.4849; forward primer: 5\'->3\' upstream of c.4730, reverse
primer: 3\'->5\' downstream of c.4950"). If you do not know a specific catalog number,
say so clearly rather than inventing one.

Return ONLY the filled protocol as a JSON array matching the template structure
provided. Do not include preamble, markdown fences, or explanation outside the JSON.
All {{PLACEHOLDER}} markers must be replaced with actual values - none should remain."""

OT2_SYSTEM_PROMPT = """You are an expert Opentrons OT-2 protocol developer. You will
receive an assay protocol context and must generate a Python script using the
opentrons library (API level 2.15) that handles the liquid-handling steps of the assay.

The script is a TEMPLATE meant for operator review - make this clear in header comments.
Only include steps the robot can physically do (pipetting, mixing, dispensing).
Do NOT include gel electrophoresis, imaging, incubation at specific temperatures
(unless using the temperature module), or manual steps.

Return ONLY the Python script as plain text. No markdown fences, no explanation."""


def _build_user_prompt(context: dict, templates: list) -> str:
    templates_json = json.dumps(templates, indent=2)
    context_json = json.dumps(context, indent=2)
    return (
        f"Context:\n{context_json}\n\n"
        f"Templates to fill (fill ALL of them and return as a JSON array in the same order):\n"
        f"{templates_json}\n\n"
        f"Fill all {{{{PLACEHOLDER}}}} fields. Return only the completed JSON array."
    )


def _build_ot2_prompt(context: dict, protocol_data: list) -> str:
    assay_types = ", ".join(p.get("assay_type", "") for p in protocol_data)
    lines = [
        "Generate an Opentrons OT-2 Python protocol for the following assay.",
        "",
        "Assay context:",
        json.dumps(context, indent=2),
        "",
        "Protocol summary:",
        f"- Assay type(s): {assay_types}",
        f"- Cell model: {context.get('cell_model', 'not specified')}",
        f"- Doses: {context.get('treatment_dose', 'not specified')}",
        f"- Timepoints: {context.get('treatment_timepoints', 'not specified')}",
        f"- Vehicle control: {context.get('vehicle_control', 'not specified')}",
        f"- Bio replicates: {context.get('bio_replicates', 3)}",
        f"- Gene: {context.get('gene', '')}",
        "",
        "Include:",
        "1. metadata block with protocolName, author='TTPicker Assay Generator', apiLevel='2.15'",
        "2. Labware definitions (tip racks, plates, reservoir)",
        "3. Pipette setup (p300_multi_gen2, p20_multi_gen2)",
        "4. Treatment dilution series preparation",
        "5. Cell treatment dispensing with replicate layout",
        "6. Vehicle control dispensing",
        "7. Clear comments marking automated vs manual steps",
        "",
        "Return only the Python script.",
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# API Call
# ---------------------------------------------------------------------------

def _call_claude_api(system: str, user_message: str, api_key: str) -> str:
    """Synchronous call to the Anthropic messages API."""
    import requests as _req
    resp = _req.post(
        "https://api.anthropic.com/v1/messages",
        headers={
            "Content-Type": "application/json",
            "x-api-key": api_key,
            "anthropic-version": "2023-06-01",
        },
        json={
            "model": "claude-opus-4-5-20251101",
            "max_tokens": 8192,
            "system": system,
            "messages": [{"role": "user", "content": user_message}],
        },
        timeout=120,
    )
    resp.raise_for_status()
    data = resp.json()
    return "".join(
        block.get("text", "")
        for block in data.get("content", [])
        if block.get("type") == "text"
    )


def _parse_json_response(raw: str) -> Any:
    """Parse JSON from Claude response, stripping markdown fences if present."""
    cleaned = raw.strip()
    cleaned = re.sub(r"^```(?:json)?\s*", "", cleaned)
    cleaned = re.sub(r"\s*```$", "", cleaned)
    return json.loads(cleaned)


# ---------------------------------------------------------------------------
# Protocol Rendering
# ---------------------------------------------------------------------------

def _esc(text: str) -> str:
    """Minimal HTML escaping."""
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
    )


def _render_protocol_html(protocols: list) -> str:
    """Convert filled protocol JSON to styled HTML."""
    html_parts = []
    for proto in protocols:
        html_parts.append('<div class="protocol-section">')
        html_parts.append(f'<h2 class="protocol-type">{_esc(proto.get("assay_type", ""))}</h2>')
        if proto.get("objective"):
            html_parts.append(f'<div class="protocol-objective"><strong>Objective:</strong> {_esc(proto["objective"])}</div>')
        if proto.get("rationale"):
            html_parts.append(f'<div class="protocol-rationale"><strong>Rationale:</strong> {_esc(proto["rationale"])}</div>')
        for section in proto.get("sections", []):
            html_parts.append('<div class="protocol-subsection">')
            html_parts.append(f'<h3>{_esc(section.get("title", ""))}</h3>')
            if "items" in section:
                html_parts.append("<ul>")
                for item in section["items"]:
                    html_parts.append(f"<li>{_esc(str(item))}</li>")
                html_parts.append("</ul>")
            if "steps" in section:
                html_parts.append("<ol>")
                for step in section["steps"]:
                    html_parts.append(f"<li>{_esc(str(step))}</li>")
                html_parts.append("</ol>")
            if "content" in section:
                html_parts.append(f'<div class="protocol-content">{_esc(str(section["content"]))}</div>')
            html_parts.append("</div>")
        html_parts.append("</div>")
    return "\n".join(html_parts)


def _render_protocol_text(protocols: list) -> str:
    """Convert filled protocol JSON to plain text."""
    lines = []
    for proto in protocols:
        lines.append("=" * 70)
        lines.append(proto.get("assay_type", ""))
        lines.append("=" * 70)
        if proto.get("objective"):
            lines.append(f"\nObjective: {proto['objective']}")
        if proto.get("rationale"):
            lines.append(f"Rationale: {proto['rationale']}")
        lines.append("")
        for section in proto.get("sections", []):
            lines.append(f"--- {section.get('title', '')} ---")
            for item in section.get("items", []):
                lines.append(f"  * {item}")
            for i, step in enumerate(section.get("steps", []), 1):
                lines.append(f"  {i}. {step}")
            if "content" in section:
                lines.append(f"  {section['content']}")
            lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# OT-2 Script Generation
# ---------------------------------------------------------------------------

def generate_ot2_script(context: dict, protocol_data: list, api_key: str) -> str:
    """Generate an Opentrons OT-2 Python script for the assay."""
    prompt = _build_ot2_prompt(context, protocol_data)
    raw = _call_claude_api(OT2_SYSTEM_PROMPT, prompt, api_key)
    script = raw.strip()
    script = re.sub(r"^```(?:python)?\s*", "", script)
    script = re.sub(r"\s*```$", "", script)
    return script


# ---------------------------------------------------------------------------
# Main Entry Point
# ---------------------------------------------------------------------------

def generate_protocol(pipeline_data: dict, user_answers: dict, api_key: str = None) -> dict:
    """
    Main function: generates validation assay protocol and OT-2 script.

    Args:
        pipeline_data: Full pipeline JSON output (lastData from frontend)
        user_answers:  Dict from the question form
        api_key:       Anthropic API key (defaults to ANTHROPIC_API_KEY env var)

    Returns:
        {
            "protocol_html":   str,
            "protocol_text":   str,
            "ot2_script":      str,
            "templates_used":  list,
        }
    """
    api_key = api_key or os.environ.get("ANTHROPIC_API_KEY", "")
    if not api_key:
        raise ValueError(
            "No Anthropic API key found. "
            "Set the ANTHROPIC_API_KEY environment variable before starting the server."
        )

    context = build_claude_context(pipeline_data, user_answers)
    equipment = user_answers.get("readout_equipment") or []
    template_keys = select_templates_with_equipment(pipeline_data, equipment)
    templates = [TEMPLATE_MAP[k] for k in template_keys]

    user_prompt = _build_user_prompt(context, templates)
    raw_response = _call_claude_api(SYSTEM_PROMPT, user_prompt, api_key)
    filled_protocols = _parse_json_response(raw_response)

    if isinstance(filled_protocols, dict):
        filled_protocols = [filled_protocols]

    protocol_html = _render_protocol_html(filled_protocols)
    protocol_text = _render_protocol_text(filled_protocols)
    ot2_script = generate_ot2_script(context, filled_protocols, api_key)

    return {
        "protocol_html":  protocol_html,
        "protocol_text":  protocol_text,
        "ot2_script":     ot2_script,
        "templates_used": template_keys,
    }

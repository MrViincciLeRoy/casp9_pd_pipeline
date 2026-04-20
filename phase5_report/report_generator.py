import json
from pathlib import Path
from datetime import datetime
from utils import get_logger
from config import RESULTS_DIR, TOP_N_CANDIDATES

log = get_logger(__name__)

def generate_report(final_candidates: list[dict]) -> Path:
    top = final_candidates[:TOP_N_CANDIDATES]
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    out = RESULTS_DIR / f"report_{ts}.md"

    lines = [
        "# Caspase-9 Inhibitor Candidates — PD Pipeline Report",
        f"\n**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"**Target:** Caspase-9 (PDB: 2AR9) — Parkinson's Disease",
        f"**Top candidates:** {len(top)}\n",
        "---\n",
        "## Top Candidates\n",
        "| Rank | SMILES | Vina Score (kcal/mol) | Predicted Ki (nM) | BBB |",
        "|------|--------|----------------------|-------------------|-----|",
    ]

    for i, c in enumerate(top, 1):
        smiles = c.get("smiles", "N/A")[:40]
        vina = c.get("vina_score", "N/A")
        ki = c.get("predicted_ki_nm", "N/A")
        bbb = "✅" if c.get("bbb_permeable") else "❌"
        lines.append(f"| {i} | `{smiles}` | {vina} | {ki} | {bbb} |")

    lines += [
        "\n---",
        "## Pipeline Summary",
        "- Phase 1: Structures downloaded from RCSB PDB, inhibitors from ChEMBL",
        "- Phase 2: Virtual screening with AutoDock Vina + GNINA validation",
        "- Phase 3: ML pre-screening with ChemBERTa-2 fine-tuned on Caspase-9 Ki data",
        "- Phase 4: ADMET filtering — Lipinski Rule of Five + BBB permeability (pkCSM)",
        "- Phase 5: Report generated for publication / wet lab collaboration",
        "\n## Next Steps",
        "- Upload to bioRxiv as preprint",
        "- Contact wet lab via ResearchGate or r/bioinformatics",
        "- Run MD simulations on top 5 hits to confirm stability",
    ]

    out.write_text("\n".join(lines))
    log.info(f"Report saved → {out}")
    return out

def run(candidates: list[dict]) -> Path:
    return generate_report(candidates)

import argparse
import json
import os
import sys
from pathlib import Path
from utils import get_logger

log = get_logger("pipeline")


def load_json(path: Path) -> list:
    if path and path.exists():
        data = json.loads(path.read_text())
        return data if isinstance(data, list) else []
    return []


def find_zinc_file(zinc_dir: Path) -> Path | None:
    matches = sorted(zinc_dir.glob("zinc_subset_*.json"))
    if matches:
        return matches[-1]
    return None


def run(phases: list[str]):
    from config import RESULTS_DIR, CHEMBL_DIR, ZINC_DIR

    log.info("=== CASP9 × PD PIPELINE STARTING ===")

    inhibitor_data = []
    zinc_smiles = []
    ml_hits = []
    admet_passed = []

    if "1" in phases:
        log.info("--- Phase 1: Data Gathering ---")
        from phase1_data import fetch_structures, fetch_chembl, fetch_zinc
        fetch_structures()
        fetch_chembl()
        fetch_zinc()

    if "2" in phases:
        log.info("--- Phase 2: Virtual Docking ---")
        log.info("Requires AutoDock Vina + GNINA installed locally.")

    if "3" in phases:
        log.info("--- Phase 3: ML Screening ---")
        inhibitor_data = load_json(CHEMBL_DIR / "casp9_inhibitors.json")
        if not inhibitor_data:
            log.error("No ChEMBL data — run Phase 1 first")
            sys.exit(1)

        zinc_file = find_zinc_file(ZINC_DIR)
        if not zinc_file:
            log.error(f"No ZINC file found in {ZINC_DIR} — run Phase 1 first")
            sys.exit(1)

        zinc_smiles = load_json(zinc_file)
        log.info(f"Loaded {len(zinc_smiles)} ZINC compounds from {zinc_file.name}")

        from phase3_ml import ml_screen
        ml_hits = ml_screen(inhibitor_data, zinc_smiles)
        log.info(f"Phase 3 complete — {len(ml_hits)} ML hits")

    if "4" in phases:
        log.info("--- Phase 4: ADMET Filtering ---")
        candidates = ml_hits or load_json(RESULTS_DIR / "ml_hits.json")
        if not candidates:
            log.warning("No ML hits — writing empty ADMET results and continuing")
            (RESULTS_DIR / "admet_passed.json").write_text("[]")
        else:
            from phase4_admet import run_admet_filter
            admet_passed = run_admet_filter(candidates)
            log.info(f"Phase 4 complete — {len(admet_passed)} candidates passed ADMET")

    if "5" in phases:
        log.info("--- Phase 5: Report Generation ---")
        final = admet_passed or load_json(RESULTS_DIR / "admet_passed.json")
        from phase5_report import generate_report
        report_path = generate_report(final)
        log.info(f"Pipeline complete. Report: {report_path}")

    log.info("=== PIPELINE DONE ===")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Caspase-9 × PD Drug Discovery Pipeline")
    parser.add_argument("--phases", nargs="+", default=["1", "3", "4", "5"])
    args = parser.parse_args()
    run(args.phases)

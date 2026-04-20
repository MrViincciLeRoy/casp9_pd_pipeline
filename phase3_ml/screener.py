import json
import torch
from pathlib import Path
from torch.utils.data import DataLoader
from torch.optim import AdamW
from utils import get_logger
from config import RESULTS_DIR, EPOCHS, LEARNING_RATE, BATCH_SIZE, ML_KI_CUTOFF_NM
from .model import CaspaseInhibitorModel
from .dataset import CaspaseDataset

log = get_logger(__name__)
MODEL_PATH = RESULTS_DIR / "casp9_model.pt"

def train(smiles_list: list[str], ki_values: list[float]):
    log.info(f"Training on {len(smiles_list)} compounds...")
    dataset = CaspaseDataset(smiles_list, ki_values)
    loader = DataLoader(dataset, batch_size=BATCH_SIZE, shuffle=True)
    model = CaspaseInhibitorModel()
    optimizer = AdamW(model.parameters(), lr=LEARNING_RATE)
    loss_fn = torch.nn.MSELoss()

    model.train()
    for epoch in range(EPOCHS):
        total_loss = 0
        for batch in loader:
            optimizer.zero_grad()
            preds = model(batch["input_ids"], batch["attention_mask"])
            loss = loss_fn(preds, batch["label"])
            loss.backward()
            optimizer.step()
            total_loss += loss.item()
        log.info(f"Epoch {epoch+1}/{EPOCHS} — loss: {total_loss/len(loader):.4f}")

    model.save(str(MODEL_PATH))
    log.info(f"Model saved → {MODEL_PATH}")
    return model

def screen(smiles_list: list[str]) -> list[dict]:
    log.info(f"Screening {len(smiles_list)} compounds via ML model...")
    model = CaspaseInhibitorModel()
    model.load(str(MODEL_PATH))
    dataset = CaspaseDataset(smiles_list, [0.0] * len(smiles_list))
    loader = DataLoader(dataset, batch_size=BATCH_SIZE)
    predictions = []

    with torch.no_grad():
        for batch in loader:
            preds = model(batch["input_ids"], batch["attention_mask"]).tolist()
            predictions.extend(preds)

    results = [
        {"smiles": s, "predicted_ki_nm": round(p, 2)}
        for s, p in zip(smiles_list, predictions)
        if p < ML_KI_CUTOFF_NM
    ]
    results.sort(key=lambda x: x["predicted_ki_nm"])
    out = RESULTS_DIR / "ml_hits.json"
    out.write_text(json.dumps(results, indent=2))
    log.info(f"{len(results)} ML hits → {out}")
    return results

def run(inhibitor_data: list[dict], zinc_smiles: list[str]):
    smiles = [d["smiles"] for d in inhibitor_data]
    ki_vals = [d["ki_nm"] for d in inhibitor_data]
    train(smiles, ki_vals)
    return screen(zinc_smiles)

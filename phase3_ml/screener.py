import gc
import json
import os
import torch
from torch.utils.data import DataLoader
from torch.optim import AdamW
from utils import get_logger
<<<<<<< HEAD
from config import RESULTS_DIR, EPOCHS, LEARNING_RATE, BATCH_SIZE
=======
from config import RESULTS_DIR, BATCH_SIZE, LEARNING_RATE
>>>>>>> 9639c8b (env vars)

log = get_logger(__name__)
MODEL_PATH = RESULTS_DIR / "casp9_model.pt"

CI = os.getenv("CI", "false").lower() == "true"
<<<<<<< HEAD
EFFECTIVE_BATCH = 16 if CI else BATCH_SIZE
KI_CUTOFF = float(os.getenv("ML_KI_CUTOFF_NM", "100"))
ZINC_LIMIT = int(os.getenv("ZINC_SCREEN_LIMIT", "0")) or None


def train(smiles_list: list[str], ki_values: list[float]):
    log.info(f"Training on {len(smiles_list)} compounds (batch={EFFECTIVE_BATCH})...")
=======
KI_CUTOFF = float(os.getenv("ML_KI_CUTOFF_NM", "100"))
ZINC_LIMIT = int(os.getenv("ZINC_SCREEN_LIMIT", "0")) or None
EPOCHS = int(os.getenv("EPOCHS", "10"))
EFFECTIVE_BATCH = int(os.getenv("BATCH_SIZE", str(16 if CI else BATCH_SIZE)))


def train(smiles_list: list[str], ki_values: list[float]):
    log.info(f"Training on {len(smiles_list)} compounds | epochs={EPOCHS} | batch={EFFECTIVE_BATCH}")
>>>>>>> 9639c8b (env vars)
    from .dataset import CaspaseDataset
    from .model import CaspaseInhibitorModel

    dataset = CaspaseDataset(smiles_list, ki_values)
    loader = DataLoader(dataset, batch_size=EFFECTIVE_BATCH, shuffle=True)
    model = CaspaseInhibitorModel()
    optimizer = AdamW(model.parameters(), lr=LEARNING_RATE)
    loss_fn = torch.nn.MSELoss()

    model.train()
    for epoch in range(EPOCHS):
        total_loss = 0.0
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
    del model, optimizer, loader, dataset
    gc.collect()


def screen(smiles_list: list[str]) -> list[dict]:
    from .dataset import CaspaseDataset
    from .model import CaspaseInhibitorModel

    if ZINC_LIMIT:
        smiles_list = smiles_list[:ZINC_LIMIT]

    log.info(f"Screening {len(smiles_list)} compounds | Ki cutoff={KI_CUTOFF}nM | batch={EFFECTIVE_BATCH}")
<<<<<<< HEAD

=======
>>>>>>> 9639c8b (env vars)
    model = CaspaseInhibitorModel()
    model.load(str(MODEL_PATH))
    dataset = CaspaseDataset(smiles_list, [0.0] * len(smiles_list))
    loader = DataLoader(dataset, batch_size=EFFECTIVE_BATCH)
    predictions = []

    with torch.no_grad():
        for batch in loader:
            preds = model(batch["input_ids"], batch["attention_mask"]).tolist()
            predictions.extend(preds)

    results = [
        {"smiles": s, "predicted_ki_nm": round(p, 2)}
        for s, p in zip(smiles_list, predictions)
        if p < KI_CUTOFF
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
import torch
from torch.utils.data import Dataset
from transformers import AutoTokenizer
from config import CHEMBERTA_MODEL

tokenizer = AutoTokenizer.from_pretrained(CHEMBERTA_MODEL)

class CaspaseDataset(Dataset):
    def __init__(self, smiles_list: list[str], ki_values: list[float]):
        self.smiles = smiles_list
        self.labels = ki_values

    def __len__(self):
        return len(self.smiles)

    def __getitem__(self, idx):
        enc = tokenizer(self.smiles[idx], return_tensors="pt", padding="max_length",
                        truncation=True, max_length=512)
        return {
            "input_ids": enc["input_ids"].squeeze(0),
            "attention_mask": enc["attention_mask"].squeeze(0),
            "label": torch.tensor(self.labels[idx], dtype=torch.float),
        }

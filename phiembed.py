"""
========================================================================
This script is for running PHIEmbed. It takes a FASTA file as input
and outputs the predicted host genus, along with a list of probabilities 
for each host genus in our dataset.

@author    Mark Edward M. Gonzales
========================================================================
"""

import argparse
import os
import re

from Bio import SeqIO
import joblib
import torch
from transformers import T5Tokenizer, T5EncoderModel

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
tokenizer = T5Tokenizer.from_pretrained(
    "Rostlab/prot_t5_xl_half_uniref50-enc", do_lower_case=False
)
model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_half_uniref50-enc").to(
    device
)


# Adapted from https://github.com/agemagician/ProtTrans
def embed(sequence):
    sequences = [sequence]
    sequences = [
        " ".join(list(re.sub(r"[UZOB]", "X", sequence))) for sequence in sequences
    ]
    ids = tokenizer(sequences, add_special_tokens=True, padding="longest")

    input_ids = torch.tensor(ids["input_ids"]).to(device)
    attention_mask = torch.tensor(ids["attention_mask"]).to(device)

    with torch.no_grad():
        embedding_repr = model(input_ids=input_ids, attention_mask=attention_mask)

    embedding = embedding_repr.last_hidden_state[0, : len(sequence)]
    return embedding.mean(dim=0).tolist()


def predict(embedding, clf):
    top_prediction = clf.predict([embedding])
    classes = clf.classes_
    proba = clf.predict_proba([embedding])

    scores = []
    for idx, class_name in enumerate(classes):
        scores.append((class_name, proba[0][idx]))

    return sorted(scores, key=lambda x: x[1], reverse=True)


def write_results(id, scores, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(f"{output_dir}/{id}.csv", "w") as f:
        for entry in scores:
            f.write(f"{entry[0]},{entry[1]}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        required=True,
        help="Path to the FASTA file containing the receptor-binding protein sequences",
    )

    parser.add_argument(
        "--model",
        required=True,
        help="Path to the trained model (recognized format: joblib, framework: scikit-learn)",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Path to the file to which the results of running PHIEmbed will be written",
    )

    args = parser.parse_args()

    clf = joblib.load(args.model)
    for record in SeqIO.parse(args.input, "fasta"):
        write_results(str(record.id), predict(embed(str(record.seq)), clf), args.output)

"""
=======================================================================
This script is for training PHIEmbed. It takes a CSV file corresponding
to the training dataset as input and outputs a trained scikit-learn
random forest classifier (serialized in joblib format).

@author    Mark Edward M. Gonzales
=======================================================================
"""

import argparse

import joblib
import pandas as pd
from sklearn.ensemble import RandomForestClassifier

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        required=True,
        help="Path to the training dataset",
    )

    parser.add_argument(
        "--threads",
        help="Number of threads to be used for training (default: -1, that is, use all threads)",
        type=int,
        default=-1,
    )

    args = parser.parse_args()

    train = pd.read_csv(
        args.input,
        header=None,
        names=["Protein ID", "Host"] + [str(i) for i in range(1, 1025)],
    )
    X_train = train.loc[:, train.columns.isin([str(i) for i in range(1, 1025)])]
    y_train = train.loc[:, train.columns.isin(["Host"])]

    assert X_train.shape[1] == 1024 and y_train.shape[1] == 1

    clf = RandomForestClassifier(
        class_weight="balanced",
        max_features="sqrt",
        min_samples_leaf=1,
        min_samples_split=2,
        n_estimators=150,
        n_jobs=args.threads,
        verbose=True,
    )

    clf.fit(X_train.values, y_train.values.ravel())
    joblib.dump(clf, "phiembed_trained.joblib.gz", compress=True)

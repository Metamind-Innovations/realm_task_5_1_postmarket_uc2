import argparse
import json
import os

from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
import pandas as pd


phenotype_map = {
    "NM": 0,
    "LNM": 1,
    "IM": 2,
    "LIM": 3,
    "PM": 4,
    "LPM": 5,
    "UM": 6,
    "LUM": 7,
    "RM": 8,
    "LRM": 9,
    "INDETERMINATE": -1,
    "NF": 10,
    "DF": 11,
    "IF": 12,
    "PF": 13,
    "PDF": 14,
}


def preprocess_groundtruth(csv_file_path):
    """
    Preprocesses a groundtruth phenotype CSV file and transforms phenotype values
    to numeric values based on a predefined mapping.

    Args:
        csv_file_path (str): Path to the groundtruth CSV file with phenotype data

    Returns:
        pd.DataFrame: DataFrame with transformed phenotype values
    """
    try:
        groundtruth_df = pd.read_csv(csv_file_path, sep=",", index_col="Unnamed: 0")
    except (ValueError, KeyError):
        groundtruth_df = pd.read_csv(csv_file_path, sep=",", index_col="Unnamed: 0")

    if "Sample ID" in groundtruth_df.columns:
        groundtruth_df = groundtruth_df.rename(columns={"Sample ID": "Sample"})

    phenotype_columns = [
        col for col in groundtruth_df.columns if col not in ["Unnamed: 0", "Sample"]
    ]

    for col in phenotype_columns:
        groundtruth_df[col] = groundtruth_df[col].map(phenotype_map)

    return groundtruth_df


def preprocess_predictions(csv_file_path):
    """
    Preprocesses a predictions phenotype CSV file and transforms phenotype values
    to numeric values based on a predefined mapping.

    Args:
        csv_file_path (str): Path to the predictions CSV file with phenotype data

    Returns:
        pd.DataFrame: DataFrame with transformed phenotype values
    """
    predictions_df = pd.read_csv(csv_file_path)

    if "Sample ID" in predictions_df.columns:
        predictions_df = predictions_df.rename(columns={"Sample ID": "Sample"})

    phenotype_columns = [col for col in predictions_df.columns if col not in ["Sample"]]

    for col in phenotype_columns:
        predictions_df[col] = predictions_df[col].map(phenotype_map)

    return predictions_df


def evaluate_phenotype_predictions(groundtruth_file, predictions_file):
    """
    Evaluates phenotype predictions against ground truth data.

    Args:
        groundtruth_file (str): Path to the ground truth CSV file
        predictions_file (str): Path to the predictions CSV file

    Returns:
        pd.DataFrame: DataFrame with evaluation metrics for each phenotype
    """
    groundtruth_df = preprocess_groundtruth(groundtruth_file)
    predictions_df = preprocess_predictions(predictions_file)

    # Ensure both dataframes have the same samples in the same order
    if "Sample" in groundtruth_df.columns and "Sample" in predictions_df.columns:
        groundtruth_df = groundtruth_df.set_index("Sample")
        predictions_df = predictions_df.set_index("Sample")
        groundtruth_df, predictions_df = groundtruth_df.align(
            predictions_df, join="inner"
        )

    phenotype_columns = [
        col
        for col in groundtruth_df.columns
        if col not in ["Unnamed: 0", "Sample"] and col in predictions_df.columns
    ]

    results = {
        "Phenotype": [],
        "Accuracy": [],
        "Precision": [],
        "Recall": [],
        "F1-Score": [],
        "Samples": [],
    }

    # Calculate metrics for each phenotype
    for col in phenotype_columns:
        y_true = groundtruth_df[col].dropna()
        y_pred = predictions_df.loc[y_true.index, col]

        valid_indices = y_true.index.intersection(y_pred.dropna().index)
        y_true = y_true.loc[valid_indices]
        y_pred = y_pred.loc[valid_indices]

        if len(y_true) == 0:
            continue

        accuracy = accuracy_score(y_true, y_pred)

        # Using 'weighted' to account for class imbalance
        precision = precision_score(y_true, y_pred, average="weighted", zero_division=0)
        recall = recall_score(y_true, y_pred, average="weighted", zero_division=0)
        f1 = f1_score(y_true, y_pred, average="weighted", zero_division=0)

        results["Phenotype"].append(col)
        results["Accuracy"].append(accuracy)
        results["Precision"].append(precision)
        results["Recall"].append(recall)
        results["F1-Score"].append(f1)
        results["Samples"].append(len(y_true))

    results_df = pd.DataFrame(results)

    return results_df


def compare_evaluation_results(rwd_results_df, synthetic_results_df):
    """
    Compare evaluation metrics between RWD and synthetic results.

    Args:
        rwd_results_df (pd.DataFrame): Results dataframe from RWD evaluation
        synthetic_results_df (pd.DataFrame): Results dataframe from synthetic data evaluation

    Returns:
        dict: Comparison of metrics and their differences
    """
    comparison = {}

    phenotypes = set(rwd_results_df["Phenotype"]).union(
        set(synthetic_results_df["Phenotype"])
    )

    metric_columns = ["Accuracy", "Precision", "Recall", "F1-Score"]

    comparison["summary"] = {
        metric: {"rwd": [], "synthetic": [], "difference": []}
        for metric in metric_columns
    }

    for phenotype in phenotypes:
        comparison[phenotype] = {}

        rwd_row = rwd_results_df[rwd_results_df["Phenotype"] == phenotype]
        synthetic_row = synthetic_results_df[
            synthetic_results_df["Phenotype"] == phenotype
        ]

        for metric in metric_columns:
            rwd_value = float(rwd_row[metric].iloc[0]) if not rwd_row.empty else None
            synthetic_value = (
                float(synthetic_row[metric].iloc[0])
                if not synthetic_row.empty
                else None
            )

            comparison[phenotype][metric] = {
                "rwd": rwd_value,
                "synthetic": synthetic_value,
            }

            if rwd_value is not None and synthetic_value is not None:
                difference = synthetic_value - rwd_value
                comparison[phenotype][metric]["difference"] = difference

                if metric != "Samples":
                    comparison["summary"][metric]["rwd"].append(rwd_value)
                    comparison["summary"][metric]["synthetic"].append(synthetic_value)
                    comparison["summary"][metric]["difference"].append(difference)
            else:
                comparison[phenotype][metric]["difference"] = None

    # Calculate averages for summary statistics
    for metric in metric_columns:
        if metric != "Samples":
            if comparison["summary"][metric]["rwd"]:
                comparison["summary"][metric]["rwd_avg"] = sum(
                    comparison["summary"][metric]["rwd"]
                ) / len(comparison["summary"][metric]["rwd"])
            else:
                comparison["summary"][metric]["rwd_avg"] = None

            if comparison["summary"][metric]["synthetic"]:
                comparison["summary"][metric]["synthetic_avg"] = sum(
                    comparison["summary"][metric]["synthetic"]
                ) / len(comparison["summary"][metric]["synthetic"])
            else:
                comparison["summary"][metric]["synthetic_avg"] = None

            if comparison["summary"][metric]["difference"]:
                comparison["summary"][metric]["difference_avg"] = sum(
                    comparison["summary"][metric]["difference"]
                ) / len(comparison["summary"][metric]["difference"])
            else:
                comparison["summary"][metric]["difference_avg"] = None

    for metric in ["Accuracy", "Precision", "Recall", "F1-Score"]:
        del comparison["summary"][metric]["rwd"]
        del comparison["summary"][metric]["synthetic"]
        del comparison["summary"][metric]["difference"]

    return comparison


def main():
    parser = argparse.ArgumentParser(
        description="Transform phenotype values in a CSV file"
    )
    parser.add_argument(
        "--groundtruth_file", required=True, help="Path to the input CSV file"
    )
    parser.add_argument(
        "--rwd_predictions_file", required=True, help="Path to the predictions CSV file"
    )
    parser.add_argument(
        "--synthetic_predictions_file",
        required=True,
        help="Path to the predictions CSV file",
    )
    parser.add_argument(
        "--output_file",
        default="artifacts/adversarial_evaluation.json",
        help="Path to the output JSON file",
    )
    args = parser.parse_args()

    os.makedirs("artifacts", exist_ok=True)

    rwd_results_df = evaluate_phenotype_predictions(
        args.groundtruth_file, args.rwd_predictions_file
    )
    synthetic_results_df = evaluate_phenotype_predictions(
        args.groundtruth_file, args.synthetic_predictions_file
    )

    comparison = compare_evaluation_results(rwd_results_df, synthetic_results_df)

    with open(args.output_file, "w") as f:
        json.dump(comparison, f, indent=4)


if __name__ == "__main__":
    main()

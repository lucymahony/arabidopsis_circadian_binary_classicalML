import pandas as pd
import numpy as np
#import xgboost as xgb
import lightgbm as lgb
from sklearn.model_selection import StratifiedKFold
import argparse 
from sklearn.feature_selection import SelectKBest, chi2
from sklearn.metrics import f1_score
from utils import load_data, preprocess_data, split_train_test, evaluate_model_by_class

# Chi2 feature selection. Fivefold cross-validation for each set of reduced features.
# Plot mean F1score from the test dataset as the feature set is reduced 


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Train a model using k-mer counts.")
    parser.add_argument("--cv", type=int, default=5, help="Number of cross-validation folds.")
    parser.add_argument("--test_size_split", type=float, default=0.1, help="Proportion of data to use for testing.")
    parser.add_argument("--input_data", type=str,
                        default="/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/arabidopsis_circadian_binary_classicalML/data/processed/promoters_transcripts_6_mers_labelled.csv",
                        help="Path to the input data file.")
    return parser.parse_args()

def mean_f1_score(y_true, y_pred):
    """Calculate the mean F1 score for two classes."""
    f1_class0 = f1_score(y_true, y_pred, pos_label=0)
    f1_class1 = f1_score(y_true, y_pred, pos_label=1)
    return (f1_class0 + f1_class1) / 2

def feature_selection(model, X_train, y_train, X_test, y_test):
    results = []
    for k in [8192, 8000, 7000, 6000, 5000, 4000, 3000, 2000, 1000, 500, 100, 50, 10]:
        selector = SelectKBest(score_func=chi2, k=k)
        X_train_new = selector.fit_transform(X_train, y_train)
        X_test_new = selector.transform(X_test)

        model.fit(X_train_new, y_train)
        y_pred = model.predict(X_test_new)
        avg_f1 = mean_f1_score(y_test, y_pred)

        print(f"Number of features: {k}")
        print(f"Mean F1 score: {avg_f1:.4f}")
        print("--------------------------------------------------------")

        results.append({
            'num_features': k,
            'mean_f1': avg_f1
        })

    results_df = pd.DataFrame(results)
    results_df.to_csv('/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/arabidopsis_circadian_binary_classicalML/data/processed/feature_selection_results.csv', index=False)


def main():
    random_seed = 0 
    args = parse_args()
    X, y = load_data(args.input_data)
    X = preprocess_data(X)
    X_train, X_test, y_train, y_test = split_train_test(X, y, args.test_size_split, random_seed)

    model = lgb.LGBMClassifier(num_leaves=100, subsample=0.5, max_depth=9, learning_rate=0.13853569656385206, n_estimators=350)

    feature_selection(model, X_train, y_train, X_test, y_test)


if __name__ == "__main__":
    main()

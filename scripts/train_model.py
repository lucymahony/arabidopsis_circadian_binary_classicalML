import pandas as pd
import numpy as np
#import xgboost as xgb
import lightgbm as lgb
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import WhiteKernel, RBF
from scipy.stats import expon
import argparse 

from utils import load_data, preprocess_data, split_train_test, evaluate_model_by_class



def parse_args():
    """Parse command-line arguments."""
    
    parser = argparse.ArgumentParser(description="Train a model using k-mer counts.")
    parser.add_argument("--models", type=str, required=True,
                        choices=['logistic_regression', 'random_forest', 'svm', 'knn', 'xgboost', 'decision_tree', 'lightgbm', 'gaussian_process', 'all'],
                        default="all", help="Which models to train.")
    parser.add_argument("--cv", type=int, default=5, help="Number of cross-validation folds.")
    parser.add_argument("--test_size_split", type=float, default=0.1, help="Proportion of data to use for testing.")
    parser.add_argument("--input_data", type=str,
                        default="/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/arabidopsis_circadian_binary_classicalML/data/processed/promoters_transcripts_6_mers_labelled.csv",
                        help="Path to the input data file.")
    return parser.parse_args()


def get_models_and_grids(model_names, random_seed=0):
    # Distributions  - The number of samples not specified, so set to 3 for all except learning_rate_values which are important therefore 10
    gamma_values = sorted(expon(scale=0.1).rvs(size=3, random_state=random_seed))
    c_values = sorted(expon(scale=100).rvs(size=3, random_state=random_seed))
    learning_rate_values = sorted(expon(scale=1).rvs(size=10, random_state=random_seed)) 
    min_child_weight_values = sorted(expon(scale=10).rvs(size=3, random_state=random_seed))

    models_and_grids = {
        'logistic_regression': (
            LogisticRegression(solver='liblinear'),
            {'penalty': ['l1', 'l2'], 'C': [1.0, 0.5, 0.1]}
        ),
        'random_forest': (
            RandomForestClassifier(),
            {
                'n_estimators': list(np.linspace(100, 4000, 10, dtype=int)),
                'max_depth': list(range(1, 11)),
                'min_samples_leaf': list(range(1, 11)),
                'min_samples_split': [2, 5, 10],
                'max_features': ['auto', 'sqrt'],
                'criterion': ['gini', 'entropy']
            }
        ),
        'svm': (
            SVC(),
            {'kernel': ['linear', 'rbf'],
             'C': c_values,
             'gamma': gamma_values}
        ),
        'knn': (
            KNeighborsClassifier(),
            {
                'n_neighbors': list(range(1, 21)),
                'leaf_size': list(range(1, 6)),
                'weights': ['uniform', 'distance'],
                'algorithm': ['auto', 'ball_tree', 'kd_tree', 'brute']
            }
        ),
        #'xgboost': (
        #    xgb.XGBClassifier(use_label_encoder=False, eval_metric='logloss'),
        #    {
        #        'max_depth': [2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20],
        #        'subsample': [0.2, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
        #        'learning_rate': learning_rate_values,
        #        'min_child_weight': min_child_weight_values,
        #        'max_delta_step': [0, 1, 2]
        #    }
        #),
        'decision_tree': (
            DecisionTreeClassifier(),
            {
                'criterion': ['gini', 'entropy'],
                'min_samples_leaf': [1, 2, 3, 4, 5],
                'max_depth': [1, 2, 3, 4, 5],
                'min_samples_split': [2, 3, 4, 5]
            }
        ),
        'lightgbm': (
            lgb.LGBMClassifier(),
            {
                'num_leaves': [10, 20, 50, 100, 200],
                'subsample': [0.2, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
                'min_data_in_leaf': [10, 25, 50, 75, 100],
                'max_depth': [3, 5, 6, 7, 8, 10, 15, 20, 25],
                'learning_rate': learning_rate_values,
                'n_estimators': [50, 150, 200, 250, 300, 350, 400, 450, 500]
            }
        ),
        'gaussian_process': (
            GaussianProcessClassifier(kernel=WhiteKernel(noise_level=1e-7) + RBF()),
            {
                'normalize_y': [True, False],
                'copy_X_train': [True, False],
                'alpha': [1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12]
            }
        )
    }
    if model_names == 'all':
        return models_and_grids
    selected = model_names.split(',')
    return {k: models_and_grids[k] for k in selected if k in models_and_grids}


def train_and_tune_model(model, param_grid, X_train, y_train, cv):
    grid = GridSearchCV(model, param_grid, cv=StratifiedKFold(n_splits=cv), scoring='f1_weighted', n_jobs=-1)
    grid.fit(X_train, y_train)
    return grid.best_estimator_, grid.best_score_



def main():
    random_seed = 0 
    args = parse_args()
    X, y = load_data(args.input_data)
    X = preprocess_data(X)
    X_train, X_test, y_train, y_test = split_train_test(X, y, args.test_size_split, random_seed)

    models_and_grids = get_models_and_grids(args.models)

    results = []
    all_cv_results = []

    for name, (model, grid) in models_and_grids.items():
        print(f"\nTraining model: {name}")
        best_model, mean_cv_f1 = train_and_tune_model(model, grid, X_train, y_train, args.cv)

        eval_scores = evaluate_model_by_class(best_model, X_train, y_train, X_test, y_test)
        eval_scores["Model"] = name
        eval_scores["Mean F1 CV"] = mean_cv_f1
        results.append(eval_scores)

        parameter_grid = eval_scores.copy()
        parameter_grid["Parameters"] = str(grid)
        all_cv_results.append(parameter_grid)

    results_df = pd.DataFrame(results)
    results_df.to_csv("/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/arabidopsis_circadian_binary_classicalML/model_f1_scores_by_class.csv", index=False)
    final_cv_results = pd.concat(all_cv_results, ignore_index=True)
    final_cv_results.to_csv("/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/arabidopsis_circadian_binary_classicalML/all_models_gridsearch_results.csv", index=False)

    results_df.sort_values(by="Test F1 Class 1", ascending=False, inplace=True)
    print("\nFinal Results (sorted by Test F1 Class 1):")
    print(results_df)


if __name__ == "__main__":
    main()

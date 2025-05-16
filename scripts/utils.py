
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import f1_score, classification_report


def load_data(input_data_file_path):
    df = pd.read_csv(input_data_file_path)
    X = df.drop(['label', 'transcript'], axis=1)
    y = df['label']
    return X, y


def preprocess_data(X):
    scaler = MinMaxScaler() # default is 0-1 
    X_scaled = scaler.fit_transform(X)
    return X_scaled


def split_train_test(X, y, test_size_split=0.1, random_seed=0):
    return train_test_split(X, y, test_size=test_size_split, stratify=y, random_state=random_seed)


def evaluate_model_by_class(model, X_train, y_train, X_test, y_test):
    y_train_pred = model.predict(X_train)
    y_test_pred = model.predict(X_test)

    train_f1_class0 = f1_score(y_train, y_train_pred, pos_label=0)
    train_f1_class1 = f1_score(y_train, y_train_pred, pos_label=1)
    test_f1_class0 = f1_score(y_test, y_test_pred, pos_label=0)
    test_f1_class1 = f1_score(y_test, y_test_pred, pos_label=1)

    return {
        "Train F1 Class 0": train_f1_class0,
        "Train F1 Class 1": train_f1_class1,
        "Test F1 Class 0": test_f1_class0,
        "Test F1 Class 1": test_f1_class1,
    }


def evaluate_model(model, X_test, y_test):
    y_pred = model.predict(X_test)
    score = f1_score(y_test, y_pred, average='weighted')
    print(classification_report(y_test, y_pred))
    return score

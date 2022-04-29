import sys
import numpy as np
sys.path.append("../../../../")
from simCRN.multivariate_reg import (
    read_eq_data_file
)
from sklearn.model_selection import train_test_split
import pandas as pd

FILENAME = "4-4-2-semi-asym-AB-AC"  # Note: Adjust this as needed.
Ci_all_array, Am_array, Cmin, Cmax, Ai = read_eq_data_file(f"../../../{FILENAME}.txt")

# I will use this script to convert things into CSV format.
# This is because my TensorFlow workflows take CSVs.

# 70% train, 10% validation, and 20% test
random_states = [0, 1, 2, 3, 4]
for my_rand in random_states: 
    X_train, X_test, y_train, y_test = train_test_split(
        Am_array, Ci_all_array, test_size=0.2, random_state=0 # Using the same test set across the random trials.
    )
    X_train, X_val, y_train, y_val = train_test_split(
        X_train, y_train, test_size=0.125, random_state=my_rand
    )

    num_train = X_train.shape[0]
    num_val = X_val.shape[0]
    num_test = X_test.shape[0]

    print(num_train, num_val, num_test)

    # Putting together a Pandas dataframe and saving it as a CSV file.
    data = {
        "A1": np.concatenate((X_train[:, 0], X_val[:, 0], X_test[:, 0])),
        "A2": np.concatenate((X_train[:, 1], X_val[:, 1], X_test[:, 1])),
        "A3": np.concatenate((X_train[:, 2], X_val[:, 2], X_test[:, 2])),
        "A4": np.concatenate((X_train[:, 3], X_val[:, 3], X_test[:, 3])),
        "C1": np.concatenate((y_train[:, 0], y_val[:, 0], y_test[:, 0])),
        "C2": np.concatenate((y_train[:, 1], y_val[:, 1], y_test[:, 1])),
        "bin": ["train"] * num_train + ["val"] * num_val + ["test"] * num_test,
    }
    my_df = pd.DataFrame(data=data)

    my_df.to_csv(f"{FILENAME}_rs_{my_rand}.csv", index=False)

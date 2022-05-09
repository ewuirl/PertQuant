import pandas as pd
import sklearn.preprocessing
import sklearn.utils
import numpy as np
import sklearn.gaussian_process.kernels
import keras
import matplotlib.pyplot as plt

def figure_formatting(): # make the plot look xmgrace-esque
    font = {'family': 'sans-serif', 'weight': 'bold', 'size': 18}
    plt.rcParams['font.sans-serif'] = ['Arial']
    plt.rc('font', **font)
    plt.rcParams['axes.linewidth'] = 2.5
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.major.width'] = 2.5
    plt.rcParams['ytick.major.size'] = 10
    plt.rcParams['ytick.major.width'] = 2.5
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['mathtext.default'] = 'regular'
    plt.rcParams['legend.fancybox'] = False # no rounded legend box

figure_formatting()

descriptors = ["A1", "A2", "A3", "A4"]  # Note: change as needed


def normalize_data(df_train, df_test, fnames, lname, unit_trans=1, debug=False):
    _df_train = df_train.copy().dropna(subset=fnames + lname)
    _df_test = df_test.copy().dropna(subset=fnames + lname)
    X_train, X_test = _df_train[fnames].values, _df_test[fnames].values
    y_train, y_test = (
        _df_train[lname].values * unit_trans,
        _df_test[lname].values * unit_trans,
    )
    if debug:
        print(
            "training data reduced from %d -> %d because of nan."
            % (len(df_train), y_train.shape[0])
        )
        print(
            "test data reduced from %d -> %d because of nan."
            % (len(df_test), y_test.shape[0])
        )
    x_scaler = sklearn.preprocessing.StandardScaler()
    x_scaler.fit(X_train)
    X_train = x_scaler.transform(X_train)
    X_test = x_scaler.transform(X_test)
    y_scaler = sklearn.preprocessing.StandardScaler()
    y_scaler.fit(y_train)
    y_train = y_scaler.transform(y_train)
    y_test = y_scaler.transform(y_test)
    for ii, ln in enumerate(lname):
        if debug:
            print(
                (
                    "mean in target: %s, train/test: %f/%f"
                    % (ln, np.mean(y_train[:, ii]), np.mean(y_test[:, ii]))
                )
            )
    return X_train, X_test, y_train, y_test, x_scaler, y_scaler


def main():

    df = pd.read_csv("4-4-2-asym-AB-AC_rs_0.csv")
    df_train = df[df["bin"] == "train"]
    df_val = df[df["bin"] == "val"]
    df_train_all = df_train.append(df_val)
    df_test = df[df["bin"] == "test"]

    fnames = descriptors
    lname = ["C1"]  # Note: change as needed

    X_train, X_test, y_train, y_test, x_scaler, y_scaler = normalize_data(
        df_train_all, df_test, fnames, lname, unit_trans=1, debug=True
    )

    model = keras.models.load_model(
        f"quant_ANN_{lname[0]}.h5", compile=False
    )  # Note: change as needed

    predictions = y_scaler.inverse_transform(model.predict(X_test))

    ### Will append the predictions to the dataframe df_test, and save to an Excel file

    df_predictions = pd.DataFrame(
        predictions, columns=["Prediction"]
    )  # making the predictions into a pandas dataframe, in order to append
    df_predictions.reset_index(drop=True, inplace=True)
    df_test.reset_index(drop=True, inplace=True)

    df_predictions = pd.concat(
        [df_predictions, df_test], axis=1
    )  # combining test data with the predictions, column wise

    # Reordering the columns of the dataframe
    column_order = ["C1", "C2", "Prediction"] + descriptors
    df_predictions = df_predictions[column_order]

    # Writing it all to a csv

    df_predictions.to_csv(f"predictions_{lname[0]}.csv", index=False)

    ### Making a parity plot

    truth_values = df_predictions[lname[0]]  # Note: change as needed
    prediction_values = df_predictions["Prediction"]

    truth_values = [item * 1E6 for item in truth_values] # Converting to micromolar
    prediction_values = [item * 1E6 for item in prediction_values]

    plt.scatter(truth_values, prediction_values)
    plt.xlabel('Truth (\u03BCM)', fontweight='bold')
    plt.ylabel('Prediction (\u03BCM)', fontweight='bold')
    plt.xticks([0.00, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50])
    plt.title('Parity plot for P$_1$', fontweight='bold') # Note: change as needed
    plt.savefig(f'{lname[0]}_parityplot.png', format='png', bbox_inches='tight', dpi=300)
    plt.savefig(f'{lname[0]}_parityplot.eps', format='eps', bbox_inches='tight') # saving as a vector plot


if __name__ == "__main__":
    main()

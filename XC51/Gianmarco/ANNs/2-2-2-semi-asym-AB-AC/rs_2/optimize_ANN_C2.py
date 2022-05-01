import pandas as pd
import sklearn.preprocessing
import sklearn.utils
import numpy as np
import sklearn.gaussian_process.kernels
import pickle
from sklearn.metrics import r2_score, mean_absolute_error
from hyperopt import hp, tpe, fmin, Trials
from functools import partial
from molSimplifyAD.retrain.nets import build_ANN
from molSimplifyAD.retrain.model_optimization import train_model_hyperopt

descriptors = ["A1", "A2"] # Note: change as needed.

# Z-normalization.
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


# This function is used for hyperparameter optimization.
def optimize(
    X,
    y,
    lname,
    regression=True,
    hyperopt_step=100,
    arch=False,
    epochs=1000,
    X_val=False,
    y_val=False,
    input_model=False,
    allarch=True,
):
    np.random.seed(1234)
    if arch is False:
        if allarch:
            architectures = [
                (128, 128),
                (256, 256),
                (512, 512),
                (128, 128, 128),
                (256, 256, 256),
                (512, 512, 512),
            ]
        else:
            architectures = [
                (512, 512),
            ]
    else:
        architectures = [arch]
    bzs = [16, 32, 64, 128, 256]
    ress = [True, False]
    bypasses = [True, False]
    if input_model is True:
        space = {
            "lr": hp.uniform("lr", 1e-5, 1e-2),
            "batch_size": hp.choice("batch_size", bzs),
            "beta_1": hp.uniform("beta_1", 0.75, 0.99),
            "decay": hp.loguniform("decay", np.log(1e-5), np.log(1e-1)),
            "amsgrad": True,
            "patience": 100,
        }
    else:
        space = {
            "lr": hp.uniform("lr", 1e-5, 1e-3),
            "drop_rate": hp.uniform("drop_rate", 0, 0.5),
            "reg": hp.loguniform("reg", np.log(1e-5), np.log(1e0)),
            "batch_size": hp.choice("batch_size", bzs),
            "hidden_size": hp.choice("hidden_size", architectures),
            "beta_1": hp.uniform("beta_1", 0.80, 0.99),
            "decay": hp.loguniform("decay", np.log(1e-5), np.log(1e0)),
            "res": hp.choice("res", ress),
            "bypass": hp.choice("bypass", bypasses),
            "amsgrad": True,
            "patience": 100,
        }
    objective_func = partial(
        train_model_hyperopt,
        X=X,
        y=y,
        lname=lname,
        regression=regression,
        epochs=epochs,
        X_val=X_val,
        y_val=y_val,
        input_model=input_model,
    )
    trials = Trials()
    best_params = fmin(
        objective_func,
        space,
        algo=tpe.suggest,
        trials=trials,
        max_evals=hyperopt_step,
        rstate=np.random.RandomState(0),
    )
    best_params.update(
        {
            "hidden_size": architectures[best_params["hidden_size"]],
            "batch_size": bzs[best_params["batch_size"]],
            "res": ress[best_params["res"]],
            "bypass": bypasses[best_params["bypass"]],
            "amsgrad": True,
            "patience": 100,
        }
    )
    res = train_model_hyperopt(
        best_params,
        X,
        y,
        lname,
        regression=regression,
        epochs=epochs,
        X_val=X_val,
        y_val=y_val,
        input_model=input_model,
    )
    best_params.update({"epochs": res["epochs"]})
    return trials, best_params


# This function is used to train the final model.
def train_one_model(
    X_train,
    X_test,
    X_val,
    y_train,
    y_test,
    y_val,
    best_params,
    lname,
    regression,
    model_name,
    y_scaler,
):
    model = build_ANN(
        best_params, input_len=X_train.shape[-1], lname=lname, regression=regression
    )
    val_data = [X_test, y_test]
    model.fit(
        X_train,
        y_train,
        epochs=best_params["epochs"],
        verbose=2,
        batch_size=best_params["batch_size"],
        validation_data=val_data,
    )
    results_train = model.evaluate(X_train, y_train)
    results_test = model.evaluate(X_test, y_test)
    results_val = model.evaluate(X_val, y_val)
    res_dict_train, res_dict_test, res_dict_val = {}, {}, {}
    for ii, key in enumerate(model.metrics_names):
        res_dict_train.update({key: results_train[ii]})
        res_dict_test.update({key: results_test[ii]})
        res_dict_val.update({key: results_val[ii]})
    if regression:
        y_train = y_scaler.inverse_transform(y_train)
        y_test = y_scaler.inverse_transform(y_test)
        y_val = y_scaler.inverse_transform(y_val)
        hat_y_train = y_scaler.inverse_transform(model.predict(X_train))
        hat_y_test = y_scaler.inverse_transform(model.predict(X_test))
        hat_y_val = y_scaler.inverse_transform(model.predict(X_val))
        res_dict_train.update({"mae_org": mean_absolute_error(y_train, hat_y_train)})
        res_dict_test.update({"mae_org": mean_absolute_error(y_test, hat_y_test)})
        res_dict_val.update({"mae_org": mean_absolute_error(y_val, hat_y_val)})

        scaled_train_mae = mean_absolute_error(y_train, hat_y_train) / (
            max(y_train) - min(y_train)
        )
        res_dict_train.update(
            {"scaled_mae_org": scaled_train_mae[0]}
        )  # scaled_train_mae is an array with one entry
        scaled_test_mae = mean_absolute_error(y_test, hat_y_test) / (
            max(y_train) - min(y_train)
        )
        res_dict_test.update(
            {"scaled_mae_org": scaled_test_mae[0]}
        )  # scaled_test_mae is an array with one entry
        scaled_val_mae = mean_absolute_error(y_val, hat_y_val) / (
            max(y_train) - min(y_train)
        )
        res_dict_val.update(
            {"scaled_mae_org": scaled_val_mae[0]}
        )  # scaled_test_mae is an array with one entry


        train_R2 = r2_score(y_train, hat_y_train)
        res_dict_train.update({"R2_org": train_R2})
        test_R2 = r2_score(y_test, hat_y_test)
        res_dict_test.update({"R2_org": test_R2})
        val_R2 = r2_score(y_val, hat_y_val)
        res_dict_val.update({"R2_org": val_R2})

    else:
        hat_y_train = model.predict(X_train)
        hat_y_test = model.predict(X_test)
    print("res_dict_train: ", res_dict_train)
    print("res_dict_test: ", res_dict_test)
    print("res_dict_val: ", res_dict_val)
    model.save(f"{model_name}.h5")

    # Writing a text file of the results
    with open(f'performance_{lname[0]}.txt', 'w') as f:
        f.write('Train metrics\n')
        f.write(f'{str(res_dict_train)}\n')
        f.write('Test metrics\n')
        f.write(f'{str(res_dict_test)}\n')
        f.write('Val metrics\n')
        f.write(f'{str(res_dict_val)}\n')


# Used to check if a column has all the same values (if so, will remove it).
def is_same(
    s,
):
    a = s.to_numpy()
    return (a[0] == a).all()


df = pd.read_csv("4-4-2-asym-AB-AC_rs_2.csv")  # Note: Change as needed.

# Removing all columns for which all data has the same values.

cols_to_remove = []

for label, content in df.iteritems():  # iterating over columns
    if is_same(content):
        cols_to_remove.append(label)

df = df.drop(
    columns=cols_to_remove
)  # removing the columns that don't give us any information

with open("dropped_cols.txt", "w") as f:
    f.write(str(cols_to_remove))

df_train = df[df["bin"] == "train"]
df_val = df[df["bin"] == "val"]
df_test = df[df["bin"] == "test"]
df_train_all = df_train.append(df_val)

fnames = [val for val in df_train.columns.values if val in descriptors]
target_property = "C2"  # Note: Can adjust this to target the other C_i
lname = [target_property]

# Normalizing with train and validation data.
_, X_test, _, y_test, x_scaler, y_scaler = normalize_data(
    df_train_all, df_test, fnames, lname, unit_trans=1, debug=True
)
_, X_val, _, y_val, _, _ = normalize_data(
    df_train_all, df_val, fnames, lname, unit_trans=1, debug=True
)
X_train_all, X_train, y_train_all, y_train, _, _ = normalize_data(
    df_train_all, df_train, fnames, lname, unit_trans=1, debug=True
)

trials, best_params = optimize(
    X_train,
    y_train,
    lname=lname,
    regression=True,
    hyperopt_step=200,
    arch=False,
    epochs=2000,
    X_val=X_val,
    y_val=y_val,
    input_model=False,
    allarch=True,
)
print("best_params: ", best_params)
with open(f"best_params_{target_property}.txt", "w") as f:
    f.write('Best hyperparameters:\n')
    f.write(f'{str(best_params)}\n')
with open(f"trial_{target_property}.pkl", "wb") as fo:
    pickle.dump(trials, fo)

train_one_model(
    X_train_all,
    X_test,
    X_val,
    y_train_all,
    y_test,
    y_val,
    best_params,
    lname,
    regression=True,
    model_name=f"quant_ANN_{target_property}",
    y_scaler=y_scaler,
)

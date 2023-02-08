from collections import namedtuple
from econml.dml import CausalForestDML
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.multioutput import MultiOutputRegressor


def run_CF_c_py(
    Y,
    T,
    X,
    W,
    model_y=GradientBoostingRegressor(
        n_estimators=100, max_depth=3, min_samples_leaf=20
    ),
    model_t=MultiOutputRegressor(
        GradientBoostingRegressor(n_estimators=100, max_depth=3, min_samples_leaf=20)
    ),
    cv=5,
    criterion="mse",
    n_estimators=1000,
    min_samples_leaf=10,
    min_impurity_decrease=0.001,
    random_state=78343,
):

    est = CausalForestDML(
        model_y=model_y,
        model_t=model_t,
        cv=int(cv),
        criterion=criterion,
        n_estimators=int(n_estimators),
        min_samples_leaf=int(min_samples_leaf),
        min_impurity_decrease=min_impurity_decrease,
        random_state=int(random_state),
    )

    est.tune(Y, T, X=X, W=W)
    est.fit(Y, T, X=X, W=W)

    return est

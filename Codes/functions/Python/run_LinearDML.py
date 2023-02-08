from econml.dml import LinearDML
from sklearn.preprocessing import PolynomialFeatures


def run_LinearDML(
    Y,
    T,
    X,
    W,
    X_test,
    model_y=GradientBoostingRegressor(
        n_estimators=100, max_depth=3, min_samples_leaf=20
    ),
    model_t=MultiOutputRegressor(
        GradientBoostingRegressor(n_estimators=100, max_depth=3, min_samples_leaf=20)
    ),
    featurizer=PolynomialFeatures(degree=2, include_bias=False),
    linear_first_stages=False,
    cv=5,
):

    est = LinearDML(
        model_y=model_y,
        model_t=model_t,
        featurizer=featurizer,
        linear_first_stages=linear_first_stages,
        cv=int(cv),
    )

    est.fit(Y, T, X=X, W=W)

    treatment_effects = est.const_marginal_effect(X_test)

    return treatment_effects

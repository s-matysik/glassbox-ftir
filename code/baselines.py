"""Black-box baselines of PROTOCOL_v1.1 Section 6, on the raw spectrum.

All four see the full 1869-point transmittance spectrum, never the ACO
attributes. Hyperparameters are tuned inside the training fold only: a
stratified inner cross-validation for SVM, random forest and PLS-DA, and a
single stratified inner split for the convolutional network, whose fit cost
makes a full inner grid intractable at twenty seeds (declared as a deviation).

Every estimator returns hard labels for the held-out rows plus the selected
hyperparameters, so the caller can record what was chosen on each fold.
"""
from __future__ import annotations

import numpy as np
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.cross_decomposition import PLSRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV, StratifiedKFold, train_test_split
from sklearn.preprocessing import LabelBinarizer, StandardScaler
from sklearn.svm import SVC

SVM_GRID = {"C": [1.0, 10.0, 100.0], "gamma": ["scale", 1e-3, 1e-4]}
RF_GRID = {"max_features": ["sqrt", 0.05]}
PLS_GRID = {"n_components": [10, 20, 30]}
CNN_LR_GRID = [1e-3, 3e-3]
CNN_EPOCHS = 40
CNN_SCREEN_EPOCHS = 20


class PLSDA(BaseEstimator, ClassifierMixin):
    """PLS regression on one-hot class indicators, arg-max decision.

    A NIPALS fit with ``n_components = m`` contains every fit with fewer
    components, so ``predict_k`` reproduces ``PLSRegression(n_components=k)``
    exactly (verified to 3e-15) without refitting. The inner search therefore
    costs one fit per inner fold instead of one per grid point.
    """

    def __init__(self, n_components: int = 10):
        self.n_components = n_components

    def fit(self, X, y):
        self._lb = LabelBinarizer().fit(y)
        self.classes_ = self._lb.classes_
        nc = int(min(self.n_components, X.shape[0] - 1, X.shape[1]))
        self.n_fitted_ = max(nc, 1)
        self._pls = PLSRegression(n_components=self.n_fitted_).fit(
            X, self._lb.transform(y).astype(float))
        return self

    def decision_k(self, X, k: int):
        p = self._pls
        k = int(min(k, self.n_fitted_))
        coef = p.x_rotations_[:, :k] @ p.y_loadings_[:, :k].T
        return ((X - p._x_mean) / p._x_std) @ coef * p._y_std + p._y_mean

    def predict_k(self, X, k: int):
        return self.classes_[np.argmax(self.decision_k(X, k), axis=1)]

    def predict(self, X):
        return self.predict_k(X, self.n_fitted_)


def _inner_cv(y, seed, n_inner=3):
    n = int(min(n_inner, np.unique(y, return_counts=True)[1].min()))
    if n < 2:
        return None
    return StratifiedKFold(n_splits=n, shuffle=True, random_state=seed)


def _grid(est, grid, Xtr, ytr, Xte, seed):
    cv = _inner_cv(ytr, seed)
    if cv is None:
        keys = list(grid)
        params = {k: grid[k][0] for k in keys}
        return est.set_params(**params).fit(Xtr, ytr).predict(Xte), params
    gs = GridSearchCV(est, grid, cv=cv, scoring="balanced_accuracy", n_jobs=1,
                      refit=True).fit(Xtr, ytr)
    return gs.predict(Xte), gs.best_params_


def svm_rbf(Xtr, ytr, Xte, seed):
    sc = StandardScaler().fit(Xtr)
    return _grid(SVC(kernel="rbf"), SVM_GRID, sc.transform(Xtr), ytr, sc.transform(Xte), seed)


def random_forest(Xtr, ytr, Xte, seed):
    est = RandomForestClassifier(n_estimators=200, random_state=seed, n_jobs=1)
    return _grid(est, RF_GRID, Xtr, ytr, Xte, seed)


def plsda(Xtr, ytr, Xte, seed):
    from sklearn.metrics import balanced_accuracy_score
    sc = StandardScaler().fit(Xtr)
    A, B = sc.transform(Xtr), sc.transform(Xte)
    grid = [k for k in PLS_GRID["n_components"] if k <= min(len(ytr) - 1, A.shape[1])]
    grid = grid or [1]
    cv = _inner_cv(ytr, seed)
    if cv is None:
        best = grid[0]
    else:
        scores = np.zeros(len(grid))
        for itr, ite in cv.split(A, ytr):
            m = PLSDA(max(grid)).fit(A[itr], ytr[itr])
            for j, k in enumerate(grid):
                scores[j] += balanced_accuracy_score(ytr[ite], m.predict_k(A[ite], k))
        best = grid[int(np.argmax(scores))]
    model = PLSDA(best).fit(A, ytr)
    return model.predict_k(B, best), {"n_components": model.n_fitted_}


def _make_cnn(n_classes):
    import torch.nn as nn
    return nn.Sequential(
        nn.Conv1d(1, 16, 15, stride=4, padding=7), nn.ReLU(), nn.MaxPool1d(2),
        nn.Conv1d(16, 32, 9, stride=2, padding=4), nn.ReLU(), nn.MaxPool1d(2),
        nn.AdaptiveMaxPool1d(4), nn.Flatten(), nn.Linear(32 * 4, n_classes))


def _fit_cnn(X, yi, n_classes, lr, epochs, seed, bs=64):
    import torch
    import torch.nn as nn
    torch.manual_seed(seed)
    net = _make_cnn(n_classes)
    opt = torch.optim.Adam(net.parameters(), lr=lr)
    loss_fn = nn.CrossEntropyLoss()
    Xt = torch.tensor(X, dtype=torch.float32).unsqueeze(1)
    yt = torch.tensor(yi, dtype=torch.long)
    gen = torch.Generator().manual_seed(seed)
    n = len(Xt)
    for _ in range(epochs):
        perm = torch.randperm(n, generator=gen)
        for i in range(0, n, bs):
            b = perm[i:i + bs]
            opt.zero_grad()
            loss_fn(net(Xt[b]), yt[b]).backward()
            opt.step()
    return net


def _cnn_predict(net, X, classes):
    import torch
    with torch.no_grad():
        logits = net(torch.tensor(X, dtype=torch.float32).unsqueeze(1))
    return classes[logits.argmax(1).numpy()]


def cnn1d(Xtr, ytr, Xte, seed):
    sc = StandardScaler().fit(Xtr)
    A, B = sc.transform(Xtr), sc.transform(Xte)
    classes = np.array(sorted(set(map(str, ytr))))
    look = {c: i for i, c in enumerate(classes)}
    yi = np.array([look[str(v)] for v in ytr])
    counts = np.unique(yi, return_counts=True)[1]
    if counts.min() >= 2 and len(CNN_LR_GRID) > 1:
        ia, ib = train_test_split(np.arange(len(yi)), test_size=0.25,
                                  stratify=yi, random_state=seed)
        best, best_lr = -1.0, CNN_LR_GRID[0]
        for lr in CNN_LR_GRID:
            net = _fit_cnn(A[ia], yi[ia], len(classes), lr, CNN_SCREEN_EPOCHS, seed)
            pred = _cnn_predict(net, A[ib], classes)
            score = float(np.mean([(pred[yi[ib] == j] == classes[j]).mean()
                                   for j in np.unique(yi[ib])]))
            if score > best:
                best, best_lr = score, lr
    else:
        best_lr = CNN_LR_GRID[0]
    net = _fit_cnn(A, yi, len(classes), best_lr, CNN_EPOCHS, seed)
    return _cnn_predict(net, B, classes), {"lr": best_lr}


METHODS = {"svm_rbf": svm_rbf, "random_forest": random_forest,
           "plsda": plsda, "cnn1d": cnn1d}

"""Epigenetic-clock machinery: Horvath-style transformation + elastic-net clocks.

The design follows Horvath (2013, Genome Biology): chronological age is
transformed with a log-linear function, CpG beta values are regressed on
transformed age with a penalized (elastic-net) linear model, and predictions
are mapped back to years with the inverse transform.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNetCV
from sklearn.model_selection import KFold


def horvath_transform(age: np.ndarray, adult_age: float = 20.0) -> np.ndarray:
    """Log transform below `adult_age`, linear above it (Horvath 2013)."""
    age = np.asarray(age, dtype=float)
    return np.where(
        age <= adult_age,
        np.log((age + 1.0) / (adult_age + 1.0)),
        (age - adult_age) / (adult_age + 1.0),
    )


def inverse_horvath_transform(t: np.ndarray, adult_age: float = 20.0) -> np.ndarray:
    t = np.asarray(t, dtype=float)
    return np.where(
        t <= 0,
        np.exp(t) * (adult_age + 1.0) - 1.0,
        t * (adult_age + 1.0) + adult_age,
    )


@dataclass
class EpigeneticClock:
    """A fitted clock: per-probe coefficients on transformed age."""

    probes: list[str]
    coef: np.ndarray
    intercept: float
    adult_age: float = 20.0

    def predict(self, beta: pd.DataFrame) -> np.ndarray:
        """Predict DNAm age for samples (rows) x probes (cols).

        Missing clock probes are silently skipped and predictions rescaled to
        the fraction of clock weight actually observed — standard practice when
        array versions differ (450k vs EPIC).
        """
        present = [p for p in self.probes if p in beta.columns]
        if not present:
            raise ValueError("No clock probes found in the beta matrix.")
        idx = [self.probes.index(p) for p in present]
        coef = self.coef[idx]
        x = beta[present].to_numpy(dtype=float)
        x = np.nan_to_num(x, nan=0.5)
        weight_present = np.sum(np.abs(coef))
        weight_total = np.sum(np.abs(self.coef)) + 1e-12
        raw = x @ coef + self.intercept * (weight_present / weight_total)
        return inverse_horvath_transform(raw, self.adult_age)

    @classmethod
    def from_coefficient_csv(cls, path: str, adult_age: float = 20.0) -> "EpigeneticClock":
        """Load published coefficients: CSV with columns probe,coef (+ optional intercept row)."""
        df = pd.read_csv(path)
        intercept = 0.0
        if "intercept" in df.columns:
            intercept = float(df["intercept"].iloc[0])
        return cls(df["probe"].tolist(), df["coef"].to_numpy(float), intercept, adult_age)


def train_clock(
    beta: pd.DataFrame,
    ages: np.ndarray,
    cv_folds: int = 5,
    seed: int = 0,
    l1_ratios: tuple[float, ...] = (0.5, 0.7, 0.9, 0.95, 0.99),
) -> EpigeneticClock:
    """Train an elastic-net DNAm-age clock with cross-validated regularization.

    Parameters
    ----------
    beta : samples x CpG probes (values in [0, 1])
    ages : chronological ages per sample
    """
    x = beta.to_numpy(dtype=float)
    x = np.nan_to_num(x, nan=0.5)
    y = horvath_transform(np.asarray(ages, dtype=float))
    cv = KFold(n_splits=cv_folds, shuffle=True, random_state=seed)
    model = ElasticNetCV(
        l1_ratio=list(l1_ratios), cv=cv, random_state=seed, max_iter=10_000, n_jobs=-1
    )
    model.fit(x, y)
    return EpigeneticClock(beta.columns.tolist(), model.coef_, float(model.intercept_))


def age_acceleration_residual(dnam_age: np.ndarray, chrono_age: np.ndarray) -> np.ndarray:
    """IEAA-style residual: DNAm age regressed on chronological age."""
    dnam_age = np.asarray(dnam_age, dtype=float)
    chrono_age = np.asarray(chrono_age, dtype=float)
    slope, intercept = np.polyfit(chrono_age, dnam_age, 1)
    return dnam_age - (slope * chrono_age + intercept)

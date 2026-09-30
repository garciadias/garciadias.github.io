"""Figures for the two validation-metric slides of the IAA-SO chemical tagging deck.

Generates, into public/presentations/iaa-so-chemical-tagging-2026/:

  silhouette_explained.png   the a(i) / b(i) geometry beside a real silhouette plot
  homogeneity_explained.png  three clusterings of one labelled set, scored

Every number printed on a figure is computed here by scikit-learn, never typed by
hand, so the slides can be checked by re-running this script. The printed summary
at the end is what the slide text must agree with.

Run (no project dependency on matplotlib; use a throwaway environment):

    uv venv /tmp/figvenv --python 3.12
    uv pip install --python /tmp/figvenv/bin/python matplotlib scikit-learn numpy
    /tmp/figvenv/bin/python scripts/figures/validation_metrics.py

Style follows the deck's existing figures: white ground, no ticks, hairline
spines, the astro-theme palette, a centred title and an italic muted caption.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import (
    completeness_score,
    homogeneity_score,
    silhouette_samples,
    silhouette_score,
    v_measure_score,
)

OUT = (
    Path(__file__).resolve().parents[2]
    / "public"
    / "presentations"
    / "iaa-so-chemical-tagging-2026"
)

# The deck's astro-theme palette (src/components/RevealDeck.vue).
INDIGO = "#3730a3"   # --accent
BLUE = "#0369a1"     # --accent-cyan
ROSE = "#be123c"     # --accent-pink
GREEN = "#15803d"    # --accent-green
AMBER = "#b45309"    # --accent-orange
INK = "#1c2333"      # --r-main-color
MUTED = "#55607a"    # --comment
HAIR = "#d5dae6"     # --line
GHOST = "#c3cbdb"    # unhighlighted data points

SEED = 11


def style_axes(ax, keep=("left", "bottom")):
    """Hairline spines, no ticks: the deck's figures carry no axis furniture."""
    for name, spine in ax.spines.items():
        if name in keep:
            spine.set_color(HAIR)
            spine.set_linewidth(0.8)
        else:
            spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])


def caption(fig, text):
    fig.text(
        0.012,
        0.022,
        text,
        fontsize=9.5,
        style="italic",
        color=MUTED,
        ha="left",
        va="bottom",
    )


def three_blobs(rng):
    """Three 2-D blobs with clearly different separations, so s(i) varies."""
    a = rng.normal([0.0, 0.0], [0.62, 0.62], size=(70, 2))
    b = rng.normal([3.4, 0.35], [0.62, 0.62], size=(70, 2))
    c = rng.normal([1.7, 3.05], [0.70, 0.70], size=(70, 2))
    X = np.vstack([a, b, c])
    y = np.repeat([0, 1, 2], 70)
    return X, y


# --------------------------------------------------------------------------- #
# Figure 1 - silhouette                                                        #
# --------------------------------------------------------------------------- #
def silhouette_figure():
    rng = np.random.default_rng(SEED)
    X, y = three_blobs(rng)
    colours = [BLUE, AMBER, GREEN]

    sil = silhouette_samples(X, y)
    mean_s = silhouette_score(X, y)

    # Pick a well-placed point in cluster 0 to carry the a(i) / b(i) annotation.
    own = np.where(y == 0)[0]
    i = own[np.argmin(np.abs(sil[own] - np.percentile(sil[own], 62)))]

    d = np.linalg.norm(X - X[i], axis=1)
    a_i = d[(y == 0) & (np.arange(len(X)) != i)].mean()
    b_by_cluster = {k: d[y == k].mean() for k in (1, 2)}
    nearest = min(b_by_cluster, key=lambda k: b_by_cluster[k])
    b_i = b_by_cluster[nearest]
    s_i = (b_i - a_i) / max(a_i, b_i)
    assert np.isclose(s_i, sil[i], atol=1e-9), (s_i, sil[i])

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(13.4, 4.9), gridspec_kw={"width_ratios": [1.06, 1.0]}
    )
    fig.patch.set_facecolor("white")

    # -- left: the geometry ------------------------------------------------- #
    for k in (0, 1, 2):
        ax1.scatter(
            *X[y == k].T, s=17, color=GHOST, alpha=0.75, linewidths=0, zorder=1
        )

    # Spokes from i: to its own cluster (a) and to the nearest other one (b).
    for p in X[(y == 0) & (np.arange(len(X)) != i)]:
        ax1.plot(
            [X[i, 0], p[0]], [X[i, 1], p[1]],
            color=BLUE, alpha=0.20, linewidth=0.55, zorder=2,
        )
    for p in X[y == nearest]:
        ax1.plot(
            [X[i, 0], p[0]], [X[i, 1], p[1]],
            color=ROSE, alpha=0.13, linewidth=0.55, zorder=2,
        )

    ax1.scatter(
        *X[i], s=300, marker="*", color=INDIGO,
        edgecolors="white", linewidths=1.1, zorder=6,
    )

    # Mean-distance rings: a(i) and b(i) drawn to scale about i.
    for radius, colour, label in ((a_i, BLUE, "a"), (b_i, ROSE, "b")):
        ax1.add_patch(
            plt.Circle(
                X[i], radius, fill=False, color=colour,
                linewidth=1.5, linestyle=(0, (5, 3)), alpha=0.85, zorder=5,
            )
        )
        ax1.annotate(
            f"{label}(i) = {radius:.2f}",
            xy=(X[i, 0], X[i, 1] + radius),
            xytext=(X[i, 0] - 0.1, X[i, 1] + radius + 0.22),
            fontsize=12, color=colour, fontweight="bold", ha="center",
        )

    ax1.text(
        0.5, -0.085,
        f"s(i) = (b − a) / max(a, b) = ({b_i:.2f} − {a_i:.2f}) / {b_i:.2f} = {s_i:.2f}",
        transform=ax1.transAxes, fontsize=13, color=INK,
        ha="center", va="top", fontweight="bold",
    )
    ax1.text(
        0.017, 0.965, "own cluster", transform=ax1.transAxes,
        fontsize=11.5, color=BLUE, fontweight="bold", va="top",
    )
    ax1.text(
        0.017, 0.898, "nearest other cluster", transform=ax1.transAxes,
        fontsize=11.5, color=ROSE, fontweight="bold", va="top",
    )
    ax1.set_title(
        "One point: how close is home, how close is the next street?",
        fontsize=13.5, color=INK, pad=12,
    )
    ax1.set_aspect("equal")
    ax1.margins(0.13)
    style_axes(ax1)

    # -- right: the silhouette plot ----------------------------------------- #
    offset = 0
    gap = 9
    for k in (0, 1, 2):
        vals = np.sort(sil[y == k])
        ypos = np.arange(offset, offset + len(vals))
        ax2.barh(ypos, vals, height=1.0, color=colours[k], alpha=0.82, linewidth=0)
        ax2.text(
            -0.045, offset + len(vals) / 2, f"cluster {k}",
            fontsize=11, color=colours[k], ha="right", va="center", fontweight="bold",
        )
        ax2.text(
            0.012, offset + len(vals) / 2, f"mean {vals.mean():.2f}",
            fontsize=10.5, color="white", ha="left", va="center", fontweight="bold",
        )
        offset += len(vals) + gap

    ax2.axvline(mean_s, color=INDIGO, linestyle="--", linewidth=1.6, zorder=5)
    ax2.text(
        mean_s + 0.012, offset - 4,
        f"overall mean {mean_s:.2f}",
        fontsize=11.5, color=INDIGO, fontweight="bold", va="top",
    )
    ax2.axvline(0, color=MUTED, linewidth=0.9)

    # Kaufman & Rousseeuw's published reading of the coefficient.
    for thresh, label in ((0.26, "weak"), (0.51, "reasonable"), (0.71, "strong")):
        ax2.axvline(thresh, color=HAIR, linewidth=1.0, zorder=0)
        ax2.text(
            thresh, -gap * 1.5, label, fontsize=9.5,
            color=MUTED, ha="center", va="top",
        )

    ax2.set_xlim(-0.12, 1.0)
    ax2.set_ylim(-gap * 3.2, offset)
    ax2.set_title(
        "Every point, sorted inside its cluster: the silhouette plot",
        fontsize=13.5, color=INK, pad=12,
    )
    for name, spine in ax2.spines.items():
        spine.set_visible(False)
    ax2.set_yticks([])
    ax2.set_xticks([0.0, 0.26, 0.51, 0.71, 1.0])
    ax2.tick_params(axis="x", colors=MUTED, labelsize=10, length=0, pad=16)

    caption(
        fig,
        "s(i) runs from −1 to 1: 1 is comfortably home, 0 is on the border, "
        "below 0 means the point sits closer to another cluster than to its own. "
        "Thresholds: Kaufman & Rousseeuw (1990).",
    )
    fig.tight_layout(rect=(0, 0.055, 1, 1))
    fig.savefig(OUT / "silhouette_explained.png", dpi=170, facecolor="white")
    plt.close(fig)

    return {
        "a_i": a_i, "b_i": b_i, "s_i": s_i, "mean_s": mean_s,
        "per_cluster": {k: sil[y == k].mean() for k in (0, 1, 2)},
        "frac_negative": float((sil < 0).mean()),
    }


# --------------------------------------------------------------------------- #
# Figure 2 - homogeneity                                                       #
# --------------------------------------------------------------------------- #
def homogeneity_figure():
    """One labelled set, three clusterings: the point is that h alone is gameable."""
    rng = np.random.default_rng(SEED)
    X, truth = three_blobs(rng)
    n = len(X)

    # (1) recovers the truth, (2) shatters every class, (3) merges two classes.
    perfect = truth.copy()

    shattered = np.empty(n, dtype=int)
    for k in range(3):
        idx = np.where(truth == k)[0]
        # split each true class by angle about its own centroid: 3 wedges
        ang = np.arctan2(*(X[idx] - X[idx].mean(0)).T[::-1])
        shattered[idx] = k * 3 + np.digitize(ang, [-np.pi / 3, np.pi / 3])

    merged = np.where(truth == 2, 1, truth)  # classes 1 and 2 collapsed into one

    panels = [
        ("Recovers the truth", perfect, "3 clusters, 3 classes"),
        ("Shattered: 9 clusters", shattered, "every class cut into 3"),
        ("Merged: 2 clusters", merged, "two classes in one cluster"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(12.0, 5.5))
    fig.patch.set_facecolor("white")
    marks = ["o", "s", "^"]
    scores = {}

    for ax, (title, pred, sub) in zip(axes, panels):
        h = homogeneity_score(truth, pred)
        c = completeness_score(truth, pred)
        v = v_measure_score(truth, pred)
        scores[title] = (h, c, v)

        palette = plt.get_cmap("tab20")(np.linspace(0, 1, 20))
        for cl in np.unique(pred):
            sel = pred == cl
            # marker = the TRUE class, colour = the cluster the algorithm gave
            for t in np.unique(truth[sel]):
                m = sel & (truth == t)
                ax.scatter(
                    *X[m].T, s=30, marker=marks[t],
                    color=palette[(cl * 3) % 20],
                    alpha=0.9, linewidths=0.4, edgecolors="white",
                )

        ax.set_title(
            f"{title}\n{sub}", fontsize=13.5, color=INK, pad=10, linespacing=1.55,
        )

        hcol = ROSE if h > 0.995 and c < 0.995 else INK
        ccol = ROSE if c > 0.995 and h < 0.995 else INK
        ax.text(
            0.5, -0.115,
            f"homogeneity {h:.2f}",
            transform=ax.transAxes, fontsize=13.5, color=hcol,
            ha="center", va="top", fontweight="bold",
        )
        ax.text(
            0.5, -0.225,
            f"completeness {c:.2f}      V {v:.2f}",
            transform=ax.transAxes, fontsize=11.5, color=ccol,
            ha="center", va="top",
        )
        ax.set_aspect("equal")
        ax.margins(0.14)
        style_axes(ax)

    axes[1].text(
        0.5, -0.325, "perfect purity, bought by cutting the truth up",
        transform=axes[1].transAxes, fontsize=11, color=ROSE,
        ha="center", va="top", style="italic",
    )
    axes[2].text(
        0.5, -0.325, "perfect coverage, bought by lumping it together",
        transform=axes[2].transAxes, fontsize=11, color=ROSE,
        ha="center", va="top", style="italic",
    )

    caption(
        fig,
        "Marker shape = the true class, colour = the cluster the algorithm produced.\n"
        "h = 1 − H(C|K)/H(C) (Rosenberg & Hirschberg 2007); V is the harmonic mean "
        "of homogeneity and completeness.",
    )
    # Explicit margins, not tight_layout: the score/annotation text below each
    # panel is placed at negative axes coordinates, which tight_layout cannot
    # see, so it kept clipping the two-line titles and colliding the notes
    # with the caption.
    fig.subplots_adjust(left=0.02, right=0.98, top=0.88, bottom=0.36, wspace=0.12)
    fig.savefig(OUT / "homogeneity_explained.png", dpi=170, facecolor="white")
    plt.close(fig)

    return scores


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)

    s = silhouette_figure()
    print("silhouette_explained.png")
    print(f"  a(i)={s['a_i']:.3f}  b(i)={s['b_i']:.3f}  s(i)={s['s_i']:.3f}")
    print(f"  overall mean silhouette = {s['mean_s']:.3f}")
    for k, v in s["per_cluster"].items():
        print(f"  cluster {k} mean = {v:.3f}")
    print(f"  fraction of points with s < 0 = {s['frac_negative']:.3f}")

    h = homogeneity_figure()
    print("homogeneity_explained.png")
    for title, (hh, cc, vv) in h.items():
        print(f"  {title:24} h={hh:.3f}  c={cc:.3f}  V={vv:.3f}")

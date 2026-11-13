"""Joint SAMap UMAP of panels a-c of the SAMap supplementary figure: the nuclei drawn, their order, and the frame.

Imported by the notebooks of panels a, b and c, which add their own color layer and key. load_umap_cells() reads
the three tables it is given; nothing is written.
"""
import matplotlib as mpl
import pandas as pd
from matplotlib.lines import Line2D

from sfig_samap_style import PSEUDOTIME_CMAP, PSEUDOTIME_NORM, new_panel, panel_title

DOT = 0.22                 # pt^2: about 90,000 nuclei share each 57 mm panel
OTHER_SPECIES = "#3a3a3a"  # nuclei of the species whose values a panel does not show (as Fig. 6, MOUSE_CONTEXT)
AXIS_GREY = "#555555"      # UMAP axis arrows and their labels


def load_umap_cells(joint_path, mouse_cells_path, human_cells_path) -> dict:
    """Nuclei drawn on the UMAPs, each species in a seeded random order, and the frame around them.

    Every nucleus is drawn. The random order keeps early or late nuclei from always sitting on top; the seeds and the
    joint object's nucleus order fix it. Returns {"mouse", "human"} DataFrames (umap1, umap2 and the values), "xlim",
    "ylim" and "aspect".
    """
    umap = (
        pd.read_parquet(joint_path, columns=["species", "X_umap__0", "X_umap__1"])
        .rename(columns={"X_umap__0": "umap1", "X_umap__1": "umap2"})
    )
    umap["species"] = umap["species"].astype(str)
    if not umap.index.is_unique:
        raise ValueError(f"{joint_path} has duplicated nucleus names")

    mouse = pd.read_csv(mouse_cells_path, index_col="cell")
    human = pd.read_csv(human_cells_path, index_col="cell")

    # Stored global rank (average rank / N over the non-apoptotic mouse nuclei; NaN for apoptotic ones).
    mouse_values = mouse[["cell_type", "cell_subtype"]].assign(pt_rank=mouse["pt_harmonized_rank"])

    # Human pseudotime (pt_human) and the transfer of record (pt_ranknorm), one row per human nucleus.
    human_values = human[["pt_transferred", "pt_human", "cell_type"]].rename(columns={"pt_transferred": "pt_ranknorm"})

    mouse_drawn = umap[umap.species == "mm"].join(mouse_values[["pt_rank", "cell_type", "cell_subtype"]], how="left")
    human_drawn = umap[umap.species == "hs"].join(human_values[["pt_ranknorm", "pt_human", "cell_type"]], how="left")

    # Frame: every nucleus plus 2% of the range on each side.
    xy = umap[["umap1", "umap2"]].values
    pad = 0.02 * (xy.max(0) - xy.min(0))
    xlim = (xy[:, 0].min() - pad[0], xy[:, 0].max() + pad[0])
    ylim = (xy[:, 1].min() - pad[1], xy[:, 1].max() + pad[1])

    return {
        "mouse": mouse_drawn.sample(frac=1, random_state=1),
        "human": human_drawn.sample(frac=1, random_state=3),
        "xlim": xlim,
        "ylim": ylim,
        "aspect": (xlim[1] - xlim[0]) / (ylim[1] - ylim[0]),
    }


def umap_panel(letter: str, title: str, cells: dict, draw, key=None, colorbar_title=None):
    """One UMAP panel: equal-aspect axes without a box, short axis arrows, and an optional key and colorbar.

    The key and the colorbar sit in a band between the title and the UMAP, so the UMAP takes the panel's width and
    has the same size in every panel of the row. draw(ax) adds the panel's own scatter layers. key is a list of
    (color, label) for the grey nuclei; the colorbar is the fixed 0-1 pseudotime scale.
    """
    fig, width_mm, height_mm = new_panel(letter)
    title_mm, band_mm, bottom_mm, left_mm, right_mm = 10.0, 11.0, 4.0, 6.0, 1.5
    # The largest equal-aspect axes that fit under the band, at the bottom left.
    axes_h_mm = height_mm - title_mm - band_mm - bottom_mm
    axes_w_mm = axes_h_mm * cells["aspect"]
    if axes_w_mm > width_mm - left_mm - right_mm:  # too wide: shrink, keeping the aspect
        axes_w_mm = width_mm - left_mm - right_mm
        axes_h_mm = axes_w_mm / cells["aspect"]
    ax = fig.add_axes([left_mm / width_mm, bottom_mm / height_mm, axes_w_mm / width_mm, axes_h_mm / height_mm])
    xlim, ylim = cells["xlim"], cells["ylim"]
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    panel_title(fig, width_mm, height_mm, title)

    draw(ax)

    # UMAP coordinates have no units, so two short arrows replace a box and ticks.
    dx, dy = xlim[1] - xlim[0], ylim[1] - ylim[0]
    length = 0.16 * dx
    arrow = dict(
        arrowprops=dict(arrowstyle="-|>", lw=0.8, color=AXIS_GREY, mutation_scale=6,
                        shrinkA=0, shrinkB=0, capstyle="projecting", joinstyle="miter"),
        annotation_clip=False,
    )
    ax.annotate("", xy=(xlim[0] + length, ylim[0]), xytext=(xlim[0], ylim[0]), **arrow)
    ax.annotate("", xy=(xlim[0], ylim[0] + length), xytext=(xlim[0], ylim[0]), **arrow)
    ax.text(xlim[0] + 0.02 * dx, ylim[0] - 0.012 * dy, "UMAP 1", ha="left", va="top", fontsize=7, color=AXIS_GREY)
    ax.text(xlim[0] - 0.03 * dx, ylim[0] + 0.02 * dy, "UMAP 2", ha="right", va="bottom", rotation=90,
            fontsize=7, color=AXIS_GREY)

    band_top = 1 - title_mm / height_mm
    if colorbar_title:
        # A short horizontal bar at the band's right end: its title above, its ticks below.
        bar_w_mm, bar_h_mm = 15.0, 1.6
        cax = fig.add_axes([(width_mm - right_mm - 1.0 - bar_w_mm) / width_mm, band_top - (4.0 + bar_h_mm) / height_mm,
                            bar_w_mm / width_mm, bar_h_mm / height_mm])
        bar = fig.colorbar(mpl.cm.ScalarMappable(norm=PSEUDOTIME_NORM, cmap=PSEUDOTIME_CMAP), cax=cax,
                           orientation="horizontal", ticks=[0, 0.5, 1])
        bar.ax.set_xticklabels(["0", "0.5", "1"])
        bar.outline.set_visible(False)
        bar.ax.tick_params(length=1.5, pad=1.5)
        cax.set_title(colorbar_title, fontsize=8, loc="left", pad=2.5)
    if key:
        kax = fig.add_axes([left_mm / width_mm, band_top - band_mm / height_mm,
                            (width_mm - left_mm - right_mm) / width_mm, band_mm / height_mm])
        kax.axis("off")
        kax.legend(
            handles=[Line2D([], [], ls="none", marker="o", ms=3.2, color=color, label=label) for color, label in key],
            loc="upper left",
            frameon=False,
            fontsize=8,
            handletextpad=0.2,
            labelspacing=0.5,
            borderaxespad=0,
            handlelength=0.9,
        )
    return fig

#!/usr/bin/env python3
"""Overview of the seeded multiscale potential and its BDM vector field."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import plot_decompositions as base


def scalar_p2_image(q, n, tri, te, resolution=600):
    coords = (np.arange(resolution)+.5)/resolution
    xx, yy = np.meshgrid(coords, coords)
    i = np.minimum((xx*n).astype(int), n-1)
    j = np.minimum((yy*n).astype(int), n-1)
    fx, fy = xx*n-i, yy*n-j
    upper = fy > fx
    tid = 2*(j*n+i)+upper.astype(int)
    lam0 = np.where(upper, 1-fy, 1-fx)
    lam1 = np.where(upper, fx, fx-fy)
    lam2 = np.where(upper, fy-fx, fy)
    lam = (lam0, lam1, lam2)
    ids = np.concatenate((tri, len(set(tri.ravel()))+te), axis=1)
    vals = q[ids[tid]]
    shape = [a*(2*a-1) for a in lam] + [4*lam[0]*lam[1],
            4*lam[1]*lam[2], 4*lam[2]*lam[0]]
    return sum(vals[..., k]*shape[k] for k in range(6))


def main():
    xy, tri, edges, te, _ = base.mesh(10)
    normal, td, basis = base.bdm_data(xy, tri, edges, te)
    smap, qmap = base.prolongations(xy, edges, normal)
    fields, components = base.source_fields(xy, edges, normal, qmap,
                                            return_components=True)
    smooth = smap @ components["noisy_smooth_coefficients"]
    noise = qmap @ components["noisy_potential_coefficients"]
    source = fields["multiscale_noisy_potential"]
    qimage = scalar_p2_image(components["noisy_potential_coefficients"],
                             10, tri, te)
    smooth_im, _, _, _ = base.evaluate_bdm_image(smooth, 10, td, basis)
    noise_im, _, _, _ = base.evaluate_bdm_image(noise, 10, td, basis)
    source_im, _, _, _ = base.evaluate_bdm_image(source, 10, td, basis)
    panels = [qimage, smooth_im[0], smooth_im[1], source_im[0],
              source_im[1], noise_im[0], noise_im[1]]
    titles = ["Noisy P2 potential", "Smooth vector · x", "Smooth vector · y",
              "BDM sum · x", "BDM sum · y", "Rotgrad noise · x",
              "Rotgrad noise · y"]
    scales = [np.max(np.abs(qimage)), np.max(np.abs(smooth_im)),
              np.max(np.abs(source_im)), np.max(np.abs(noise_im))]
    plt.rcParams.update({"font.size": 11, "font.family": "DejaVu Sans"})
    fig, axes = plt.subplots(2, 4, figsize=(16, 9), dpi=120)
    axes = axes.ravel()
    for k, (arr, title) in enumerate(zip(panels, titles)):
        scale = scales[0 if k == 0 else 1 if k in (1, 2) else
                       2 if k in (3, 4) else 3]
        axes[k].imshow(arr, origin="lower", extent=(0, 1, 0, 1),
                       cmap="RdBu_r", vmin=-scale, vmax=scale,
                       interpolation="nearest")
        axes[k].set_title(title, pad=5)
        axes[k].set_xticks([]); axes[k].set_yticks([])
        axes[k].set_aspect("equal")
    axes[7].axis("off")
    axes[7].text(.04, .96,
        "Seed: 20260919\nFour Fourier bands: 1–2, 3–4,\n5–6, 7–9 cycles/unit\n"
        "Eight random modes per band\n"
        f"Smooth H¹ norm: {components['smooth_h1_norm']:.4g}\n"
        f"Potential H¹ norm: {components['potential_h1_norm']:.4g}\n"
        "Continuous P2 potential\nBDM1 vector on 10×10 cells",
        ha="left", va="top", transform=axes[7].transAxes,
        fontsize=13, linespacing=1.6)
    axes[7].text(.04, .08,
        f"Color ranges: q ±{scales[0]:.3g}; smooth ±{scales[1]:.3g};\n"
        f"BDM sum ±{scales[2]:.3g}; rotgrad ±{scales[3]:.3g}",
        ha="left", va="bottom", transform=axes[7].transAxes, fontsize=11)
    fig.suptitle("Multiscale noisy potential · source construction", fontsize=19, y=.986)
    fig.subplots_adjust(left=.018, right=.985, bottom=.035, top=.925,
                        wspace=.065, hspace=.16)
    path = Path(__file__).parent / "output" / "multiscale_noisy_potential_source.png"
    path.parent.mkdir(exist_ok=True)
    fig.savefig(path, dpi=120)
    plt.close(fig)
    print(path)


if __name__ == "__main__":
    main()

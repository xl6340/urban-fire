#!/usr/bin/env python3
"""
Redraw Fig 1A map (CalFire WUI vs wildland perimeters).

- No green California / landcover background (reviewer request)
- Light-gray footprint = WUI2020 ∪ MODIS MCD12Q1 LC_Type2 (2020) urban (class 13)
"""

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from rasterio import features
from rasterio.mask import mask as rio_mask
from rasterio.warp import Resampling, calculate_default_transform, reproject

BASE = Path(__file__).resolve().parents[2]
LC_PATH = BASE / "dataSrc/landcover/data/MCD12Q1.061_LC_Type2_doy2020001000000_aid0001.tif"
# UMD LC_Type2: 13 = Urban and Built-up Lands
URBAN_CLASS = 13


def urban_mask_in_crs(lc_path, ca_gdf, dst_crs):
    """Return (urban_bool array, transform, extent) in dst_crs, clipped to CA."""
    ca_dst = ca_gdf.to_crs(dst_crs)
    with rasterio.open(lc_path) as src:
        # reproject full LC to dst_crs at ~500 m
        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds, resolution=500
        )
        dst = np.zeros((height, width), dtype=np.uint8)
        reproject(
            source=rasterio.band(src, 1),
            destination=dst,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=transform,
            dst_crs=dst_crs,
            resampling=Resampling.nearest,
        )

        # clip to CA polygon via rasterize
        ca_mask = features.geometry_mask(
            [ca_dst.union_all()],
            out_shape=(height, width),
            transform=transform,
            invert=True,
        )
        urban = (dst == URBAN_CLASS) & ca_mask

    # extent for imshow: [left, right, bottom, top]
    left = transform.c
    top = transform.f
    right = left + width * transform.a
    bottom = top + height * transform.e
    extent = [left, right, bottom, top]
    return urban, extent


def main():
    ca = gpd.read_file(BASE / "dataPrc/boundary/CA_State.shp")
    wui = gpd.read_file(BASE / "dataPrc/WUI/WUI2020.shp").to_crs(ca.crs)
    wild = gpd.read_file(BASE / "dataPrc/firePrmt/fires/Wildland.shp").to_crs(ca.crs)
    wui_fire = gpd.read_file(BASE / "dataPrc/firePrmt/fires/Urban-edge.shp").to_crs(ca.crs)

    urban, extent = urban_mask_in_crs(LC_PATH, ca, ca.crs)
    print(f"Urban pixels (MCD12Q1 2020 class {URBAN_CLASS}): {int(urban.sum()):,}")

    c_wui_fire = np.array([216, 118, 89]) / 255
    c_wild = np.array([41, 157, 143]) / 255
    c_wui_area = np.array([214, 214, 214]) / 255
    c_outline = np.array([0.35, 0.35, 0.35])

    fig, ax = plt.subplots(figsize=(5.2, 6.2), facecolor="white")

    # CA outline only (no green fill)
    ca.boundary.plot(ax=ax, color=c_outline, linewidth=0.8, zorder=1)

    # Urban (MODIS) then WUI — same light gray ("WUI & urban")
    rgba = np.zeros((*urban.shape, 4), dtype=float)
    rgba[urban, :3] = c_wui_area
    rgba[urban, 3] = 1.0
    # aspect='auto' so imshow does not override map limits / equal aspect
    ax.imshow(
        rgba,
        extent=extent,
        origin="upper",
        zorder=2,
        interpolation="nearest",
        aspect="auto",
    )

    wui.plot(ax=ax, color=c_wui_area, edgecolor="none", linewidth=0, zorder=3)

    # Fires on top (clip_on=False: northern perimeters may extend past CA)
    wild.plot(ax=ax, color=c_wild, edgecolor="none", linewidth=0, zorder=4)
    wui_fire.plot(ax=ax, color=c_wui_fire, edgecolor="none", linewidth=0, zorder=5)
    for coll in ax.collections:
        coll.set_clip_on(False)

    ax.set_axis_off()
    # Fit extent to CA + all fire perimeters so northern events aren't clipped
    bounds = np.vstack(
        [ca.total_bounds, wild.total_bounds, wui_fire.total_bounds]
    )
    minx, miny = bounds[:, 0].min(), bounds[:, 1].min()
    maxx, maxy = bounds[:, 2].max(), bounds[:, 3].max()
    pad_x = (maxx - minx) * 0.03
    pad_y = (maxy - miny) * 0.03
    # Extra headroom on the north (largest wildland events sit above CA outline)
    ax.set_xlim(minx - pad_x, maxx + pad_x)
    ax.set_ylim(miny - pad_y, maxy + pad_y * 1.5)
    ax.set_aspect("equal", adjustable="box")

    handles = [
        Patch(facecolor=c_wui_fire, edgecolor="none", label="WUI fire"),
        Patch(facecolor=c_wild, edgecolor="none", label="Wildland fire"),
        Patch(facecolor=c_wui_area, edgecolor="0.7", label="WUI & urban"),
    ]
    ax.legend(
        handles=handles,
        loc="lower left",
        frameon=False,
        fontsize=12,
        handlelength=1.0,
        handleheight=1.0,
        borderaxespad=0.4,
    )

    fig.tight_layout(pad=0.3)
    # Re-apply limits after tight_layout (can reset view)
    ax.set_xlim(minx - pad_x, maxx + pad_x)
    ax.set_ylim(miny - pad_y, maxy + pad_y * 1.5)
    for pth in [
        BASE / "Fig/Fig1A.png",
        BASE / "Fig/Fig1A.pdf",
        BASE / "Fig/map/Fig1A_noGreen.png",
        BASE / "dataPrc/Fig1A.png",
    ]:
        fig.savefig(pth, dpi=400, bbox_inches="tight", facecolor="white", pad_inches=0.15)
        print("Saved", pth)
    plt.close(fig)


if __name__ == "__main__":
    main()

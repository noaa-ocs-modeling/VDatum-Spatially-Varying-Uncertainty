#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Plot shortest water-path distance PNGs from station*.nc files over an ADCIRC fort.14 mesh.


1) Usage examples:

# All stations, with custom MATLAB colormap (variable 'cm' in the .mat file)+ HTML
python S1_Wdist_Pacific-plot.py --mesh fort.14 --nc-dir . --pattern "station*.nc" --outdir figs \
  --cmap /full/path/to/cm_dist2.mat --make-html\
  --pac-center 180 --lon-min 130 --lon-max 230 --lat-min -25 --lat-max 55 \
  --max-val 7500 --levels-step 100


# Single station (e.g., station17.nc)
python S1_Wdist_Pacific-plot.py --mesh fort.14 --nc-dir . --pattern "station17.nc" --outdir figs \
  --station-id 17 --cmap /full/path/to/cm_dist2.mat \
  --pac-center 180 --lon-min 130 --lon-max 230 --lat-min -25 --lat-max 55 \
  --max-val 7500 --levels-step 100
  
 2)gageinfo.csv   text file with columns: station_id  node_id
 
 3) Here’s the full list you need to install into your environment:
Standard library (comes with Python, no need to install)
      Os, re, argparse, glob

External packages (must install with pip or conda)
      Numpy, pandas, scipy, netCDF4, matplotlib

conda install numpy pandas scipy netcdf4 matplotlib 

"""


import os
import re
import glob
import argparse
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation


# ------------------------ helpers ------------------------
def pacific_wrap(lon_deg, center=180.0):
    """
    Wrap longitudes so that 'center' is in the middle (e.g., 180 or 150).
    Output range: [center-180, center+180).
    """
    lon = np.asarray(lon_deg, dtype=float)
    wrapped = (lon - center + 180.0) % 360.0 - 180.0 + center
    return wrapped


def get_station_id_from_filename(fn):
    """
    Extract the digits after 'station' in a name like 'station123.nc'.
    Falls back to the filename stem if not found.
    """
    base = os.path.basename(fn)
    m = re.search(r"station(\d+)\.nc$", base, re.IGNORECASE)
    return m.group(1) if m else os.path.splitext(base)[0]


def read_f14_elems_first(path):
    """
    Read fort.14 (triangles only). Returns:
      node_ids, x, y, z, tris(0-based), sorted_ids, order, title
    """
    with open(path, "rt", newline="") as f:
        title = f.readline().strip()
        nums  = f.readline().strip().split()
        numele, numnod = int(nums[0]), int(nums[1])


        node_ids = np.empty(numnod, dtype=np.int64)
        x = np.empty(numnod, dtype=np.float64)
        y = np.empty(numnod, dtype=np.float64)
        z = np.empty(numnod, dtype=np.float64)
        for i in range(numnod):
            parts = f.readline().split()
            node_ids[i] = int(float(parts[0]))
            x[i], y[i], z[i] = float(parts[1]), float(parts[2]), float(parts[3])


        tri_ids = np.empty((numele, 3), dtype=np.int64)
        k = 0
        while k < numele:
            parts = f.readline().split()
            if not parts:
                continue
            if len(parts) == 4:
                n1, n2, n3 = int(parts[1]), int(parts[2]), int(parts[3])
            else:
                n1, n2, n3 = int(parts[-3]), int(parts[-2]), int(parts[-1])
            tri_ids[k] = (n1, n2, n3)
            k += 1


    # Map 1-based node ids in triangles -> 0-based indices in x/y/z via sorting
    order      = np.argsort(node_ids)
    sorted_ids = node_ids[order]
    pos        = np.searchsorted(sorted_ids, tri_ids)
    ok         = (pos >= 0) & (pos < sorted_ids.size) & (sorted_ids[np.clip(pos,0,sorted_ids.size-1)] == tri_ids)
    tris       = order[pos[np.all(ok, axis=1)]]
    valid_tri  = np.all((tris >= 0) & (tris < node_ids.size), axis=1)
    return node_ids, x, y, z, tris[valid_tri], sorted_ids, order, title


def ensure_cmap():
    """
    Use user's custom_cmap if defined at runtime; else fall back to viridis.
    """
    try:
        return custom_cmap  # set dynamically when --cmap is used
    except NameError:
        return plt.cm.viridis


# ------------------------ plotting core ------------------------
def plot_one_station(tri_full, x_wrapped, y, dis, node0_id, sorted_ids, order,
                     station_label, out_png,
                     lon_min, lon_max, lat_min, lat_max,
                     max_val, levels_step):
    # Convert sentinel to NaN
    dis = np.asarray(dis, dtype=float)
    dis = np.where(dis >= 9.999e5, np.nan, dis)


    # Build an all-tri triangulation to compute a mask for NaNs
    tri_all = Triangulation(x_wrapped, y, triangles=tri_full.triangles)
    mask_nan = np.any(~np.isfinite(dis[tri_all.triangles]), axis=1)
    tri_all.set_mask(mask_nan)


    # Crop to bbox by keeping triangles with any vertex inside the box
    tp = tri_all.triangles
    xv, yv = x_wrapped[tp], y[tp]
    in_box = (xv >= lon_min) & (xv <= lon_max) & (yv >= lat_min) & (yv <= lat_max)
    keep_tri = np.any(in_box, axis=1)
    tri = Triangulation(x_wrapped, y, triangles=tp[keep_tri])
    tri.set_mask(np.any(~np.isfinite(dis[tri.triangles]), axis=1))  # re-mask NaNs in cropped mesh


    # Levels & colormap
    levels = np.arange(0, max_val + levels_step, levels_step)
    cmap = ensure_cmap()


    # Find source node index
    p0 = np.searchsorted(sorted_ids, int(node0_id))
    node0_idx = order[p0] if (p0 < sorted_ids.size and sorted_ids[p0] == int(node0_id)) else None


    # Plot
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_facecolor("black")
    cf = ax.tricontourf(tri, dis, levels=levels, cmap=cmap, vmin=0, vmax=max_val)
    cbar = plt.colorbar(cf, ax=ax)
    cbar.set_label("Water-path distance (km)")
    cbar.set_ticks(np.arange(0, max_val + 1000, 1000))


    if node0_idx is not None:
        gx, gy = x_wrapped[node0_idx], y[node0_idx]
        ax.plot(gx, gy, marker="+", markersize=10, markeredgewidth=1.8, color="red", zorder=8)
        ax.text(gx + 6, gy + 6, str(station_label), fontsize=28, fontweight="bold", color="black",
                bbox=dict(facecolor="white", edgecolor="none", pad=4.0), zorder=9)


    ax.set_xlim(lon_min, lon_max)
    ax.set_ylim(lat_min, lat_max)
    ax.set_xlabel("Longitude (deg)")
    ax.set_ylabel("Latitude (deg)")
    ax.set_title("Shortest water-path distance (Pacific-centered)")


    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ------------------------ CLI ------------------------
def main():
    ap = argparse.ArgumentParser(description="Plot water-path distance PNGs for station NetCDFs.")
    ap.add_argument("--mesh", default="fort.14", help="ADCIRC mesh (fort.14)")
    ap.add_argument("--nc-dir", default=".", help="Directory containing station*.nc files")
    ap.add_argument("--pattern", default="station*.nc", help="Glob pattern for station files")
    ap.add_argument("--outdir", default="figs", help="Directory to save PNGs")
    ap.add_argument("--pac-center", type=float, default=180.0, help="Pacific wrap center (e.g., 180 or 150)")
    ap.add_argument("--lon-min", type=float, default=130.0)
    ap.add_argument("--lon-max", type=float, default=230.0)
    ap.add_argument("--lat-min", type=float, default=-25.0)
    ap.add_argument("--lat-max", type=float, default=55.0)
    ap.add_argument("--max-val", type=float, default=7500.0)
    ap.add_argument("--levels-step", type=float, default=100.0, help="Contour interval (km)")
    ap.add_argument("--station-id", type=str,
                    help="If provided, only plot this station ID (e.g., 17).")
    ap.add_argument("--cmap", type=str,
                    help="Path to a MATLAB .mat file containing 'cm' colormap array")
    ap.add_argument("--make-html", action="store_true",
                    help="If set, create an HTML gallery (index.html) in the output directory")
    args = ap.parse_args()


    os.makedirs(args.outdir, exist_ok=True)


    # Optional: load MATLAB colormap (variable name must be 'cm')
    if args.cmap:
        try:
            import scipy.io as sio
            from matplotlib.colors import ListedColormap
            mat_data = sio.loadmat(args.cmap)
            cm_array = mat_data["cm"]
            if cm_array.max() > 1:
                cm_array = cm_array / 255.0
            globals()["custom_cmap"] = ListedColormap(cm_array, name="my_colormap")
            print(f"Loaded custom colormap from {args.cmap}")
        except Exception as e:
            print(f"Warning: failed to load colormap from {args.cmap}: {e}")
            print("Proceeding with default colormap (viridis).")


    # 1) Read mesh once
    node_ids, x_raw, y, z, tris, sorted_ids, order, title = read_f14_elems_first(args.mesh)
    x = pacific_wrap(x_raw, args.pac_center)


    # Reusable triangulation object (full domain)
    tri_full = Triangulation(x, y, triangles=tris)


    # 2) Gather station files
    files = sorted(glob.glob(os.path.join(args.nc_dir, args.pattern)))
    if not files:
        print("No station NetCDFs found with the given --nc-dir and --pattern.")
        return


    # Optional: restrict to one station by ID
    if args.station_id is not None:
        target = str(int(args.station_id))  # normalize '0017' -> '17'
        files = [f for f in files if get_station_id_from_filename(f) == target]
        if not files:
            print(f"No NetCDF found matching station{target}.nc")
            return


    print(f"Mesh title: {title}")
    print(f"Found {len(files)} station NetCDF(s) in {args.nc_dir}")


    # 3) Plot each file
    pngs = []
    for i, fn in enumerate(files, 1):
        station_label = get_station_id_from_filename(fn)
        with nc.Dataset(fn) as ds:
            dis = np.array(ds.variables["dis"][:], dtype=float)
            node0_id = int(ds.variables["node0"][0])


        out_png = os.path.join(args.outdir, f"{station_label}_Pacific_WD.png")
        print(f"[{i:02d}/{len(files)}] {os.path.basename(fn)} → {os.path.basename(out_png)}")


        plot_one_station(
            tri_full=tri_full,
            x_wrapped=x,
            y=y,
            dis=dis,
            node0_id=node0_id,
            sorted_ids=sorted_ids,
            order=order,
            station_label=station_label,
            out_png=out_png,
            lon_min=args.lon_min,
            lon_max=args.lon_max,
            lat_min=args.lat_min,
            lat_max=args.lat_max,
            max_val=args.max_val,
            levels_step=args.levels_step,
        )
        pngs.append(os.path.basename(out_png))


    # 4) Optional HTML gallery
    if args.make_html:
        html_path = os.path.join(args.outdir, "index.html")
        with open(html_path, "w") as f:
            f.write(
                "<!doctype html><html><head><meta charset='utf-8'>"
                "<meta name='viewport' content='width=device-width,initial-scale=1'>"
                "<title>Station Distance Maps</title>"
                "<style>"
                "body{font-family:system-ui,Segoe UI,Roboto,Arial,sans-serif;margin:20px;}"
                ".grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(320px,1fr));gap:16px;}"
                ".card{border:1px solid #ddd;border-radius:8px;overflow:hidden;box-shadow:0 1px 3px rgba(0,0,0,.08);}"
                ".card h3{margin:0;padding:10px 12px;font-size:16px;background:#f7f7f7;border-bottom:1px solid #eee;}"
                ".card img{width:100%;height:auto;display:block;}"
                "</style></head><body>"
                "<h1>Station Distance Maps</h1><div class='grid'>"
            )
            for png in sorted(pngs, key=lambda s: (len(s), s)):
                f.write(f"<div class='card'><h3>{png}</h3><img src='{png}' alt='{png}'></div>")
            f.write("</div></body></html>")
        print(f"HTML gallery written to {html_path}")


if __name__ == "__main__":
    main()

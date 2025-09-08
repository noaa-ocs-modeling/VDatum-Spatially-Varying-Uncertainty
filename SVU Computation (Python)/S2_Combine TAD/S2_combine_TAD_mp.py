#!/usr/bin/env python3
"""
Combine all  xx_hh_ll_mp*.nc into 3 NetCDF files with fixed shapes:

  TADxx.nc  : datums(datum=7, station=N_total)
  TADhh.nc  : Htime_Hval(high_low=480, station=N_total)   # [0:239]=times, [240:479]=values
  TADll.nc  : Ltime_Lval(high_low=480, station=N_total)   # [0:239]=times, [240:479]=values

Notes
- We PRESERVE the original encoding:
  * Highs: 'higher-highs' are stored as (value + 100)
  * Lows : 'lower-lows'  are stored as (value - 100)
  * Unused slots are filled with -9999.9
  * station are mesh nodes
- This script DOES NOT decode markers or strip padding; correlation coefficient code(S3) should decode before analysis.

How to run
python S2_combine_TAD_mp.py --pattern *xx_hh_ll_mp*.nc" --out-prefix ""   --chunk-stations 8192
Chunk stations is related to Pacific grid.

Created in Matlab by Liujua.Tang@noaa.gov 01/27/2020
Converted to python and changed the output to the netcdf format instead of .mat by  mojgan.rostaminia@noaa.gov, 09/08/2025
"""


import argparse, glob, sys
import numpy as np
from netCDF4 import Dataset

FILLV = -9999.9
CAP   = 240           # split index (0..239 times, 240..479 values)

def scan_total_stations(file_list):
    total = 0
    widths = set()
    for fn in file_list:
        with Dataset(fn, "r") as ds:
            if "stationN" not in ds.variables:
                raise ValueError(f"{fn}: missing 'stationN'")
            nsta = ds.variables["stationN"].shape[0]
            total += nsta
            wH = ds.variables["Htime_Hval"].shape[1]
            wL = ds.variables["Ltime_Lval"].shape[1]
            widths.add((wH, wL))
    for (wH, wL) in widths:
        if (wH != 2*CAP) or (wL != 2*CAP):
            raise ValueError(f"Expected width {2*CAP}; got H={wH}, L={wL}")
    return total

def create_out_TADxx(path, N_total, chunk_stations):
    nc = Dataset(path, "w", format="NETCDF4")
    nc.createDimension("datum", 7)
    nc.createDimension("station", None)  # unlimited
    v_sta = nc.createVariable("stationN", "i4", ("station",))
    v_dat = nc.createVariable(
        "datums", "f4", ("datum", "station"),
        zlib=True, complevel=4,
        chunksizes=(7, min(chunk_stations, max(1, N_total)))
    )
    nc.title = "Combined tidal datums (7 x N_stations)"
    v_dat.long_name = "Datums [MHHW, MHW, MSL, (MHW+MLW)/2, (MHHW+MLLW)/2, MLW, MLLW]"
    return nc, v_sta, v_dat

def create_out_hl(path, var_name, title, N_total, chunk_stations):
    nc = Dataset(path, "w", format="NETCDF4")
    nc.createDimension("high_low", 2*CAP)   # 480
    nc.createDimension("station", None)
    v_sta = nc.createVariable("stationN", "i4", ("station",))
    v_hl  = nc.createVariable(
        var_name, "f4", ("high_low", "station"),
        zlib=True, complevel=4, fill_value=FILLV,
        chunksizes=(2*CAP, min(chunk_stations, max(1, N_total)))
    )
    nc.title = title
    v_hl.long_name = (
        "First half (0..239)=times (s since model epoch), "
        "second half (240..479)=values; highs have +100 on higher-highs, "
        "lows have -100 on lower-lows. Unused slots = -9999.9"
    )
    return nc, v_sta, v_hl

def normalize_datums(dat):
    """
    Ensure datums are (station,7). Accepts (station,7) or (7,station).
    Returns (station,7) float32.
    """
    dat = np.array(dat)
    if dat.ndim != 2 or 7 not in dat.shape:
        raise ValueError("datums must be 2-D and have a size-7 axis")
    if dat.shape[1] == 7:   # (station,7)
        out = dat.astype(np.float32, copy=False)
    else:                   # (7,station) -> transpose
        out = dat.T.astype(np.float32, copy=False)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pattern", default="xx_hh_ll_mp*.nc", help="glob for input part files")
    ap.add_argument("--out-prefix", default="", help="prefix for output files (e.g., R58_)")
    ap.add_argument("--chunk-stations", type=int, default=8192, help="station chunk size for output")
    args = ap.parse_args()

    files = sorted(glob.glob(args.pattern))
    if not files:
        print(f"No files match: {args.pattern}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(files)} files:")
    for fn in files:
        print("  ", fn)

    # Pass 1: total station count + sanity
    N_total = scan_total_stations(files)
    print(f"Total stations across all parts: {N_total}")

    # Create outputs
    path_xx = f"{args.out_prefix}TADxx.nc"
    path_hh = f"{args.out_prefix}TADhh.nc"
    path_ll = f"{args.out_prefix}TADll.nc"

    nc_xx, sta_xx, dat_xx = create_out_TADxx(path_xx, N_total, args.chunk_stations)
    nc_hh, sta_hh, hl_hh  = create_out_hl(path_hh, "Htime_Hval",
                                          "Combined highs (times/values packed)", N_total, args.chunk_stations)
    nc_ll, sta_ll, hl_ll  = create_out_hl(path_ll, "Ltime_Lval",
                                          "Combined lows  (times/values packed)", N_total, args.chunk_stations)

    try:
        # Append along station dimension as we go
        offset = 0
        for fn in files:
            print(f"Reading {fn} …")
            with Dataset(fn, "r") as ds:
                stationN = np.array(ds.variables["stationN"][:], dtype=np.int32)   # (nsta,)
                datums   = normalize_datums(ds.variables["datums"][:])             # (nsta,7) float32
                H = np.array(ds.variables["Htime_Hval"][:], dtype=np.float32)      # (nsta,480)
                L = np.array(ds.variables["Ltime_Lval"][:], dtype=np.float32)      # (nsta,480)

            nsta = stationN.shape[0]
            if H.shape != (nsta, 2*CAP) or L.shape != (nsta, 2*CAP):
                raise ValueError(f"{fn}: unexpected H/L shapes; got {H.shape} and {L.shape}")

            # write station ids
            sta_slice = slice(offset, offset + nsta)
            sta_xx[sta_slice] = stationN
            sta_hh[sta_slice] = stationN
            sta_ll[sta_slice] = stationN

            # datums: (datum, station) → transpose from (nsta,7)
            dat_xx[:, sta_slice] = datums.T

            # highs/lows: (high_low, station) → transpose from (nsta,480)
            hl_hh[:, sta_slice] = H.T
            hl_ll[:, sta_slice] = L.T

            offset += nsta
            print(f"  wrote stations {sta_slice.start}..{sta_slice.stop-1}")

        print(f"Done. Wrote:\n  {path_xx}\n  {path_hh}\n  {path_ll}")
    finally:
        nc_xx.close()
        nc_hh.close()
        nc_ll.close()

if __name__ == "__main__":
    main()



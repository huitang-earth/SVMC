#!/usr/bin/env python3
"""
tabpfn_alloc.py  -  TabPFN-based carbon allocation ratio predictor

Predicts daily cratio_leaf (leaf allocation ratio) using TabPFNRegressor,
replacing the fixed parameter in the SVMC Fortran model.

Three modes
-----------
train:
    Load multi-site historical data, train TabPFN, save model to a pickle file.
    python tabpfn_alloc.py train --sites BE-Lon CH-Oe2 CZ-KrP ... \\
                                  --svm_output_pattern "/path/SVM_{site}_hr-day*.nc" \\
                                  --climate_pattern  "/path/FieldObs_{site}.*.hr.timeshift_era.nc" \\
                                  --manage_pattern   "/path/FieldObs_{site}.*.management.nc" \\
                                  --index_best 0 \\
                                  --model_out tabpfn_cratio.pkl

predict:
    Load saved model, read input NetCDF files for ONE simulation run, write daily
    cratio_leaf / cratio_root to a NetCDF file.
    python tabpfn_alloc.py predict \\
                                   --model tabpfn_cratio.pkl \\
                                   --svm_output  /path/SVM_{site}_hr-day*.nc \\
                                   --climate_dir /path/FieldObs_{site}.*.hr.timeshift_era.nc \\
                                   --manage_dir  /path/FieldObs_{site}.*.management.nc \\
                                   --index_best 0 \\
                                   --output cratio_tabpfn.nc

predict_from_csv:
    Load saved model, read the daily feature CSV written by SVMC Fortran after each run,
    and write cratio_leaf / cratio_root to cratio_tabpfn.csv for the next run.
    python tabpfn_alloc.py predict_from_csv \\
                                   --model tabpfn_cratio.pkl \\
                                   --input tabpfn_features.csv \\
                                   --output cratio_tabpfn.csv \\
                                   --latitude 51.0

The features used (matching the notebook df_rowstack columns):
  water_vapor_saturation_deficit, precipitation_flux,
  surface_downwelling_shortwave_flux_in_air, air_temperature,
  air_temperature_7 (7-day rolling mean), air_temperature_14 (14-day rolling mean),
  Dpsi, Chi, SoilMoistPot, GPP,
  harvest_seq, fertilizer_seq, grazing_seq, organic_material_seq, mowing_seq,
  day_length
"""

import argparse
import glob
import os
import pickle
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from huggingface_hub import login, get_token
from scipy.signal import savgol_filter

# Migrate deprecated env var to new name before PyTorch is imported
if "PYTORCH_CUDA_ALLOC_CONF" in os.environ and "PYTORCH_ALLOC_CONF" not in os.environ:
    os.environ["PYTORCH_ALLOC_CONF"] = os.environ.pop("PYTORCH_CUDA_ALLOC_CONF")

from tabpfn import TabPFNRegressor
from tabpfn_time_series import TimeSeriesDataFrame, FeatureTransformer
from tabpfn_time_series.data_preparation import generate_test_X
from tabpfn_time_series.features import RunningIndexFeature, CalendarFeature, AutoSeasonalFeature
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from sklearn.model_selection import train_test_split

def _hf_login():
    """
    Ensure a HuggingFace token is available for TabPFN weight access.

    First run:
      Set HF_TOKEN in your environment (export HF_TOKEN=hf_...) before running.
      login() is called once, which writes the token to the local cache file
      (~/.cache/huggingface/token).

    All subsequent runs (including Fortran CALL SYSTEM invocations):
      The cached token is found by get_token() and huggingface_hub uses it
      automatically — login() is skipped entirely, so there is no network
      round-trip overhead.
    """
    if get_token():
        # Token already cached from a prior login — nothing to do.
        print("HuggingFace Hub: using cached token.")
        return

    token = os.environ.get("HF_TOKEN")
    if not token:
        raise RuntimeError(
            "No HuggingFace token found. On first use:\n"
            "  export HF_TOKEN=hf_...   # your token from huggingface.co/settings/tokens\n"
            "  python tabpfn_alloc.py train ...  (or predict / predict_from_csv)\n"
            "The token will be saved to cache; subsequent runs need no HF_TOKEN."
        )
    login(token=token, add_to_git_credential=False)
    print("HuggingFace Hub: token saved to cache. Future runs will skip login.")


FEATURE_COLS = [
    "water_vapor_saturation_deficit",
    "precipitation_flux",
    "surface_downwelling_shortwave_flux_in_air",
    "air_temperature",
    "air_temperature_7",
    "air_temperature_14",
    "Dpsi",
    "Chi",
    "SoilMoistPot",
    "GPP",
    "harvest_seq",
    "fertilizer_seq",
    "grazing_seq",
    "organic_material_seq",
    "mowing_seq",
    "day_length",
]


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def _day_length(dates, latitude_deg=51.0):
    """Approximate astronomical day length (hours) from date and latitude."""
    lat = np.radians(latitude_deg)
    doy = dates.day_of_year.values
    # Solar declination (Spencer 1971)
    B = 2 * np.pi * (doy - 1) / 365.0
    decl = (0.006918 - 0.399912 * np.cos(B) + 0.070257 * np.sin(B)
            - 0.006758 * np.cos(2 * B) + 0.000907 * np.sin(2 * B))
    cos_hour = -np.tan(lat) * np.tan(decl)
    cos_hour = np.clip(cos_hour, -1.0, 1.0)
    return 2 * np.degrees(np.arccos(cos_hour)) / 15.0


def _management_sequences(manage_da):
    """
    Convert management_type DataArray (daily) into five sequence columns,
    matching the original notebook logic:

    - Days 1–13 after an event are marked with values 2–14 (m+1 where m=1..13).
    - The event day itself stays 0.
    - The sequence stops early if another event of the same type occurs.
    - All other days are 0.

    Management type codes: 1=harvest, 2=fertilizer, 3=grazing, 4=organic, 5=mowing
    """
    types = {1: "harvest_seq", 2: "fertilizer_seq",
             3: "grazing_seq", 4: "organic_material_seq", 5: "mowing_seq"}
    mt = manage_da.values.astype(float)
    n = len(mt)
    result = {}
    for code, col in types.items():
        event = (mt == code)
        seq = np.zeros(n, dtype=float)
        for i in range(n):
            if event[i]:
                for m in range(1, 14):
                    if i + m >= n:
                        break
                    if event[i + m]:
                        break  # another event starts — stop
                    seq[i + m] = m + 1
        result[col] = seq
    return result


def _read_latitude_from_nc(clim_files):
    """Read site latitude from the first climate NetCDF file."""
    ds = xr.open_dataset(clim_files[0], decode_times=False)
    for name in ("lat", "latitude", "Lat", "Latitude"):
        if name in ds.coords or name in ds:
            lat = float(ds[name].values.flat[0])
            ds.close()
            return lat
    ds.close()
    raise KeyError(
        f"Could not find latitude coordinate in {clim_files[0]}. "
        "Expected one of: lat, latitude, Lat, Latitude."
    )


def load_features_for_site(climate_pattern, manage_pattern, svm_output_pattern,
                             index_best=0):
    """
    Load all 16 features from NetCDF files for one site.
    Latitude is read directly from the climate NetCDF file.

    Returns
    -------
    df_X : pd.DataFrame with FEATURE_COLS columns, daily index
    """

    # --- Climate (hourly → daily) ---
    clim_files = sorted(glob.glob(climate_pattern))
    if not clim_files:
        raise FileNotFoundError(f"No climate files found: {climate_pattern}")
    dclim = xr.open_mfdataset(clim_files, combine="nested", concat_dim="time",
                               decode_times=True)
    vpd_da   = dclim["water_vapor_saturation_deficit"][:, 0, 0].resample(time="1D").mean()
    prec_da  = dclim["precipitation_flux"][:, 0, 0].resample(time="1D").mean()
    rad_da   = dclim["surface_downwelling_shortwave_flux_in_air"][:, 0, 0].resample(time="1D").mean()
    temp_da  = dclim["air_temperature"][:, 0, 0].resample(time="1D").mean()
    clim_idx = pd.DatetimeIndex(temp_da.time.values)
    vpd    = pd.Series(vpd_da.values,  index=clim_idx, name="vpd")
    prec   = pd.Series(prec_da.values, index=clim_idx, name="prec")
    rad    = pd.Series(rad_da.values,  index=clim_idx, name="rad")
    temp_s = pd.Series(temp_da.values, index=clim_idx, name="air_temperature")
    dclim.close()

    # Read latitude from the NetCDF file (same source as Fortran's netCDF_readlonlat)
    latitude_deg = _read_latitude_from_nc(clim_files)
    print(f"    Latitude read from NetCDF: {latitude_deg:.4f} deg")

    temp_r7  = temp_s.rolling(7,  min_periods=1).mean().rename("air_temperature_7")
    temp_r14 = temp_s.rolling(14, min_periods=1).mean().rename("air_temperature_14")

    # --- SVM model output (GPP, Chi, Dpsi, SoilMoistPot) ---
    svm_files = sorted(glob.glob(svm_output_pattern))
    if not svm_files:
        raise FileNotFoundError(f"No SVM output files found: {svm_output_pattern}")
    svm = xr.open_dataset(svm_files[0], decode_times=True)
    # Select one or more ensemble members and average if multiple indices given
    def _extract(var):
        da = svm[var]
        if da.ndim == 4:          # (expid, time, y, x)
            da = da[index_best, :, 0, 0].mean(dim=da.dims[0])
        elif da.ndim == 3:        # (time, y, x)
            da = da[:, 0, 0]
        return pd.Series(da.values, index=pd.DatetimeIndex(svm.time.values), name=var)
    dpsi      = _extract("Dpsi")
    chi       = _extract("Chi")
    gpp       = _extract("GPP")
    soilmoist = _extract("SoilMoistPot")
    svm.close()

    # Use dpsi.index as the master index; reindex everything else to it
    common_idx = clim_idx

    # --- Management ---
    ma_files = sorted(glob.glob(manage_pattern))
    if not ma_files:
        raise FileNotFoundError(f"No management files found: {manage_pattern}")
    dma = xr.open_mfdataset(ma_files, combine="nested", concat_dim="time", decode_times=True)
    ma_type = dma["management_type"][:, 0, 0]
    ma_idx  = pd.DatetimeIndex(ma_type.time.values)
    dma.close()
    seq_dict = _management_sequences(ma_type)

    # --- Day length ---
    dl = _day_length(pd.DatetimeIndex(common_idx), latitude_deg=latitude_deg)

    # --- Assemble DataFrame ---
    df = pd.DataFrame(index=common_idx)
    df["water_vapor_saturation_deficit"]              = vpd.values
    df["precipitation_flux"]                          = prec.values
    df["surface_downwelling_shortwave_flux_in_air"]   = rad.values
    df["air_temperature"]                             = temp_s.values
    df["air_temperature_7"]                           = temp_r7.values
    df["air_temperature_14"]                          = temp_r14.values
    df["Dpsi"]                                        = dpsi.values
    df["Chi"]                                         = chi.values
    df["SoilMoistPot"]                                = soilmoist.values
    df["GPP"]                                         = gpp.values

    for col, arr in seq_dict.items():
        df[col] = pd.Series(arr, index=common_idx).values

    df["day_length"] = dl
    
    #df = df.dropna()


    return df[FEATURE_COLS]


def load_target_for_site(svm_output_pattern, index_best=0, savgol_window=31, savgol_order=2):
    """
    Load and smooth the Jmax / cratio_leaf (Allocation_ratio) target variable
    from SVM daily output NetCDF.

    In the notebook, the variable stored as 'Jmax' in the daily SVM output
    is the inverted cratio_leaf (= Allocation_ratio).
    """
    svm_files = sorted(glob.glob(svm_output_pattern))
    if not svm_files:
        raise FileNotFoundError(f"No SVM output files found: {svm_output_pattern}")
    svm = xr.open_dataset(svm_files[0], decode_times=True)
    da = svm["Jmax"]
    if da.ndim == 4:
        da = da[index_best, :, 0, 0].mean(dim=da.dims[0])
    elif da.ndim == 3:
        da = da[:, 0, 0]
    idx = pd.DatetimeIndex(svm.time.values)
    svm.close()

    series = pd.Series(da.values, index=idx, name="Allocation_ratio")
    # Savitzky-Golay smoothing (same as notebook)
    smoothed = savgol_filter(series.values, window_length=savgol_window, polyorder=savgol_order)
    smoothed = np.clip(smoothed, 0.1, 0.9)
    return pd.Series(smoothed, index=idx, name="Allocation_ratio")


# ---------------------------------------------------------------------------
# Time-series feature helpers
# ---------------------------------------------------------------------------

def _ts_features_train(df_X, df_y, item_id):
    """
    Build a TimeSeriesDataFrame from one site's features + target,
    apply RunningIndexFeature / CalendarFeature / AutoSeasonalFeature,
    and return the transformed training slice (last step held out for
    prediction_length=1, matching the notebook pattern).
    """
    df = df_X.copy()
    df.index = df.index.rename("timestamp")
    df["target"] = df_y.values
    df["item_id"] = item_id
    df = df.set_index("item_id", append=True).reorder_levels(["item_id", "timestamp"])

    tsdf = TimeSeriesDataFrame(df)
    train_ts, _ = tsdf.train_test_split(prediction_length=1)
    test_ts = generate_test_X(train_ts, 1)

    features = [RunningIndexFeature(), CalendarFeature()] #, AutoSeasonalFeature()]
    tsdf_feat, _ = FeatureTransformer(features).transform(tsdf, test_ts)
    return tsdf_feat


def _ts_features_predict(df_X):
    """
    Apply the same time-series features to a prediction-only DataFrame
    (no target available).  The full time series is used as context so
    that AutoSeasonalFeature can detect seasonality; all n rows are
    returned ready for model.predict().
    """
    df = df_X.copy()
    df.index = df.index.rename("timestamp")
    df["target"] = 0.0          # dummy — not used for prediction
    df["item_id"] = 0
    df = df.set_index("item_id", append=True).reorder_levels(["item_id", "timestamp"])

    tsdf = TimeSeriesDataFrame(df)
    train_ts, _ = tsdf.train_test_split(prediction_length=1)
    test_ts = generate_test_X(train_ts, 1)

    features = [RunningIndexFeature(), CalendarFeature()] #, AutoSeasonalFeature()]
    tsdf_feat, _ = FeatureTransformer(features).transform(tsdf, test_ts)
    return tsdf_feat.drop(columns=["target"])
    
    #features = [RunningIndexFeature(), CalendarFeature()] #, AutoSeasonalFeature()]
    #train_ts, _ = FeatureTransformer(features).transform(train_ts, test_ts)
    #return train_ts.drop(columns=["target"])
    
    


# ---------------------------------------------------------------------------
# check_input
# ---------------------------------------------------------------------------

def check_input(df_X, df_y=None, site="site", out_dir="."):
    """
    Diagnostic checks and plots for TabPFN input preparation.

    Checks
    ------
    1. Missing values per feature column (printed + flagged if any).
    2. Value-range sanity (min / mean / max per column).
    3. Index continuity — flags gaps larger than 1 day.
    4. df_X / df_y length mismatch.

    Plots saved to out_dir
    ----------------------
    {site}_features.png    — one subplot per feature column
    {site}_target.png      — target time series (if df_y provided)
    {site}_nan_heatmap.png — red = missing value
    """
    import matplotlib.dates as mdates

    os.makedirs(out_dir, exist_ok=True)
    issues = []

    # Flatten multi-level index (item_id, timestamp) produced by ts feature transforms
    if isinstance(df_X.index, pd.MultiIndex):
        df_X = df_X.copy()
        df_X.index = df_X.index.get_level_values("timestamp")

    cols = df_X.columns.tolist()

    # 1. Missing values
    nan_counts = df_X.isna().sum()
    print(f"\n[check_input] {site} — NaN counts per feature:")
    print(nan_counts.to_string())
    if nan_counts.any():
        issues.append(f"  NaNs in: {list(nan_counts[nan_counts > 0].index)}")

    # 2. Value ranges
    print(f"\n[check_input] {site} — feature statistics:")
    print(df_X.describe().T[["min", "mean", "max"]].to_string())

    # 3. Index continuity
    gaps = pd.Series(df_X.index).diff().dt.days
    large_gaps = gaps[gaps > 1]
    if not large_gaps.empty:
        issues.append(f"  Time gaps > 1 day at positions: {list(large_gaps.index)}")

    # 4. Alignment with df_y
    if df_y is not None and len(df_X) != len(df_y):
        issues.append(f"  df_X length ({len(df_X)}) != df_y length ({len(df_y)})")

    if issues:
        print(f"\n[check_input] WARNINGS for {site}:")
        for w in issues:
            print(w)
    else:
        print(f"\n[check_input] {site} — all checks passed.")

    # Plot 1: feature time series (all columns)
    ncols = 2
    nrows = (len(cols) + 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, nrows * 2.5), sharex=True)
    axes = axes.flatten()
    for i, col in enumerate(cols):
        axes[i].plot(df_X.index, df_X[col], lw=0.6)
        axes[i].set_title(col, fontsize=8)
        axes[i].xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)
    fig.suptitle(f"{site} — input features", fontsize=10)
    fig.tight_layout()
    path1 = os.path.join(out_dir, f"{site}_features.png")
    fig.savefig(path1, dpi=120)
    plt.close(fig)
    print(f"  Saved: {path1}")

    # Plot 2: target
    if df_y is not None:
        fig, ax = plt.subplots(figsize=(10, 3))
        ax.plot(df_y.index, df_y.values, lw=0.8, color="tab:orange")
        ax.set_title(f"{site} — target (Allocation_ratio)")
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        fig.tight_layout()
        path2 = os.path.join(out_dir, f"{site}_target.png")
        fig.savefig(path2, dpi=120)
        plt.close(fig)
        print(f"  Saved: {path2}")

    # Plot 3: NaN heatmap (all columns)
    fig, ax = plt.subplots(figsize=(10, max(4, len(cols) * 0.35)))
    nan_matrix = df_X[cols].isna().astype(int).T
    ax.imshow(nan_matrix, aspect="auto", cmap="Reds", interpolation="none")
    ax.set_yticks(range(len(cols)))
    ax.set_yticklabels(cols, fontsize=7)
    ax.set_xlabel("time index")
    ax.set_title(f"{site} — NaN heatmap (red = missing)")
    fig.tight_layout()
    path3 = os.path.join(out_dir, f"{site}_nan_heatmap.png")
    fig.savefig(path3, dpi=120)
    plt.close(fig)
    print(f"  Saved: {path3}")


# ---------------------------------------------------------------------------
# train
# ---------------------------------------------------------------------------

def train(args):
    _hf_login()
    print("=== Training TabPFN allocation model ===")
    all_ts = []   # transformed TimeSeriesDataFrames, one per site

    for site_idx, site in enumerate(args.sites):
        print(f"  Loading site: {site}")
        try:
            df_X = load_features_for_site(
                climate_pattern=args.climate_pattern.replace("{site}", site),
                manage_pattern=args.manage_pattern.replace("{site}", site),
                svm_output_pattern=args.svm_output_x_pattern.replace("{site}", site),
                index_best=args.index_best,
            )
            df_y = load_target_for_site(
                svm_output_pattern=args.svm_output_y_pattern.replace("{site}", site),
                index_best=args.index_best,
            )
            df_X = _filter_date_ranges(df_X, args.train_dates)
            df_y = _filter_date_ranges(df_y, args.train_dates)
            print(f"    After date filtering: {len(df_X)} rows (X), {len(df_y)} rows (y).")
            if args.debug:
                check_input(df_X, df_y, site=site, out_dir=args.debug_dir)
            if args.ts_features:
                train_ts = _ts_features_train(df_X,
                                              df_y,
                                              item_id=site_idx)
                if args.debug:
                    check_input(pd.DataFrame(train_ts.drop(columns=["target"])).reset_index(level="item_id", drop=True),
                                train_ts["target"].reset_index(level="item_id", drop=True),
                                site=f"{site}_ts_features", out_dir=args.debug_dir)
            else:
                train_ts = df_X.copy()
                train_ts["target"] = df_y.values
            all_ts.append(train_ts)
            print(f"    {len(train_ts)} rows loaded.")
        except Exception as e:
            print(f"    WARNING: {e} — skipping {site}")

    if not all_ts:
        sys.exit("No training data could be loaded. Aborting.")

    crop_train = pd.concat(all_ts, axis=0, ignore_index=False)
    X = crop_train.drop(columns=["target"])
    y = crop_train["target"]

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )

    print(f"Training on {len(X_train)} samples, testing on {len(X_test)} samples.")

    if args.debug:
        csv_path = os.path.join(args.debug_dir, "tabpfn_train_input.csv")
        crop_train.to_csv(csv_path)
        print(f"Debug: training input written to {csv_path}")

    model = TabPFNRegressor(device="auto", ignore_pretraining_limits=True)
    model.fit(X_train, y_train)

    preds = model.predict(X_test)
    print(f"Test  MSE : {mean_squared_error(y_test, preds):.4f}")
    print(f"Test  MAE : {mean_absolute_error(y_test, preds):.4f}")
    print(f"Test  R²  : {r2_score(y_test, preds):.4f}")

    with open(args.model_out, "wb") as f:
        pickle.dump(model, f)
    print(f"Model saved to {args.model_out}")


# ---------------------------------------------------------------------------
# predict
# ---------------------------------------------------------------------------

def predict(args):
    _hf_login()
    print("=== Predicting cratio_leaf with TabPFN ===")

    with open(args.model, "rb") as f:
        model = pickle.load(f)
    print(f"Model loaded from {args.model}")

    site = args.sites[0]
    df_X = load_features_for_site(
        climate_pattern=args.climate_pattern.replace("{site}", site),
        manage_pattern=args.manage_pattern.replace("{site}", site),
        svm_output_pattern=args.svm_output_x_pattern.replace("{site}", site),
        index_best=args.index_best,
    )
    df_X = _filter_date_ranges(df_X, args.predict_dates)
    print(f"Features loaded: {len(df_X)} daily time steps (after date filtering).")

    X = _ts_features_predict(df_X) if args.ts_features else df_X

    if args.debug:
        X_debug = pd.DataFrame(X)
        if isinstance(X_debug.index, pd.MultiIndex):
            X_debug.index = X_debug.index.get_level_values("timestamp")
        check_input(X_debug, out_dir=args.debug_dir)
        csv_path = os.path.join(args.debug_dir, "tabpfn_predict_input.csv")
        X.to_csv(csv_path)
        print(f"Debug: predict input written to {csv_path}")

    preds = model.predict(X, output_type="quantiles", quantiles=[0.1, 0.5, 0.9])
    cratio_leaf_q10 = np.clip(preds[0], 0.1, 0.9)
    cratio_leaf     = np.clip(preds[1], 0.1, 0.9)  # median = main prediction
    cratio_leaf_q90 = np.clip(preds[2], 0.1, 0.9)
    cratio_root = 1.0 - cratio_leaf

    # --- Write output CSV ---
    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    times = X.index.get_level_values("timestamp") if isinstance(X.index, pd.MultiIndex) else X.index
    out_df = pd.DataFrame({
        "cratio_leaf_q10": cratio_leaf_q10,
        "cratio_leaf":     cratio_leaf,
        "cratio_leaf_q90": cratio_leaf_q90,
        "cratio_root":     cratio_root,
    }, index=times)

    # --- Optional validation against ground truth ---
    if args.svm_output_y_pattern is not None:
        df_y = load_target_for_site(
            svm_output_pattern=args.svm_output_y_pattern.replace("{site}", site),
            index_best=args.index_best,
        )
        df_y = _filter_date_ranges(df_y, args.predict_dates)
        n = min(len(times), len(df_y))
        y_true = df_y.values[:n]
        y_pred = cratio_leaf[:n]
        y_q10  = cratio_leaf_q10[:n]
        y_q90  = cratio_leaf_q90[:n]
        common = times[:n]
        true_col = np.full(len(times), np.nan)
        true_col[:n] = y_true
        out_df["cratio_leaf_true"] = true_col

        print("\n  ── Validation metrics ──────────────────────────────")
        print(f"  N points : {len(y_true)}")
        print(f"  MSE      : {mean_squared_error(y_true, y_pred):.4f}")
        print(f"  MAE      : {mean_absolute_error(y_true, y_pred):.4f}")
        print(f"  R²       : {r2_score(y_true, y_pred):.4f}")
        print("  ────────────────────────────────────────────────────\n")

        fig, axes = plt.subplots(2, 1, figsize=(12, 7))

        axes[0].fill_between(common, y_q10, y_q90, alpha=0.2, color="tab:blue", label="10–90% range")
        axes[0].plot(common, y_true, label="Ground truth", color="tab:orange", lw=0.8)
        axes[0].plot(common, y_pred, label="TabPFN median", color="tab:blue", lw=0.8, alpha=0.8)
        axes[0].set_ylabel("Allocation ratio")
        axes[0].set_title("cratio_leaf — prediction vs ground truth")
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)

        axes[1].scatter(y_true, y_pred, s=4, alpha=0.4, color="tab:blue")
        lims = [min(y_true.min(), y_pred.min()), max(y_true.max(), y_pred.max())]
        axes[1].plot(lims, lims, "k--", lw=0.8, label="1:1 line")
        axes[1].set_xlabel("Ground truth")
        axes[1].set_ylabel("Prediction (median)")
        axes[1].set_title(f"R²={r2_score(y_true, y_pred):.3f}  MAE={mean_absolute_error(y_true, y_pred):.3f}")
        axes[1].legend()
        axes[1].grid(True, alpha=0.3)

        fig.tight_layout()
        val_plot = os.path.join(os.path.dirname(args.output) or ".", "tabpfn_validation.png")
        fig.savefig(val_plot, dpi=120)
        plt.close(fig)
        print(f"  Validation plot saved to {val_plot}")

    out_df.to_csv(args.output)
    print(f"Predictions written to {args.output}")
    print(f"  cratio_leaf range: [{cratio_leaf.min():.3f}, {cratio_leaf.max():.3f}]")
    print(f"  cratio_root range: [{cratio_root.min():.3f}, {cratio_root.max():.3f}]")


# ---------------------------------------------------------------------------
# predict_from_csv  (called by SVMC Fortran via CALL SYSTEM)
# ---------------------------------------------------------------------------

def predict_from_csv(args):
    """
    Read the daily feature CSV written by SVMC.f90 after the main time loop.
    The CSV columns written by Fortran are:
      air_temperature, water_vapor_saturation_deficit,
      surface_downwelling_shortwave_flux_in_air, precipitation_flux,
      Dpsi, Chi, SoilMoistPot, GPP, management_type, doy

    This function derives the remaining features needed by the model
    (rolling temperature means, management sequences, day_length) and
    writes cratio_tabpfn.csv with columns: cratio_leaf, cratio_root.
    """
    print("=== predict_from_csv: building features from Fortran CSV ===")

    _hf_login()
    df = pd.read_csv(args.input)
    print(f"  Loaded {len(df)} daily rows from {args.input}")

    # Reconstruct proper DatetimeIndex from year + doy written by Fortran
    year = int(df["year"].iloc[0])
    df.index = pd.to_datetime(year * 1000 + df["doy"], format="%Y%j")
    df.index.name = "timestamp"

    # Rolling temperature means are computed in Fortran and written directly to the CSV.
    # Management sequences are computed in Fortran and written directly to the CSV.
    # Expected columns: harvest_seq, fertilizer_seq, grazing_seq, organic_material_seq, mowing_seq

    # Day length from doy and latitude
    doy = df["doy"].values
    lat = np.radians(args.latitude)
    B = 2 * np.pi * (doy - 1) / 365.0
    decl = (0.006918 - 0.399912 * np.cos(B) + 0.070257 * np.sin(B)
            - 0.006758 * np.cos(2 * B) + 0.000907 * np.sin(2 * B))
    cos_hour = np.clip(-np.tan(lat) * np.tan(decl), -1.0, 1.0)
    df["day_length"] = 2 * np.degrees(np.arccos(cos_hour)) / 15.0

    # Apply time-series features (same types as used during training)
    X = _ts_features_predict(df[FEATURE_COLS]) if args.ts_features else df[FEATURE_COLS]

    # Load model and predict
    with open(args.model, "rb") as f:
        model = pickle.load(f)
    print(f"  Model loaded from {args.model}")

    if args.debug:
        X_debug = pd.DataFrame(X)
        if isinstance(X_debug.index, pd.MultiIndex):
            X_debug.index = X_debug.index.get_level_values("timestamp")
        check_input(X_debug, out_dir=args.debug_dir)
        csv_path = os.path.join(args.debug_dir, "tabpfn_predict_from_csv_input.csv")
        X.to_csv(csv_path)
        print(f"  Debug: predict_from_csv input written to {csv_path}")

    preds = model.predict(X, output_type="quantiles", quantiles=[0.1, 0.5, 0.9])
    cratio_leaf_q10 = np.clip(preds[0], 0.1, 0.9)
    cratio_leaf     = np.clip(preds[1], 0.1, 0.9)  # median = main prediction
    cratio_leaf_q90 = np.clip(preds[2], 0.1, 0.9)
    cratio_root = 1.0 - cratio_leaf

    print("\n  ── TabPFN prediction summary ──────────────────────────")
    print(f"  Days predicted : {len(cratio_leaf)}")
    print(f"  cratio_leaf    : mean={cratio_leaf.mean():.3f}  "
          f"std={cratio_leaf.std():.3f}  "
          f"range=[{cratio_leaf.min():.3f}, {cratio_leaf.max():.3f}]")
    print(f"  cratio_root    : mean={cratio_root.mean():.3f}  "
          f"std={cratio_root.std():.3f}  "
          f"range=[{cratio_root.min():.3f}, {cratio_root.max():.3f}]")
    print("  ────────────────────────────────────────────────────────\n")

    # Write output CSV (read by SVMC.f90 immediately after this call)
    out = pd.DataFrame({
        "cratio_leaf_q10": cratio_leaf_q10,
        "cratio_leaf":     cratio_leaf,
        "cratio_leaf_q90": cratio_leaf_q90,
        "cratio_root":     cratio_root,
    })
    out.to_csv(args.output, index=False)
    print(f"  Written to {args.output}  ({len(out)} rows)")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_indices(value):
    """Accept a single int or comma-separated ints: '0' or '12,13,14'."""
    return [int(v) for v in value.split(",")]


def _parse_date_ranges(value):
    """
    Parse a single date-range string 'START:END' into a (Timestamp, Timestamp) tuple.
    Each CLI token is one range; argparse nargs="+" collects multiple tokens into a list.
    Example: '--train_dates 2017-01-01:2019-12-31 2020-06-01:2020-12-31'
    """
    parts = value.split(":")
    if len(parts) != 2:
        raise argparse.ArgumentTypeError(
            f"Date range must be 'YYYY-MM-DD:YYYY-MM-DD', got: {value!r}"
        )
    return (pd.Timestamp(parts[0]), pd.Timestamp(parts[1]))


def _filter_date_ranges(obj, ranges):
    """
    Keep only rows whose DatetimeIndex falls within any of the given (start, end) ranges.
    Ranges are inclusive on both ends.
    If ranges is None or empty, return obj unchanged.
    """
    if not ranges:
        return obj
    mask = pd.Series(False, index=obj.index)
    idx_dates = obj.index.normalize()   # strip time-of-day; compare dates only
    for start, end in ranges:
        mask |= (idx_dates >= start.normalize()) & (idx_dates <= end.normalize())
    return obj[mask]


def main():
    parser = argparse.ArgumentParser(description="TabPFN allocation ratio predictor")
    sub = parser.add_subparsers(dest="mode", required=True)

    # --- train ---
    p_train = sub.add_parser("train", help="Train model on multi-site data")
    p_train.add_argument("--sites", nargs="+", required=True,
                         help="Site codes, e.g. BE-Lon CH-Oe2 CZ-KrP")
    p_train.add_argument("--svm_output_x_pattern", required=True,
                         help="Glob pattern for SVM output NC used for input features (X); use {site} placeholder")
    p_train.add_argument("--svm_output_y_pattern", required=True,
                         help="Glob pattern for SVM output NC used for target variable (Y); use {site} placeholder")
    p_train.add_argument("--climate_pattern", required=True,
                         help="Glob pattern for climate input NC files; use {site} placeholder")
    p_train.add_argument("--manage_pattern", required=True,
                         help="Glob pattern for management NC files; use {site} placeholder")
    p_train.add_argument("--index_best", type=_parse_indices, default=[0],
                         help="Ensemble index/indices from SVM output, e.g. 0 or 12,13,14 (default: 0)")
    p_train.add_argument("--train_dates", nargs="+", type=_parse_date_ranges, default=None,
                         help="One or more date ranges to use for training, each as START:END "
                              "(e.g. 2017-01-01:2019-12-31 2020-06-01:2020-12-31). "
                              "Ranges are inclusive. If omitted, all available dates are used.")
    p_train.add_argument("--model_out", default="tabpfn_cratio.pkl",
                         help="Output pickle path for trained model")
    p_train.add_argument("--ts_features", action="store_true",
                         help="Apply RunningIndex/Calendar/AutoSeasonal time-series features before training")
    p_train.add_argument("--debug", action="store_true",
                         help="Run input checks and save diagnostic plots/CSV")
    p_train.add_argument("--debug_dir", default=".",
                         help="Directory for debug plots and CSV (default: current dir)")

    # --- predict ---
    p_pred = sub.add_parser("predict", help="Predict cratio_leaf for a simulation run")
    p_pred.add_argument("--model", required=True,
                        help="Path to trained model pickle")
    p_pred.add_argument("--svm_output_x_pattern", required=True,
                        help="Glob pattern for SVM output NC used for input features (X); use {site} placeholder")
    p_pred.add_argument("--svm_output_y_pattern", required=False, default=None,
                        help="Glob pattern for SVM output NC used for target variable (Y); use {site} placeholder")
    p_pred.add_argument("--climate_pattern", required=True,
                        help="Glob pattern for climate input NC files; use {site} placeholder")
    p_pred.add_argument("--manage_pattern", required=True,
                        help="Glob pattern for management NC files; use {site} placeholder")
    p_pred.add_argument("--sites", nargs="+", required=True,
                        help="Site codes, e.g. Qvidja")
    p_pred.add_argument("--index_best", type=_parse_indices, default=[0],
                        help="Ensemble index/indices from SVM output, e.g. 0 or 12,13,14 (default: 0)")
    p_pred.add_argument("--predict_dates", nargs="+", type=_parse_date_ranges, default=None,
                        help="One or more date ranges to use for prediction/validation, each as START:END "
                             "(e.g. 2020-01-01:2021-12-31 2022-06-01:2022-09-30). "
                             "Ranges are inclusive. If omitted, all available dates are used.")
    p_pred.add_argument("--output", default="cratio_tabpfn.csv",
                        help="Output CSV with daily cratio_leaf/cratio_root (and optionally cratio_leaf_true)")
    p_pred.add_argument("--ts_features", action="store_true",
                        help="Apply RunningIndex/Calendar/AutoSeasonal time-series features before predicting")
    p_pred.add_argument("--debug", action="store_true",
                        help="Run input checks and save diagnostic plots/CSV")
    p_pred.add_argument("--debug_dir", default=".",
                        help="Directory for debug plots and CSV (default: current dir)")

    # --- predict_from_csv ---
    p_csv = sub.add_parser("predict_from_csv",
                           help="Predict from daily feature CSV written by SVMC Fortran")
    p_csv.add_argument("--model", required=True,
                       help="Path to trained model pickle")
    p_csv.add_argument("--input", default="tabpfn_features.csv",
                       help="Feature CSV written by SVMC.f90 (default: tabpfn_features.csv)")
    p_csv.add_argument("--output", default="cratio_tabpfn.csv",
                       help="Output CSV with cratio_leaf, cratio_root (default: cratio_tabpfn.csv)")
    p_csv.add_argument("--latitude", type=float, default=51.0,
                       help="Site latitude in degrees for day-length calculation (default: 51.0)")
    p_csv.add_argument("--ts_features", action="store_true",
                       help="Apply RunningIndex/Calendar/AutoSeasonal time-series features before predicting")
    p_csv.add_argument("--debug", action="store_true",
                       help="Run input checks and save diagnostic plots/CSV")
    p_csv.add_argument("--debug_dir", default=".",
                       help="Directory for debug plots and CSV (default: current dir)")

    args = parser.parse_args()
    if args.mode == "train":
        train(args)
    elif args.mode == "predict":
        predict(args)
    elif args.mode == "predict_from_csv":
        predict_from_csv(args)


if __name__ == "__main__":
    main()

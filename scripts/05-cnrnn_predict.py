#!/usr/bin/env python3
"""
Script: 05-cnrnn_predict.py
Purpose: Run trained LSTM models to predict CNVs from CORRseq segmentation calls

Usage:
    python scripts/05-cnrnn_predict.py \
        --h5_file    <smoothed_L2R.h5> \
        --calls_dir  <path/to/CORRseq_results/> \
        --del_model  <path/to/deletion_model.keras> \
        --dup_model  <path/to/duplication_model.keras> \
        --output_dir <path/to/output/> \
        [--min_markers 5] \
        [--max_markers 120] \
        [--alpha 0.05] \
        [--batch_size 256]

Inputs:
    h5_file   : HDF5 from Step 4 (smoothed_L2R, ref_qc, case_sample_ids)
    calls_dir : Directory containing CORRseq CSV files (chr*.alpha.X.XX.csv)
    del_model : Trained Keras model for deletion calls
    dup_model : Trained Keras model for duplication calls

Outputs (written to output_dir):
    Final_CNV_Predictions.csv — all calls with prediction probabilities and labels
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import h5py
import tensorflow as tf
from sklearn.preprocessing import StandardScaler
from tensorflow.keras.models import load_model
import tensorflow.keras.backend as K


# =============================================================================
# Custom Metric (required to reload trained models)
# =============================================================================

def f1_score(y_true, y_pred):
    y_pred     = K.cast(y_pred > 0.5, "float32")
    tp         = K.sum(y_true * y_pred)
    fp         = K.sum((1 - y_true) * y_pred)
    fn         = K.sum(y_true * (1 - y_pred))
    precision  = tp / (tp + fp + K.epsilon())
    recall     = tp / (tp + fn + K.epsilon())
    return 2 * (precision * recall) / (precision + recall + K.epsilon())


# =============================================================================
# Helpers
# =============================================================================

def symmetric_pad(sequences, maxlen):
    """Pad or centre-crop each sequence to a fixed length."""
    padded = []
    for seq in sequences:
        curr = len(seq)
        if curr >= maxlen:
            start = (curr - maxlen) // 2
            padded.append(seq[start: start + maxlen])
        else:
            diff  = maxlen - curr
            pad_l = diff // 2
            pad_r = diff - pad_l
            padded.append(np.pad(seq, (pad_l, pad_r), mode="constant", constant_values=0))
    return np.array(padded)[..., np.newaxis]   # shape: (N, maxlen, 1)


def natural_chrom_key(chrom):
    """Sort chromosomes in natural order (chr1 < chr2 < ... < chr10 < chrX)."""
    chrom = str(chrom).replace("chr", "")
    try:
        return (0, int(chrom))
    except ValueError:
        return (1, chrom)


# =============================================================================
# Argument Parsing
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="CN-RNN Step 5: LSTM-based CNV prediction"
    )
    parser.add_argument("--h5_file",     required=True,
                        help="HDF5 file with smoothed log2 ratio matrix (from Step 4)")
    parser.add_argument("--calls_dir",   required=True,
                        help="Directory containing CORRseq CSV files (from Step 4)")
    parser.add_argument("--del_model",   required=True,
                        help="Path to trained Keras deletion model (.keras)")
    parser.add_argument("--dup_model",   required=True,
                        help="Path to trained Keras duplication model (.keras)")
    parser.add_argument("--output_dir",  required=True,
                        help="Directory to write final predictions")
    parser.add_argument("--min_markers", type=int,   default=5,
                        help="Discard CNV calls with fewer than this many markers (default: 5)")
    parser.add_argument("--max_markers", type=int,   default=120,
                        help="Auto-accept CNV calls with more than this many markers (default: 120)")
    parser.add_argument("--alpha",       type=float, default=0.05,
                        help="Alpha level of CORRseq files to use (default: 0.05)")
    parser.add_argument("--batch_size",  type=int,   default=256,
                        help="Batch size for model inference (default: 256)")
    return parser.parse_args()


# =============================================================================
# Main
# =============================================================================

def main():
    args = parse_args()

    for path in [args.h5_file, args.calls_dir, args.del_model, args.dup_model]:
        if not os.path.exists(path):
            sys.exit(f"ERROR: Path not found: {path}")

    os.makedirs(args.output_dir, exist_ok=True)
    output_csv = os.path.join(args.output_dir, "Final_CNV_Predictions.csv")

    seq_len = 160   # fixed input length expected by the trained models

    print("=== CN-RNN | Step 5: LSTM Prediction ===")
    print(f"  H5 file    : {args.h5_file}")
    print(f"  Calls dir  : {args.calls_dir}")
    print(f"  Del model  : {args.del_model}")
    print(f"  Dup model  : {args.dup_model}")
    print(f"  Output     : {args.output_dir}")
    print()

    # -------------------------------------------------------------------------
    # 1. Load HDF5 Reference Data
    # -------------------------------------------------------------------------
    print("Loading HDF5 reference data...")
    with h5py.File(args.h5_file, "r") as hf:
        # Matrix stored as (samples × regions); transpose to (regions × samples)
        L2R_matrix = hf["smoothed_L2R"][:].T

        case_samples = [
            x.decode("utf-8") if isinstance(x, bytes) else x
            for x in hf["case_sample_ids"][:]
        ]

        ref_data = hf["ref_qc"][:]
        ref_qc   = pd.DataFrame(ref_data)
        for col in ref_qc.columns:
            if ref_qc[col].dtype == object and len(ref_qc) > 0:
                if isinstance(ref_qc[col].iloc[0], bytes):
                    ref_qc[col] = ref_qc[col].str.decode("utf-8")

    print(f"  Matrix shape : {L2R_matrix.shape}")
    print(f"  Case samples : {len(case_samples)}")

    # -------------------------------------------------------------------------
    # 2. Scale L2R Matrix
    # -------------------------------------------------------------------------
    print("Scaling L2R matrix...")
    scaler_l2r       = StandardScaler()
    L2R_matrix_scaled = scaler_l2r.fit_transform(L2R_matrix)

    # Coordinate lookup: (chrom, start) → row index
    ref_qc["coord_key"] = list(zip(ref_qc["seqnames"].astype(str), ref_qc["start"]))
    coord_to_idx = {k: i for i, k in enumerate(ref_qc["coord_key"])}

    # Mappability column name may vary
    mapp_col = "mapp" if "mapp" in ref_qc.columns else "mappability"
    gc_col   = "gc"   if "gc"   in ref_qc.columns else "GC"

    # -------------------------------------------------------------------------
    # 3. Load CORRseq Calls
    # -------------------------------------------------------------------------
    alpha_str = f"{args.alpha:.2f}"
    print(f"Loading CORRseq calls (alpha={alpha_str}) from: {args.calls_dir}")

    call_files = [
        f for f in os.listdir(args.calls_dir)
        if f.endswith(f"{alpha_str}.csv")
    ]
    if not call_files:
        sys.exit(f"ERROR: No CORRseq CSV files matching alpha={alpha_str} found in {args.calls_dir}")

    df_list = []
    for fname in call_files:
        try:
            tmp = pd.read_csv(os.path.join(args.calls_dir, fname), header=None)
            tmp.columns = ["sample_idx", "chrom", "start", "end",
                           "length_kb", "num_markers", "type"]
            df_list.append(tmp)
        except Exception as e:
            print(f"  WARNING: Skipping {fname}: {e}")

    if not df_list:
        sys.exit("ERROR: No valid CORRseq CSV files could be loaded.")

    full_df = pd.concat(df_list, ignore_index=True)
    print(f"  Total calls loaded: {len(full_df)}")

    # Map 1-based R sample index → 0-based column index and sample ID
    full_df["matrix_col_idx"] = full_df["sample_idx"] - 1
    full_df["sample_id"] = full_df["matrix_col_idx"].apply(
        lambda i: case_samples[i] if 0 <= i < len(case_samples) else "Unknown"
    )

    # -------------------------------------------------------------------------
    # 4. QC Filtering
    # -------------------------------------------------------------------------
    print("Applying QC thresholds...")

    discard_mask = full_df["num_markers"] < args.min_markers
    large_mask   = full_df["num_markers"] > args.max_markers
    rnn_mask     = (~discard_mask) & (~large_mask)

    print(f"  Discarded  (< {args.min_markers} markers)  : {discard_mask.sum()}")
    print(f"  Auto-kept  (> {args.max_markers} markers)  : {large_mask.sum()}")
    print(f"  RNN candidates                        : {rnn_mask.sum()}")

    large_cnvs = full_df[large_mask].copy()
    large_cnvs["pred_prob"]          = 1.0
    large_cnvs["final_label"]        = 1
    large_cnvs["prediction_status"]  = "Kept_Size_TooLarge"

    rnn_df = full_df[rnn_mask].reset_index(drop=True)

    if len(rnn_df) == 0:
        print("  No RNN candidates remain after marker filtering.")
        final_df = large_cnvs.copy()
        if final_df.empty:
            final_df = pd.DataFrame(columns=[
                "sample_idx", "chrom", "start", "end", "length_kb",
                "num_markers", "type", "matrix_col_idx", "sample_id",
                "pred_prob", "final_label", "prediction_status"
            ])
        final_df.to_csv(output_csv, index=False)
        print(f"  {len(final_df)} calls written to: {output_csv}")
        print("\nStep 5 complete.")
        return

    # -------------------------------------------------------------------------
    # 5. Feature Engineering
    # -------------------------------------------------------------------------
    print("Extracting features for RNN candidates...")

    meta_features = []
    sequences     = []
    valid_indices = []

    for idx, row in rnn_df.iterrows():
        start_key = (str(row["chrom"]), row["start"])
        s_idx = coord_to_idx.get(start_key)
        if s_idx is None:
            continue

        e_idx   = s_idx + int(row["num_markers"]) - 1
        col_idx = int(row["matrix_col_idx"])

        signal_slice = L2R_matrix_scaled[s_idx: e_idx + 1, col_idx]
        gc_slice     = ref_qc.iloc[s_idx: e_idx + 1][gc_col]
        map_slice    = ref_qc.iloc[s_idx: e_idx + 1][mapp_col]

        meta_features.append([
            np.log10(row["length_kb"] * 1000),
            np.log1p(row["num_markers"]),
            np.mean(signal_slice),
            np.std(signal_slice),
            np.median(gc_slice),
            np.median(map_slice),
        ])

        # Sequence with ±20 bin context window
        vec_start = max(0, s_idx - 20)
        vec_end   = min(L2R_matrix_scaled.shape[0], e_idx + 21)
        sequences.append(L2R_matrix_scaled[vec_start:vec_end, col_idx])
        valid_indices.append(idx)

    rnn_df = rnn_df.loc[valid_indices].reset_index(drop=True)

    if len(rnn_df) == 0:
        print("  No valid RNN candidates after coordinate lookup. Skipping prediction.")
    else:
        X_meta    = np.array(meta_features)
        X_seq     = symmetric_pad(sequences, maxlen=seq_len)

        meta_scaler    = StandardScaler()
        X_meta_scaled  = meta_scaler.fit_transform(X_meta)

        print(f"  Feature shapes — Seq: {X_seq.shape}, Meta: {X_meta_scaled.shape}")

        # ---------------------------------------------------------------------
        # 6. Prediction
        # ---------------------------------------------------------------------
        rnn_df["pred_prob"]         = 0.0
        rnn_df["final_label"]       = 0
        rnn_df["prediction_status"] = "RNN_Rejected"

        custom_objects = {"f1_score": f1_score}

        # Deletions
        del_mask = rnn_df["type"] == "del"
        if del_mask.sum() > 0:
            print(f"Predicting {del_mask.sum()} deletions...")
            with tf.keras.utils.custom_object_scope(custom_objects):
                model_del = load_model(args.del_model)
            probs = model_del.predict(
                [X_seq[del_mask.values], X_meta_scaled[del_mask.values]],
                batch_size=args.batch_size, verbose=1
            ).flatten()
            rnn_df.loc[del_mask, "pred_prob"]    = probs
            rnn_df.loc[del_mask, "final_label"]  = (probs > 0.5).astype(int)
            del model_del
            tf.keras.backend.clear_session()
        else:
            print("  No deletion candidates.")

        # Duplications
        dup_mask = rnn_df["type"] == "dup"
        if dup_mask.sum() > 0:
            print(f"Predicting {dup_mask.sum()} duplications...")
            with tf.keras.utils.custom_object_scope(custom_objects):
                model_dup = load_model(args.dup_model)
            probs = model_dup.predict(
                [X_seq[dup_mask.values], X_meta_scaled[dup_mask.values]],
                batch_size=args.batch_size, verbose=1
            ).flatten()
            rnn_df.loc[dup_mask, "pred_prob"]   = probs
            rnn_df.loc[dup_mask, "final_label"] = (probs > 0.5).astype(int)
            del model_dup
            tf.keras.backend.clear_session()
        else:
            print("  No duplication candidates.")

        rnn_df.loc[rnn_df["final_label"] == 1, "prediction_status"] = "RNN_Confirmed"

    # -------------------------------------------------------------------------
    # 7. Merge and Save
    # -------------------------------------------------------------------------
    print("Merging results and saving...")

    final_df = pd.concat([rnn_df, large_cnvs], ignore_index=True)

    # Natural chromosome sort
    final_df["_chrom_key"] = final_df["chrom"].apply(natural_chrom_key)
    final_df.sort_values(
        by=["sample_idx", "_chrom_key", "start"],
        inplace=True
    )
    final_df.drop(columns=["_chrom_key"], inplace=True)
    final_df.reset_index(drop=True, inplace=True)

    final_df.to_csv(output_csv, index=False)
    print(f"  {len(final_df)} calls written to: {output_csv}")
    print("\nStep 5 complete.")


if __name__ == "__main__":
    main()

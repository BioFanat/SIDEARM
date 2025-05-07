import os
import pandas as pd
import numpy as np
import argparse
import re
# Import the FDR correction function
from statsmodels.stats.multitest import fdrcorrection

def p_value_qc_adapted(cell_df, naked_df, alpha=0.05):
    """
    Analyzes and compares damage patterns between cellular and naked DNA,
    adapting the provided QC logic for specific column names and using
    the 'Comparative' mode logic. Adds 'sig_enrich' and 'sig_deplete'
    columns to the cell_df in place.

    Args:
        cell_df (pd.DataFrame): DataFrame containing cellular DNA damage data.
                                Must contain 'qval_top', 'qval_bottom', 'zscore_scal'.
        naked_df (pd.DataFrame): DataFrame containing naked DNA damage data.
                                 Must contain 'qval_top', 'qval_bottom', 'zscore_scal'.
                                 Must have the same index as cell_df.
        alpha (float): Significance threshold for q-values.

    Returns:
        None (modifies cell_df DataFrame in place)
    """
    # Ensure dataframes have the same index for proper alignment
    if not cell_df.index.equals(naked_df.index):
        # Attempt to align based on 'pos' and 'strand' if index doesn't match
        try:
            # Make sure index setting doesn't happen if columns already used as index
            if not isinstance(cell_df.index, pd.MultiIndex) or not cell_df.index.names == ['pos', 'strand']:
                 cell_df = cell_df.set_index(['pos', 'strand'], drop=False)
            if not isinstance(naked_df.index, pd.MultiIndex) or not naked_df.index.names == ['pos', 'strand']:
                 naked_df = naked_df.set_index(['pos', 'strand'], drop=False)

            # Reindex naked_df to match cell_df's index, filling missing with NaN
            naked_df = naked_df.reindex(cell_df.index)
        except KeyError:
             raise ValueError("Input DataFrames must have 'pos' and 'strand' columns if indices don't match.")
        except Exception as e:
            raise ValueError(f"Error aligning DataFrames: {e}")

    # Add new columns for significance flags
    cell_df['sig_enrich'] = False
    cell_df['sig_deplete'] = False

    # --- Vectorized Approach for Speed ---
    # Check for NaN values that could cause issues in comparisons
    cell_q_top = cell_df['qval_top'].fillna(1.0)
    cell_q_bottom = cell_df['qval_bottom'].fillna(1.0)
    naked_q_top = naked_df['qval_top'].fillna(1.0)
    naked_q_bottom = naked_df['qval_bottom'].fillna(1.0)
    cell_z = cell_df['zscore_scal'].fillna(0)
    naked_z = naked_df['zscore_scal'].fillna(0)


    # Initial significance checks
    cell_enrich_init_sig = cell_q_top < alpha
    cell_deplete_init_sig = cell_q_bottom < alpha
    naked_enrich_sig = naked_q_top < alpha
    naked_deplete_sig = naked_q_bottom < alpha

    # Z-score comparison where both are significant
    enrich_z_comp = cell_z > naked_z
    deplete_z_comp = cell_z < naked_z

    # Enrichment QC ('Comparative' logic)
    # Case 1: Significant in cell, NOT significant in naked
    sig_enrich_case1 = cell_enrich_init_sig & (~naked_enrich_sig)
    # Case 2: Significant in cell AND significant in naked -> check Z-score
    sig_enrich_case2 = cell_enrich_init_sig & naked_enrich_sig & enrich_z_comp
    # Combine cases
    cell_df.loc[sig_enrich_case1 | sig_enrich_case2, 'sig_enrich'] = True

    # Depletion QC ('Comparative' logic)
    # Case 1: Significant in cell, NOT significant in naked
    sig_deplete_case1 = cell_deplete_init_sig & (~naked_deplete_sig)
    # Case 2: Significant in cell AND significant in naked -> check Z-score
    sig_deplete_case2 = cell_deplete_init_sig & naked_deplete_sig & deplete_z_comp
    # Combine cases
    cell_df.loc[sig_deplete_case1 | sig_deplete_case2, 'sig_deplete'] = True

    # Reset index if it was set temporarily
    if 'pos' in cell_df.index.names:
        cell_df = cell_df.reset_index(drop=True)


def apply_fdr_correction(df, pval_col, qval_col, alpha):
    """Applies FDR correction if qval column is missing and pval column exists."""
    if qval_col not in df.columns:
        if pval_col in df.columns:
            print(f"    '{qval_col}' not found. Calculating from '{pval_col}' using FDR (BH method).")
            # Ensure p-values are numeric and handle NaNs before correction
            p_values = pd.to_numeric(df[pval_col], errors='coerce').fillna(1.0)
            # Handle case where all p-values might be NaN after coercion
            if p_values.isnull().all():
                 print(f"    Warning: All values in '{pval_col}' are non-numeric. Cannot calculate q-values.")
                 df[qval_col] = 1.0 # Assign default non-significant q-value
                 return False # Indicate failure
            # Perform FDR correction (Benjamini/Hochberg)
            reject, qvals = fdrcorrection(p_values, alpha=alpha, method='indep', is_sorted=False)
            df[qval_col] = qvals
            return True # Indicate success
        else:
            print(f"    Error: Cannot calculate '{qval_col}'. Missing both '{qval_col}' and '{pval_col}'.")
            return False # Indicate failure
    return True # Indicate qval column already exists or correction succeeded


def analyze_tf_data(directory, output_filename="tf_analysis_summary.xlsx", alpha=0.05):
    """
    Analyzes TF damage data from CSV files in a directory.
    Performs FDR correction if q-value columns are missing.

    Args:
        directory (str): Path to the directory containing the CSV files.
        output_filename (str): Name for the output Excel file.
        alpha (float): Significance threshold for q-values (used for QC and FDR).
    """
    results = []
    tf_names = set()
    files = os.listdir(directory)

    # Identify unique TF names
    for f in files:
        match = re.match(r"(.+?)_(cell|naked)\.csv", f)
        if match:
            tf_names.add(match.group(1))

    print(f"Found {len(tf_names)} unique TF identifiers.")

    # Process each TF
    for tf_name in sorted(list(tf_names)):
        print(f"Processing TF: {tf_name}...")
        cell_file = os.path.join(directory, f"{tf_name}_cell.csv")
        naked_file = os.path.join(directory, f"{tf_name}_naked.csv")

        try:
            # Load data
            cell_df = pd.read_csv(cell_file)
            naked_df = pd.read_csv(naked_file)

            # --- Apply FDR Correction if necessary ---
            correction_ok_cell = True
            correction_ok_naked = True

            # Correct cell data
            if not apply_fdr_correction(cell_df, 'pval_top', 'qval_top', alpha):
                correction_ok_cell = False
            if not apply_fdr_correction(cell_df, 'pval_bottom', 'qval_bottom', alpha):
                 correction_ok_cell = False

            # Correct naked data
            if not apply_fdr_correction(naked_df, 'pval_top', 'qval_top', alpha):
                 correction_ok_naked = False
            if not apply_fdr_correction(naked_df, 'pval_bottom', 'qval_bottom', alpha):
                 correction_ok_naked = False

            if not correction_ok_cell or not correction_ok_naked:
                 print(f"  Skipping {tf_name} due to missing p-value or q-value columns.")
                 continue # Skip to the next TF if correction failed

            # --- Data Validation (Post Correction) ---
            required_cols = ['strand', 'pos', 'qval_top', 'qval_bottom', 'zscore_scal']
            if not all(col in cell_df.columns for col in required_cols):
                print(f"  Skipping {tf_name}: Missing required columns in {cell_file} even after potential correction.")
                continue
            if not all(col in naked_df.columns for col in required_cols):
                 print(f"  Skipping {tf_name}: Missing required columns in {naked_file} even after potential correction.")
                 continue
            if cell_df.empty or naked_df.empty:
                 print(f"  Skipping {tf_name}: One or both input files are empty.")
                 continue


            # Apply QC - modifies cell_df in place
            p_value_qc_adapted(cell_df, naked_df, alpha=alpha)

            # --- Calculate Metrics ---
            # Filter for significant sites AFTER QC
            sig_enrich_df = cell_df[cell_df['sig_enrich']].copy()
            sig_deplete_df = cell_df[cell_df['sig_deplete']].copy()

            # Helper function for position string
            def format_pos(row):
                # Handle potential NaN or unexpected values in strand
                strand_val = row.get('strand', '')
                strand_code = 'c' if strand_val == 'Diff' else ('s' if strand_val == 'Same' else '?')
                pos_val = row.get('pos', 'N/A')
                return f"{pos_val}{strand_code}"

            # Enrichment metrics
            num_sig_enrich = len(sig_enrich_df)
            min_qval_enrich = sig_enrich_df['qval_top'].min() if num_sig_enrich > 0 else np.nan
            max_z_enrich = np.nan
            pos_max_z_enrich = ""
            if num_sig_enrich > 0:
                # Find the row with the maximum z-score among significant enrichment sites
                # Use dropna() to handle potential NaNs in zscore_scal before finding idxmax
                idx_max_z = sig_enrich_df['zscore_scal'].dropna().idxmax()
                max_z_enrich_row = sig_enrich_df.loc[idx_max_z]
                max_z_enrich = max_z_enrich_row['zscore_scal']
                pos_max_z_enrich = format_pos(max_z_enrich_row)


            # Depletion metrics
            num_sig_deplete = len(sig_deplete_df)
            min_qval_deplete = sig_deplete_df['qval_bottom'].min() if num_sig_deplete > 0 else np.nan
            min_z_deplete = np.nan
            pos_min_z_deplete = ""
            if num_sig_deplete > 0:
                 # Find the row with the minimum z-score among significant depletion sites
                 # Use dropna() to handle potential NaNs in zscore_scal before finding idxmin
                idx_min_z = sig_deplete_df['zscore_scal'].dropna().idxmin()
                min_z_deplete_row = sig_deplete_df.loc[idx_min_z]
                min_z_deplete = min_z_deplete_row['zscore_scal']
                pos_min_z_deplete = format_pos(min_z_deplete_row)


            # Store results
            results.append({
                "TF": tf_name,
                "Num Sig Sites (Enrichment)": num_sig_enrich,
                "Min Q-Value (Enrichment)": min_qval_enrich,
                "Max Z-Score (Enrichment)": max_z_enrich,
                "Position (Max Z, Enrichment)": pos_max_z_enrich,
                "Num Sig Sites (Depletion)": num_sig_deplete,
                "Min Q-Value (Depletion)": min_qval_deplete,
                "Min Z-Score (Depletion)": min_z_deplete,
                "Position (Min Z, Depletion)": pos_min_z_deplete,
            })

        except FileNotFoundError:
            print(f"  Skipping {tf_name}: Missing cell or naked file.")
        except KeyError as e:
             print(f"  Skipping {tf_name}: Missing expected column during processing - {e}")
        except Exception as e:
            print(f"  Error processing {tf_name}: {type(e).__name__} - {e}")

    # Create DataFrame and save to Excel
    if results:
        output_df = pd.DataFrame(results)
        # Reorder columns for clarity
        output_df = output_df[[
            "TF",
            "Num Sig Sites (Enrichment)",
            "Min Q-Value (Enrichment)",
            "Max Z-Score (Enrichment)",
            "Position (Max Z, Enrichment)",
            "Num Sig Sites (Depletion)",
            "Min Q-Value (Depletion)",
            "Min Z-Score (Depletion)",
            "Position (Min Z, Depletion)",
        ]]
        # Ensure output directory exists if specified in the filename path
        output_dir = os.path.dirname(output_filename)
        if output_dir and not os.path.exists(output_dir):
             os.makedirs(output_dir)
             print(f"Created output directory: {output_dir}")

        output_df.to_excel(output_filename, index=False, engine='openpyxl')
        print(f"\nAnalysis complete. Results saved to '{output_filename}'")
    else:
        print("\nNo TFs processed successfully. No output file generated.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze TF damage data from cell and naked DNA CSV files. Performs FDR correction if q-value columns are missing.")
    parser.add_argument("directory", help="Path to the directory containing the TF CSV files.")
    parser.add_argument("-o", "--output", default="tf_analysis_summary.xlsx",
                        help="Output Excel file name (default: tf_analysis_summary.xlsx)")
    parser.add_argument("-a", "--alpha", type=float, default=0.05,
                        help="Significance threshold (q-value) for analysis and FDR correction (default: 0.05)")

    args = parser.parse_args()

    if not os.path.isdir(args.directory):
        print(f"Error: Directory not found: {args.directory}")
    else:
        # Check for statsmodels installation
        try:
             __import__('statsmodels')
        except ImportError:
             print("Error: The 'statsmodels' library is required for FDR correction.")
             print("Please install it using: pip install statsmodels")
             exit(1) # Exit if dependency is missing

        analyze_tf_data(args.directory, args.output, args.alpha)
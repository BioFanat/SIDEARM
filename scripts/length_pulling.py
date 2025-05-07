import os
import csv
import sys

def get_length_from_bed(bed_file_path):
    """
    Reads the first line of a BED file, calculates length (col3 - col2),
    and returns the length or None if an error occurs.
    """
    try:
        with open(bed_file_path, 'r') as f_bed:
            first_line = f_bed.readline()

            if not first_line:
                print(f"Warning: BED file is empty: {bed_file_path}. Skipping.")
                return None

            # BED files are typically tab-separated
            parts = first_line.strip().split('\t')

            # Ensure the line has at least 3 columns (chrom, start, end)
            if len(parts) >= 3:
                try:
                    # Columns are 0-indexed, so start is parts[1], end is parts[2]
                    start_pos = int(parts[1])
                    end_pos = int(parts[2])
                    motif_length = end_pos - start_pos
                    return motif_length
                except ValueError:
                    print(f"Warning: Could not parse coordinates as integers in first line of {bed_file_path}. Line: '{first_line.strip()}'. Skipping.")
                    return None
                except Exception as e:
                    print(f"Warning: An unexpected error occurred reading coordinates from {bed_file_path}: {e}. Skipping.")
                    return None
            else:
                print(f"Warning: First line in {bed_file_path} has fewer than 3 columns. Line: '{first_line.strip()}'. Skipping.")
                return None

    except IOError as e:
        print(f"Error: Could not open or read file {bed_file_path}: {e}. Skipping.")
        return None
    except Exception as e:
        print(f"Error: An unexpected error occurred processing file {bed_file_path}: {e}. Skipping.")
        return None


def calculate_motif_lengths_plus_minus(main_directory, output_csv_file):
    """
    Processes TF cluster directories, finds motif lengths from either the
    '_plus.bed' or '_minus.bed' file, and writes results to a CSV.

    Args:
        main_directory (str): The path to the main directory containing
                              TF cluster subfolders.
        output_csv_file (str): The path where the output CSV file will be saved.
    """
    results = []
    processed_base_clusters = set() # Keep track of base clusters processed

    print(f"Scanning directory: {main_directory}")
    print("Looking for '_plus.bed' or '_minus.bed' files for length calculation.")

    # Check if the main directory exists
    if not os.path.isdir(main_directory):
        print(f"Error: Main directory not found: {main_directory}")
        return

    # Iterate through all items (files and directories) in the main directory
    for item_name in os.listdir(main_directory):
        item_path = os.path.join(main_directory, item_name)

        # Check if the item is a directory
        if os.path.isdir(item_path):

            # Determine the base cluster name, stripping potential _plus/_minus suffixes
            # from the DIRECTORY name itself. This becomes the canonical name.
            base_cluster_name = item_name
            if base_cluster_name.endswith('_plus'):
                base_cluster_name = base_cluster_name.rsplit('_plus', 1)[0]
            elif base_cluster_name.endswith('_minus'):
                 base_cluster_name = base_cluster_name.rsplit('_minus', 1)[0]

            # Skip if we've already successfully processed this base cluster via another directory
            # (e.g., processed ETS_1 via ETS_1 folder, now skipping ETS_1_plus folder)
            if base_cluster_name in processed_base_clusters:
                continue

            # Construct the expected paths to the _plus and _minus BED files
            # using the DERIVED base_cluster_name. These files are expected
            # INSIDE the current item_path directory.
            plus_bed_path = os.path.join(item_path, f"{base_cluster_name}_plus.bed")
            minus_bed_path = os.path.join(item_path, f"{base_cluster_name}_minus.bed")

            motif_length = None
            processed_successfully = False
            source_file_type = ""

            # --- Try reading the _plus.bed file first ---
            if os.path.exists(plus_bed_path):
                motif_length = get_length_from_bed(plus_bed_path)
                if motif_length is not None:
                    processed_successfully = True
                    source_file_type = "_plus.bed"

            # --- If _plus.bed didn't exist or failed, try _minus.bed ---
            if not processed_successfully and os.path.exists(minus_bed_path):
                motif_length = get_length_from_bed(minus_bed_path)
                if motif_length is not None:
                    processed_successfully = True
                    source_file_type = "_minus.bed"

            # --- Record results or warn ---
            if processed_successfully:
                print(f"  Processed: {base_cluster_name} -> Length: {motif_length} (from {source_file_type})")
                results.append([base_cluster_name, motif_length])
                processed_base_clusters.add(base_cluster_name) # Mark base cluster as done
            else:
                # Only warn if we haven't already processed this base cluster
                # And if the directory corresponding to the base name actually exists
                # (or if the current directory *is* the base name dir)
                base_dir_path = os.path.join(main_directory, base_cluster_name)
                if item_name == base_cluster_name or os.path.isdir(base_dir_path):
                     print(f"Warning: No usable BED file ('{base_cluster_name}_plus.bed' or '{base_cluster_name}_minus.bed') found or processed successfully in directory '{item_name}'.")
                # Mark as processed even on failure to avoid duplicate warnings from related folders
                processed_base_clusters.add(base_cluster_name)


    # Check if any results were found
    if not results:
        print("\nNo valid cluster BED files (_plus or _minus) found or processed.")
        return

    # Write the results to a CSV file
    try:
        with open(output_csv_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['Cluster', 'Length']) # Write header
            writer.writerows(results) # Write data
        print(f"\nSuccessfully wrote {len(results)} cluster lengths to {output_csv_file}")
    except IOError as e:
        print(f"\nError: Could not write to output CSV file {output_csv_file}: {e}")
    except Exception as e:
        print(f"\nError: An unexpected error occurred while writing the CSV file: {e}")

# --- How to use the script ---

if __name__ == "__main__":
    # --- Configuration ---
    # Option 1: Set the directory path directly in the script
    # DEFAULT_MAIN_DIRECTORY = "/path/to/your/tf_cluster_directory" # <-- *** CHANGE THIS PATH ***
    DEFAULT_MAIN_DIRECTORY = "." # Default to current directory if not specified

    # Option 2: Get the directory path from command-line arguments
    if len(sys.argv) > 1:
        main_tf_directory = sys.argv[1]
        print(f"Using directory specified in command line: {main_tf_directory}")
    else:
        main_tf_directory = DEFAULT_MAIN_DIRECTORY
        print(f"Using default/current directory: {main_tf_directory}")
        if main_tf_directory == ".":
             print(" -> You can also specify the directory as a command-line argument.")
             print(f" -> Example: python {sys.argv[0]} /path/to/your/tf_clusters")


    output_file = "tf_cluster_motif_lengths_plus_minus.csv"
    # --- End Configuration ---

    calculate_motif_lengths_plus_minus(main_tf_directory, output_file)
import os
import pandas as pd
import argparse
from scipy.stats import fisher_exact
import re

def extract_single_summary_stats(file_path):
    """
    Extract summary statistics from single format simulation results file.
    Returns a dictionary with the extracted values.
    """
    try:
        with open(file_path, 'r') as f:
            # Read first line to get TF name
            first_line = f.readline().strip()
            tf_name = first_line.split()[-2]
            
            # Read second line (column headers) to verify format
            _ = f.readline()
            
            # Read third line (headers with "Min-Top", etc)
            _ = f.readline()
            
            # Read fourth line containing the actual values
            values_line = f.readline().strip()
            values = values_line.split()
            
            # Create dictionary with extracted values
            return {
                'TF Cluster': tf_name,
                'Format': 'single',
                'Min Q (enrichment)': float(values[0]),
                '# Significant (enrichment)': int(values[1]),
                'Min Q (depletion)': float(values[2]),
                '# Significant (depletion)': int(values[3]),
                'Largest Z': float(values[4]),
                'Smallest Z': float(values[5])
            }
    except Exception as e:
        print(f"Error processing file {file_path}: {str(e)}")
        return None

def extract_comp_summary_stats(file_path):
    """
    Extract summary statistics from the comp format with contingency tables and correlations.
    Returns a dictionary with the extracted values.
    """
    try:
        with open(file_path, 'r') as f:
            content = f.read()
            
            # Try to extract TF name from the filename if it's not in the content
            tf_name = "_".join(os.path.basename(file_path).split('_')[:-1])
            
            # Extract contingency table for enrichment
            enrichment_pattern = r"Enrichment\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)"
            enrichment_match = re.search(enrichment_pattern, content)
            
            if enrichment_match:
                both_sig_enrich = int(enrichment_match.group(1))
                only_potential_enrich = int(enrichment_match.group(2))
                only_actual_enrich = int(enrichment_match.group(3))
                neither_enrich = int(enrichment_match.group(4))
                
                # Calculate Fisher's exact test for enrichment
                oddsratio_enrich, pvalue_enrich = fisher_exact([
                    [both_sig_enrich, only_potential_enrich],
                    [only_actual_enrich, neither_enrich]
                ])
                
                # Raw counts for enrichment in pipe-delimited format
                raw_counts_enrich = f"{both_sig_enrich}|{only_potential_enrich}|{only_actual_enrich}|{neither_enrich}"
            else:
                oddsratio_enrich, pvalue_enrich = None, None
                raw_counts_enrich = "NA"
            
            # Extract contingency table for depletion
            depletion_pattern = r"Depletion\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)"
            depletion_match = re.search(depletion_pattern, content)
            
            if depletion_match:
                both_sig_deplete = int(depletion_match.group(1))
                only_potential_deplete = int(depletion_match.group(2))
                only_actual_deplete = int(depletion_match.group(3))
                neither_deplete = int(depletion_match.group(4))
                
                # Calculate Fisher's exact test for depletion
                oddsratio_deplete, pvalue_deplete = fisher_exact([
                    [both_sig_deplete, only_potential_deplete],
                    [only_actual_deplete, neither_deplete]
                ])
                
                # Raw counts for depletion in pipe-delimited format
                raw_counts_deplete = f"{both_sig_deplete}|{only_potential_deplete}|{only_actual_deplete}|{neither_deplete}"
            else:
                oddsratio_deplete, pvalue_deplete = None, None
                raw_counts_deplete = "NA"
            
            # Extract correlation statistics
            correlations = {}
            
            # Pearson correlation (all positions)
            pearson_all_pattern = r"Pearson correlation \(all positions\): ([0-9.-]+)"
            pearson_all_match = re.search(pearson_all_pattern, content)
            if pearson_all_match:
                correlations['Pearson (all)'] = float(pearson_all_match.group(1))
            
            # Spearman correlation (all positions)
            spearman_all_pattern = r"Spearman correlation \(all positions\): ([0-9.-]+)"
            spearman_all_match = re.search(spearman_all_pattern, content)
            if spearman_all_match:
                correlations['Spearman (all)'] = float(spearman_all_match.group(1))
            
            # Pearson correlation (core positions)
            pearson_core_pattern = r"Pearson correlation \(core positions\): ([0-9.-]+)"
            pearson_core_match = re.search(pearson_core_pattern, content)
            if pearson_core_match:
                correlations['Pearson (core)'] = float(pearson_core_match.group(1))
            
            # Spearman correlation (core positions)
            spearman_core_pattern = r"Spearman correlation \(core positions\): ([0-9.-]+)"
            spearman_core_match = re.search(spearman_core_pattern, content)
            if spearman_core_match:
                correlations['Spearman (core)'] = float(spearman_core_match.group(1))
            
            # Create result dictionary
            result = {
                'TF Cluster': tf_name,
                'Format': 'comp',
                'Fisher P-value (enrichment)': pvalue_enrich,
                'Fisher Odds Ratio (enrichment)': oddsratio_enrich,
                'Raw Counts (enrichment)': raw_counts_enrich,
                'Fisher P-value (depletion)': pvalue_deplete,
                'Fisher Odds Ratio (depletion)': oddsratio_deplete,
                'Raw Counts (depletion)': raw_counts_deplete
            }
            
            # Add correlation stats
            result.update(correlations)
            
            return result
            
    except Exception as e:
        print(f"Error processing file {file_path}: {str(e)}")
        return None

def detect_file_format(file_path):
    """
    Detect which format the file is in.
    Returns 'single' or 'comp' as a string.
    """
    try:
        with open(file_path, 'r') as f:
            content = f.read(200)  # Read just the first bit to check
            
            if "Test\tBoth Sig." in content or "Enrichment" in content and "Depletion" in content:
                return 'comp'
            else:
                return 'single'
    except Exception:
        # Default to single if we can't determine
        return 'single'

def extract_summary_stats(file_path):
    """
    Detect file format and extract summary statistics accordingly.
    """
    file_format = detect_file_format(file_path)
    
    if file_format == 'comp':
        return extract_comp_summary_stats(file_path)
    else:
        return extract_single_summary_stats(file_path)

def process_directory(directory_path, output_file, format_type=None):
    """
    Process all files in the specified directory and compile results into an Excel file.
    
    Arguments:
        directory_path: Path to directory containing files
        output_file: Path for output Excel file
        format_type: If specified, only process files of this format ('single' or 'comp')
    """
    # List to store all results
    single_format_results = []
    comp_format_results = []
    
    # Process each file in the directory
    for filename in os.listdir(directory_path):
        file_path = os.path.join(directory_path, filename)
        
        # Skip if not a file
        if not os.path.isfile(file_path):
            continue
        
        # Skip hidden files
        if filename.startswith('.'):
            continue
            
        # Detect file format before extracting
        detected_format = detect_file_format(file_path)
        
        # Skip if format_type is specified and doesn't match
        if format_type and detected_format != format_type:
            continue
            
        # Extract stats from file
        stats = extract_summary_stats(file_path)
        if stats:
            if stats['Format'] == 'single':
                single_format_results.append(stats)
            else:
                comp_format_results.append(stats)
    
    # Write results to separate sheets in Excel
    with pd.ExcelWriter(output_file) as writer:
        if single_format_results:
            df_single = pd.DataFrame(single_format_results)
            # Ensure columns are in the correct order for single format
            column_order_single = [
                'TF Cluster',
                'Min Q (enrichment)',
                '# Significant (enrichment)',
                'Largest Z',
                'Min Q (depletion)',
                '# Significant (depletion)',
                'Smallest Z'
            ]
            # Filter to only include columns that exist
            column_order_single = [col for col in column_order_single if col in df_single.columns]
            df_single = df_single[column_order_single]
            df_single.to_excel(writer, sheet_name='Single Format', index=False)
        
        if comp_format_results:
            df_comp = pd.DataFrame(comp_format_results)
            # Ensure columns are in the correct order for comp format
            column_order_comp = [
                'TF Cluster',
                'Fisher P-value (enrichment)',
                'Fisher Odds Ratio (enrichment)',
                'Raw Counts (enrichment)',
                'Fisher P-value (depletion)',
                'Fisher Odds Ratio (depletion)',
                'Raw Counts (depletion)',
                'Pearson (all)',
                'Spearman (all)',
                'Pearson (core)',
                'Spearman (core)'
            ]
            # Filter to only include columns that exist
            column_order_comp = [col for col in column_order_comp if col in df_comp.columns]
            df_comp = df_comp[column_order_comp]
            df_comp.to_excel(writer, sheet_name='Comp Format', index=False)
    
    print(f"Results successfully written to {output_file}")
    print(f"Processed {len(single_format_results)} single format files and {len(comp_format_results)} comp format files.")

def generate_comp_format_example(output_file):
    """
    Generate an example of the comp format summary statistics file.
    """
    example_content = """Test\tBoth Sig.\tOnly Potential\tOnly Actual\tNeither
Enrichment\t5\t0\t0\t20
Depletion\t0\t0\t1\t18
Pearson correlation (all positions): 0.9290
Spearman correlation (all positions): 0.6695
Pearson correlation (core positions): 0.9104
Spearman correlation (core positions): 0.6727
"""
    with open(output_file, 'w') as f:
        f.write(example_content)
    
    print(f"Example file created at {output_file}")

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Process simulation results from a directory into an Excel file.')
    parser.add_argument('input_dir', help='Directory containing simulation result files')
    parser.add_argument('output_file', help='Path for the output Excel file (e.g., results.xlsx)')
    parser.add_argument('--format', choices=['single', 'comp'], help='Only process files of this format')
    parser.add_argument('--generate-example', action='store_true', help='Generate an example of the comp format file')
    parser.add_argument('--example-output', help='Path for the example output file')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Generate example if requested
    if args.generate_example:
        output_path = args.example_output if args.example_output else 'example_comp_format.txt'
        generate_comp_format_example(output_path)
        return
    
    # Check if directory exists
    if not os.path.isdir(args.input_dir):
        print(f"Error: Directory '{args.input_dir}' does not exist!")
        return
    
    # Process the directory
    process_directory(args.input_dir, args.output_file, args.format)

if __name__ == "__main__":
    main()
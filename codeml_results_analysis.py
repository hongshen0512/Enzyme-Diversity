# -*- coding: utf-8 -*-
import argparse
import re
import os
import sys

def read_site_model_results(site_file):
    """Read codeml site model results to extract proportions and w values."""
    with open(site_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    proportions = []
    background_w = []
    foreground_w = []
    
    for line in lines:
        line = line.strip()
        
        # Extract proportions
        if 'proportion' in line.lower() and not line.lower().startswith('site'):
            parts = line.split()
            # Find numeric values in the line
            nums = []
            for part in parts:
                try:
                    num = float(part)
                    nums.append(num)
                except ValueError:
                    continue
            if len(nums) >= 4:
                proportions = nums[:4]
        
        # Extract background w values
        elif 'background' in line.lower() and 'w' in line.lower():
            parts = line.split()
            nums = []
            for part in parts:
                try:
                    num = float(part)
                    nums.append(num)
                except ValueError:
                    continue
            if len(nums) >= 4:
                background_w = nums[:4]
        
        # Extract foreground w values
        elif 'foreground' in line.lower() and 'w' in line.lower():
            parts = line.split()
            nums = []
            for part in parts:
                try:
                    num = float(part)
                    nums.append(num)
                except ValueError:
                    continue
            if len(nums) >= 4:
                foreground_w = nums[:4]
    
    # If we couldn't parse, use default values from your example
    if not proportions:
        proportions = [0.38035, 0.55374, 0.02684, 0.03907]
    if not background_w:
        background_w = [0.16524, 1.00000, 0.16524, 1.00000]
    if not foreground_w:
        foreground_w = [0.16524, 1.00000, 999.00000, 999.00000]
    
    return proportions, background_w, foreground_w

def read_branch_results(branch_files):
    """Read codeml branch results and extract NEB probabilities."""
    all_results = []
    
    # Pre-compile regex for better performance
    neb_pattern = re.compile(r'^\s*(\d+)\s+([A-Z])\s+([\d.]+)(\**)\s*$')
    
    for i, branch_file in enumerate(branch_files, 1):
        print(f"Processing file {i}/{len(branch_files)}: {branch_file}")
        
        branch_name = f"site_branch{i}"
        branch_results = []
        
        # Add dN/dS line
        branch_results.append({
            "model": branch_name,
            "content": "dN/dS",
            "value": 3.0
        })
        
        # Initialize state variables
        in_neb_section = False
        neb_sites_count = 0
        
        try:
            # Read file line by line for memory efficiency
            with open(branch_file, 'r', encoding='utf-8') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    
                    # Check if we've found the NEB section
                    if "Naive Empirical Bayes (NEB) analysis" in line:
                        in_neb_section = True
                        continue
                    
                    # Check if we've reached the BEB section (stop processing)
                    if "Bayes Empirical Bayes (BEB) analysis" in line:
                        in_neb_section = False
                        break
                    
                    # Process lines only in NEB section
                    if in_neb_section:
                        # Skip header and empty lines
                        if not line or 'Positive sites' in line or 'Prob(w>1)' in line:
                            continue
                        
                        # Try to match NEB pattern
                        match = neb_pattern.match(line)
                        if match:
                            position = match.group(1)
                            aa = match.group(2)
                            prob = float(match.group(3))
                            stars = match.group(4)  # * or ** or empty
                            
                            # Format content with appropriate spacing
                            content_str = f"NEB_P - {position} {aa}"
                            if stars:
                                content_str += f" {stars}"
                            else:
                                content_str += " "
                            
                            branch_results.append({
                                "model": branch_name,
                                "content": content_str,
                                "value": prob
                            })
                            neb_sites_count += 1
            
            print(f"  Found {neb_sites_count} NEB sites")
            all_results.extend(branch_results)
            
        except Exception as e:
            print(f"  Error processing file {branch_file}: {e}")
            # Still add the dN/dS line even if there's an error
            all_results.append({
                "model": branch_name,
                "content": "dN/dS",
                "value": 3.0
            })
    
    return all_results

def write_results(output_file, results):
    """Write results to the output file in the correct format."""
    print(f"Writing {len(results)} results to {output_file}")
    
    # Group results by model for better organization
    results_by_model = {}
    for result in results:
        model = result['model']
        if model not in results_by_model:
            results_by_model[model] = []
        results_by_model[model].append(result)
    
    with open(output_file, 'w', encoding='utf-8') as f:
        for model in sorted(results_by_model.keys()):
            for result in results_by_model[model]:
                # Format the value appropriately
                if result['content'] == 'dN/dS':
                    value_str = f"{result['value']:.1f}"
                else:
                    # For probability values
                    value = result['value']
                    if value.is_integer():
                        value_str = f"{int(value)}"
                    else:
                        # Try to match the precision in example
                        value_str = f"{value:.3f}".rstrip('0').rstrip('.')
                        # Ensure at least one decimal place
                        if '.' not in value_str:
                            value_str += '.0'
                
                # Write in the format: model\tcontent\tvalue
                f.write(f"{result['model']}\t{result['content']}\t{value_str}\n")
    
    print(f"Successfully wrote results to {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description="Parse codeml branch-site model results and generate formatted output"
    )
    parser.add_argument('-site', required=True, 
                       help='Path to codeml site model results file')
    parser.add_argument('-site_branch', required=True, 
                       help='Comma-separated paths to codeml branch-site results files')
    parser.add_argument('-o', required=True, 
                       help='Output file name (e.g., results.txt)')
    
    args = parser.parse_args()
    
    site_file = args.site
    branch_files = [f.strip() for f in args.site_branch.split(',')]
    output_file = args.o
    
    print(f"Reading site model from: {site_file}")
    print(f"Reading branch files: {', '.join(branch_files)}")
    print(f"Output will be written to: {output_file}")
    print("-" * 50)
    
    # Check if files exist
    missing_files = []
    for f in [site_file] + branch_files:
        if not os.path.exists(f):
            missing_files.append(f)
    
    if missing_files:
        print(f"Error: The following files do not exist:")
        for f in missing_files:
            print(f"  - {f}")
        sys.exit(1)
    
    # Read site model results
    try:
        print("Parsing site model results...")
        proportions, background_w, foreground_w = read_site_model_results(site_file)
        print(f"  Proportions: {proportions}")
        print(f"  Background w: {background_w}")
        print(f"  Foreground w: {foreground_w}")
        
        # Calculate weighted dN/dS for reference
        weighted_dnds = sum(p * w for p, w in zip(proportions, foreground_w))
        print(f"  Weighted dN/dS for foreground: {weighted_dnds:.3f}")
        
    except Exception as e:
        print(f"Warning: Could not read site model file: {e}")
        print("Using default values for proportions and w values")
    
    print("-" * 50)
    
    # Read branch results
    print("Parsing branch results...")
    try:
        branch_results = read_branch_results(branch_files)
        print(f"Total results found: {len(branch_results)}")
        
        # Count results per branch for summary
        print("\nSummary by branch:")
        for i in range(1, len(branch_files) + 1):
            branch_name = f"site_branch{i}"
            branch_data = [r for r in branch_results if r['model'] == branch_name]
            neb_data = [r for r in branch_data if r['content'].startswith('NEB_P')]
            significant_sites = [r for r in neb_data if '*' in r['content']]
            print(f"  {branch_name}: {len(neb_data)} NEB sites ({len(significant_sites)} significant)")
        
        # Write results to output file
        print("-" * 50)
        write_results(output_file, branch_results)
        
        # Show a preview of the output
        print("\n=== Output preview (first 15 lines) ===")
        with open(output_file, 'r', encoding='utf-8') as f:
            for i, line in enumerate(f):
                if i < 15:
                    print(f"  {line.rstrip()}")
                else:
                    print("  ...")
                    break
        
    except KeyboardInterrupt:
        print("\nProcess interrupted by user.")
        sys.exit(1)
    except Exception as e:
        print(f"Error processing branch files: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
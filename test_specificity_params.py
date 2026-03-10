#!/usr/bin/env python3
"""
Parameter sweep for specificity thresholds.
Tests 3 bispecific pairs across a 4x4 matrix of parameters.
"""
import subprocess
import json
import csv
from pathlib import Path
from datetime import datetime

# Parameter matrix
MAX_NONSPECIFIC_VALUES = [2, 3, 4, 5]
SPECIFICITY_THRESHOLDS = [0.70, 0.75, 0.80, 0.85]
TEST_PAIRS = [
    "CEACAM5:CLDN18",
    "ERBB2:NECTIN4",
    "ITGB6:ADAM8"
]

# Output
RESULTS_DIR = Path("specificity_sweep_results")
RESULTS_DIR.mkdir(exist_ok=True)
timestamp = datetime.now().strftime("%Y%m%d_%H%M")
results_csv = RESULTS_DIR / f"param_sweep_{timestamp}.csv"

print("=" * 80)
print("SPECIFICITY PARAMETER SWEEP")
print("=" * 80)
print(f"Testing {len(MAX_NONSPECIFIC_VALUES)} × {len(SPECIFICITY_THRESHOLDS)} = "
      f"{len(MAX_NONSPECIFIC_VALUES) * len(SPECIFICITY_THRESHOLDS)} combinations")
print(f"Pairs: {', '.join(TEST_PAIRS)}")
print(f"Max non-specific residues: {MAX_NONSPECIFIC_VALUES}")
print(f"BLAST identity thresholds: {[f'{t:.0%}' for t in SPECIFICITY_THRESHOLDS]}")
print(f"Results: {results_csv}")
print("=" * 80)
print()

# Store results
results = []

# Save original config values
config_path = Path("config.py")
with open(config_path) as f:
    original_config = f.read()

try:
    for i, max_nonspec in enumerate(MAX_NONSPECIFIC_VALUES):
        for j, spec_threshold in enumerate(SPECIFICITY_THRESHOLDS):
            run_num = i * len(SPECIFICITY_THRESHOLDS) + j + 1
            total_runs = len(MAX_NONSPECIFIC_VALUES) * len(SPECIFICITY_THRESHOLDS)

            print(f"[{run_num}/{total_runs}] Testing: max_nonspecific={max_nonspec}, "
                  f"specificity_threshold={spec_threshold:.0%}")

            # Modify config.py
            modified_config = original_config
            modified_config = modified_config.replace(
                "MAX_NONSPECIFIC_PER_600A2 = 4",
                f"MAX_NONSPECIFIC_PER_600A2 = {max_nonspec}"
            )
            modified_config = modified_config.replace(
                "SPECIFICITY_IDENTITY_THRESHOLD = 0.75",
                f"SPECIFICITY_IDENTITY_THRESHOLD = {spec_threshold}"
            )

            with open(config_path, 'w') as f:
                f.write(modified_config)

            # Run bispecific pipeline
            try:
                cmd = [
                    "python", "-m", "epitope_pipeline.bispecific",
                    *TEST_PAIRS
                ]
                result = subprocess.run(
                    cmd,
                    cwd="..",  # Run from parent directory
                    capture_output=True,
                    text=True,
                    timeout=300,
                    check=True
                )

                # Parse results from the last run directory
                runs_dir = Path("runs")
                latest_run = max(runs_dir.glob("*bispecific*"), key=lambda p: p.stat().st_mtime)
                pairs_csv = latest_run / "bispecific_pairs.csv"

                if pairs_csv.exists():
                    with open(pairs_csv) as f:
                        reader = csv.DictReader(f)
                        pair_results = list(reader)

                    # Extract scores for our 3 test pairs
                    for pair in TEST_PAIRS:
                        pair_data = next((r for r in pair_results if r['pair'] == pair), None)
                        if pair_data:
                            results.append({
                                'max_nonspecific': max_nonspec,
                                'specificity_threshold': spec_threshold,
                                'pair': pair,
                                'pair_score': float(pair_data['pair_score']),
                                'status': pair_data.get('status', ''),
                                'run_dir': latest_run.name
                            })
                            print(f"  {pair}: {float(pair_data['pair_score']):.3f}")
                else:
                    print(f"  WARNING: Results CSV not found at {pairs_csv}")

            except subprocess.TimeoutExpired:
                print(f"  ERROR: Timeout after 300s")
            except subprocess.CalledProcessError as e:
                print(f"  ERROR: Pipeline failed - {e}")
            except Exception as e:
                print(f"  ERROR: {e}")

            print()

finally:
    # Restore original config
    print("Restoring original config.py...")
    with open(config_path, 'w') as f:
        f.write(original_config)

# Write results to CSV
print("=" * 80)
print("Writing results...")
with open(results_csv, 'w', newline='') as f:
    if results:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
        print(f"✓ Results written to: {results_csv}")
    else:
        print("! No results to write")

# Print summary table
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
for pair in TEST_PAIRS:
    print(f"\n{pair}:")
    print(f"  {'Max Non-Spec':<15} {'70%':<8} {'75%':<8} {'80%':<8} {'85%':<8}")
    print("  " + "-" * 47)
    for max_nonspec in MAX_NONSPECIFIC_VALUES:
        row = [f"  {max_nonspec} residues"]
        for spec_threshold in SPECIFICITY_THRESHOLDS:
            pair_result = next(
                (r for r in results if r['pair'] == pair and
                 r['max_nonspecific'] == max_nonspec and
                 r['specificity_threshold'] == spec_threshold),
                None
            )
            if pair_result:
                score = pair_result['pair_score']
                row.append(f"{score:.3f}" if score > 0 else "0.000")
            else:
                row.append("N/A")
        print(f"{row[0]:<15} {row[1]:<8} {row[2]:<8} {row[3]:<8} {row[4]:<8}")

print("\n" + "=" * 80)
print("PARAMETER SWEEP COMPLETE")
print("=" * 80)

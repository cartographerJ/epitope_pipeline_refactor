#!/usr/bin/env python3
"""Run 2D parameter sweep on ADAM8:ITGB6."""

import logging
import sys
from datetime import datetime
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from epitope_pipeline.utils import setup_logging
from epitope_pipeline.diagnostics.bispecific_parameter_sweep_2d import run_bispecific_2d_sweep


def main():
    """Run 2D parameter sweep on ADAM8:ITGB6."""
    timestamp = datetime.now().strftime("%y%m%d_%H%M")
    run_name = f"{timestamp}_param_sweep_2d_adam8_itgb6"
    run_dir = Path(__file__).parent.parent / "runs" / run_name
    run_dir.mkdir(parents=True, exist_ok=True)

    # Setup logging
    log_file = run_dir / "parameter_sweep_2d.log"
    setup_logging(verbose=True)

    # Add file handler
    logger = logging.getLogger()
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    ))
    logger.addHandler(file_handler)

    # Run 2D sweep
    pair_name = 'ADAM8:ITGB6'
    run_bispecific_2d_sweep(pair_name, run_dir)


if __name__ == '__main__':
    main()

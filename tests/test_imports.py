"""Smoke test: every public module imports cleanly.

Catches broken imports introduced by file moves and refactors. Cheap to run
(~0.5 s) — run after every reorg commit.
"""
import importlib

PUBLIC_MODULES = [
    "epitope_pipeline",
    "epitope_pipeline.config",
    "epitope_pipeline.utils",
    "epitope_pipeline.run",
    "epitope_pipeline.bispecific",
    "epitope_pipeline.spatial",
    "epitope_pipeline.surface",
    "epitope_pipeline.conservation",
    "epitope_pipeline.specificity",
    "epitope_pipeline.scoring",
    "epitope_pipeline.target_input",
    "epitope_pipeline.structure",
    "epitope_pipeline.membrane",
    "epitope_pipeline.export",
    "epitope_pipeline.export_bispecific",
    "epitope_pipeline.visualize",
    "epitope_pipeline.visualize_bispecific",
]


def test_all_modules_import():
    for name in PUBLIC_MODULES:
        importlib.import_module(name)

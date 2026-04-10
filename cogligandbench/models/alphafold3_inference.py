"""AlphaFold3 docking wrapper.

Runs AF3 in no-MSA / single-sequence mode via a subprocess in a dedicated
conda env, then extracts ligand SDFs from the predicted mmCIFs using the
known input SMILES for bond-order recovery.

Public surface (used by cogligandbench.engine):
  - run_single(protein, ligand, output_dir, config=None, prefix=None, **kwargs)
  - run_dataset(config)
"""

from __future__ import annotations

import os
from typing import Optional


def run_single(
    protein: str,
    ligand: str,
    output_dir: str,
    config: Optional[dict] = None,
    prefix: Optional[str] = None,
    **kwargs,
) -> str:
    """Dock a single protein+ligand pair with AlphaFold3.

    Implementation lands in a later task; this stub exists so the engine
    registration tests pass.
    """
    raise NotImplementedError("alphafold3 run_single not implemented yet")


def run_dataset(config: dict) -> None:
    """Run AlphaFold3 over an entire dataset directory.

    Implementation lands in a later task; this stub exists so the engine
    registration tests pass.
    """
    raise NotImplementedError("alphafold3 run_dataset not implemented yet")

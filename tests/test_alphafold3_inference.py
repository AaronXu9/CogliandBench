"""Unit tests for the AlphaFold3 wrapper helpers (subprocess-free)."""

import json
import os

import pytest

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIXTURE_DIR = os.path.join(PROJECT_ROOT, "data", "runsNposes", "8gkf__1__1.A__1.J")
FIXTURE_PROTEIN = os.path.join(FIXTURE_DIR, "8gkf__1__1.A__1.J_protein.pdb")
FIXTURE_LIGAND = os.path.join(FIXTURE_DIR, "8gkf__1__1.A__1.J_ligand.sdf")


def _require_fixture():
    if not (os.path.exists(FIXTURE_PROTEIN) and os.path.exists(FIXTURE_LIGAND)):
        pytest.skip(f"Fixture not found at {FIXTURE_DIR}")


class TestSmilesFromSdf:
    def test_returns_canonical_smiles(self):
        _require_fixture()
        from rdkit import Chem
        from cogligandbench.models.alphafold3_inference import _smiles_from_sdf

        smiles = _smiles_from_sdf(FIXTURE_LIGAND)
        assert isinstance(smiles, str)
        assert len(smiles) > 0
        # Canonicalization should be a fixed point
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None
        assert Chem.MolToSmiles(mol) == smiles


class TestBuildAf3InputJson:
    def test_structure_has_required_top_level_keys(self):
        _require_fixture()
        from cogligandbench.models.alphafold3_inference import _build_af3_input_json

        d = _build_af3_input_json("8gkf_test", FIXTURE_PROTEIN, FIXTURE_LIGAND)
        assert d["dialect"] == "alphafold3"
        assert d["version"] == 2
        assert d["name"] == "8gkf_test"
        assert isinstance(d["modelSeeds"], list)
        assert len(d["modelSeeds"]) >= 1
        assert "sequences" in d and isinstance(d["sequences"], list)

    def test_sequences_have_one_protein_per_chain_and_one_ligand(self):
        _require_fixture()
        from cogligandbench.models.alphafold3_inference import _build_af3_input_json

        d = _build_af3_input_json("8gkf_test", FIXTURE_PROTEIN, FIXTURE_LIGAND)
        proteins = [s for s in d["sequences"] if "protein" in s]
        ligands = [s for s in d["sequences"] if "ligand" in s]
        assert len(proteins) >= 1
        assert len(ligands) == 1
        for p in proteins:
            assert p["protein"]["sequence"]
            assert p["protein"]["id"]                      # chain id, e.g. "A"
        assert ligands[0]["ligand"]["id"] == "L"
        assert ligands[0]["ligand"]["smiles"]

    def test_protein_entries_use_single_sequence_a3m(self):
        """The no-MSA invariant: every protein has unpairedMsa = '>query\\n{seq}\\n'."""
        _require_fixture()
        from cogligandbench.models.alphafold3_inference import _build_af3_input_json

        d = _build_af3_input_json("8gkf_test", FIXTURE_PROTEIN, FIXTURE_LIGAND)
        for entry in d["sequences"]:
            if "protein" not in entry:
                continue
            seq = entry["protein"]["sequence"]
            assert entry["protein"]["unpairedMsa"] == f">query\n{seq}\n"
            assert entry["protein"]["pairedMsa"] == ""
            assert entry["protein"]["templates"] == []

    def test_json_round_trips(self):
        """The output must be JSON-serializable as-is."""
        _require_fixture()
        from cogligandbench.models.alphafold3_inference import _build_af3_input_json

        d = _build_af3_input_json("8gkf_test", FIXTURE_PROTEIN, FIXTURE_LIGAND)
        encoded = json.dumps(d)
        decoded = json.loads(encoded)
        assert decoded == d

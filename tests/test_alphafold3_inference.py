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

    def test_raises_value_error_for_missing_file(self, tmp_path):
        from cogligandbench.models.alphafold3_inference import _smiles_from_sdf

        missing = tmp_path / "does_not_exist.sdf"
        with pytest.raises(ValueError, match="Cannot open SDF file"):
            _smiles_from_sdf(str(missing))


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


def _write_tiny_cif_from_sdf(
    sdf_path, cif_path, chain_id="L", resname="LIG", use_auth_asym_id=True
):
    """Build a minimal mmCIF containing chain ``chain_id`` HETATM records
    with the heavy-atom coordinates from the first molecule in ``sdf_path``.

    When ``use_auth_asym_id`` is True (default), the file includes BOTH
    ``_atom_site.auth_asym_id`` and ``_atom_site.label_asym_id`` columns.
    When False, only ``_atom_site.label_asym_id`` is written, exercising the
    fallback code path in ``_extract_ligand_from_cif``.
    """
    from rdkit import Chem

    mol = next(Chem.SDMolSupplier(str(sdf_path), removeHs=True), None)
    assert mol is not None, f"Cannot read SDF: {sdf_path}"
    conf = mol.GetConformer()

    lines = [
        "data_test",
        "loop_",
        "_atom_site.group_PDB",
        "_atom_site.id",
        "_atom_site.type_symbol",
        "_atom_site.label_atom_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_seq_id",
    ]
    if use_auth_asym_id:
        lines.append("_atom_site.auth_asym_id")
        lines.append("_atom_site.auth_seq_id")
    lines += [
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.occupancy",
        "_atom_site.B_iso_or_equiv",
    ]
    for i, atom in enumerate(mol.GetAtoms(), start=1):
        pos = conf.GetAtomPosition(atom.GetIdx())
        elem = atom.GetSymbol()
        # Atom name must be unique within the residue; use element + index
        atom_name = f"{elem}{i}"
        if use_auth_asym_id:
            lines.append(
                f"HETATM {i} {elem} {atom_name} {resname} {chain_id} 1 "
                f"{chain_id} 1 "
                f"{pos.x:.3f} {pos.y:.3f} {pos.z:.3f} 1.00 0.00"
            )
        else:
            lines.append(
                f"HETATM {i} {elem} {atom_name} {resname} {chain_id} 1 "
                f"{pos.x:.3f} {pos.y:.3f} {pos.z:.3f} 1.00 0.00"
            )
    cif_path.write_text("\n".join(lines) + "\n")


class TestExtractLigandFromCif:
    def test_recovers_canonical_smiles_from_synthetic_cif(self, tmp_path):
        _require_fixture()
        from rdkit import Chem
        from cogligandbench.models.alphafold3_inference import (
            _extract_ligand_from_cif, _smiles_from_sdf,
        )

        cif_path = tmp_path / "tiny.cif"
        _write_tiny_cif_from_sdf(FIXTURE_LIGAND, cif_path)

        template_smiles = _smiles_from_sdf(FIXTURE_LIGAND)
        template = Chem.MolFromSmiles(template_smiles)

        mol = _extract_ligand_from_cif(cif_path, template_mol=template)
        assert mol is not None
        # After bond-order recovery, the canonical SMILES should match
        assert Chem.MolToSmiles(Chem.RemoveHs(mol)) == template_smiles

    def test_raises_when_chain_l_missing(self, tmp_path):
        _require_fixture()
        from rdkit import Chem
        from cogligandbench.models.alphafold3_inference import (
            _extract_ligand_from_cif, _smiles_from_sdf,
        )

        # Write a CIF where the ligand is on chain "Z" instead of "L"
        cif_path = tmp_path / "wrong_chain.cif"
        _write_tiny_cif_from_sdf(FIXTURE_LIGAND, cif_path, chain_id="Z")

        template = Chem.MolFromSmiles(_smiles_from_sdf(FIXTURE_LIGAND))
        with pytest.raises(ValueError, match="No ligand atoms"):
            _extract_ligand_from_cif(cif_path, template_mol=template)

    def test_uses_label_asym_id_when_auth_asym_id_missing(self, tmp_path):
        """Verify the fallback to label_asym_id when auth_asym_id column is absent."""
        _require_fixture()
        from rdkit import Chem
        from cogligandbench.models.alphafold3_inference import (
            _extract_ligand_from_cif, _smiles_from_sdf,
        )

        cif_path = tmp_path / "label_only.cif"
        _write_tiny_cif_from_sdf(FIXTURE_LIGAND, cif_path, use_auth_asym_id=False)

        template = Chem.MolFromSmiles(_smiles_from_sdf(FIXTURE_LIGAND))
        mol = _extract_ligand_from_cif(cif_path, template_mol=template)
        assert mol is not None
        assert Chem.MolToSmiles(Chem.RemoveHs(mol)) == _smiles_from_sdf(FIXTURE_LIGAND)

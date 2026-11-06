import copy
import json
from pathlib import Path
import unittest

from polychemmd.records import MaterialSpec

EXAMPLE = Path(__file__).resolve().parents[1] / 'examples' / 'pe_record.json'


class RecordTests(unittest.TestCase):
    def setUp(self):
        self.data = json.loads(EXAMPLE.read_text())

    def test_roundtrip_preserves_uninterpreted_smiles(self):
        record = MaterialSpec.from_mapping(self.data)
        self.assertEqual(record.to_dict(), self.data)
        # A graph toolkit has not checked valence: the API preserves this fact.
        self.assertEqual(record.repeat_units[0].smiles, '[*:1]CC[*:2]')

    def test_invalid_state_and_numeric_values(self):
        for field, value in [('temperature_k', 0), ('temperature_k', float('nan')),
                             ('temperature_k', True), ('temperature_k', 10**400),
                             ('pressure_bar', float('inf'))]:
            with self.subTest(field=field, value=value):
                data = copy.deepcopy(self.data)
                data['state'][field] = value
                with self.assertRaises(ValueError):
                    MaterialSpec.from_mapping(data)

    def test_fractions_cannot_be_renormalized_silently(self):
        self.data['repeat_units'][0]['mole_fraction'] = .8
        with self.assertRaisesRegex(ValueError, 'sum to 1'):
            MaterialSpec.from_mapping(self.data)

    def test_unknown_keys_and_schema_version(self):
        self.data['density'] = .91
        with self.assertRaisesRegex(ValueError, 'unknown keys'):
            MaterialSpec.from_mapping(self.data)
        del self.data['density']
        self.data['schema_version'] = '9'
        with self.assertRaisesRegex(ValueError, 'schema_version'):
            MaterialSpec.from_mapping(self.data)

    def test_sequence_and_composition_consistency(self):
        self.data['repeat_units'][0]['mole_fraction'] = .5
        self.data['repeat_units'].append({'label': 'B', 'smiles': '[*:1]CO[*:2]', 'mole_fraction': .5})
        with self.assertRaisesRegex(ValueError, 'homopolymer'):
            MaterialSpec.from_mapping(self.data)
        self.data['microstructure']['sequence_type'] = 'alternating'
        self.assertEqual(len(MaterialSpec.from_mapping(self.data).repeat_units), 2)

    def test_duplicate_labels_and_boolean_dp(self):
        self.data['microstructure']['degree_of_polymerization'] = True
        with self.assertRaisesRegex(ValueError, 'positive integer'):
            MaterialSpec.from_mapping(self.data)
        self.data['microstructure']['degree_of_polymerization'] = 10
        self.data['repeat_units'][0]['mole_fraction'] = .5
        self.data['repeat_units'].append(dict(self.data['repeat_units'][0]))
        with self.assertRaisesRegex(ValueError, 'unique'):
            MaterialSpec.from_mapping(self.data)

"""Explicit storage bounds, preserving legacy scientific recipes exactly."""

import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from polychemmd.builders.packing import pack_cell
from polychemmd.engines.periodic_alkane import _neighbor_capacity_commands, build_periodic_alkane_setup
from polychemmd.engines.periodic_thermo import plan_periodic_alkane_thermo_singlepoint
from polychemmd.forcefields import assign_compact_periodic_alkanes, assign_periodic_alkanes
from polychemmd.protocols.alkane import PeriodicAlkaneProtocolSpec, plan_periodic_alkane_protocol
from polychemmd.protocols.restart import (
    PeriodicAlkaneContinuationSpec, _digest_record, plan_periodic_alkane_continuation,
)
from test_periodic_assignment import analytic_ethane


ROOT=Path(__file__).resolve().parents[1]


def synthetic_checkpoint(assigned):
    """Non-executable manifest fixture for plan-only schema checks."""
    record=assigned.to_dict()
    # Preserve the established complete checkpoint schema, replace its identity
    # with this strictly assigned source and label the producer as a fixture.
    checkpoint=json.loads((ROOT/'evidence/restart_v06/preparation/protocol.json').read_text())['checkpoint']
    checkpoint.update(assignment_canonical_sha256=hashlib.sha256(
        json.dumps(record,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
        model_id=record['model_id'],atom_count=len(record['atoms']),
        box={key:record['box'][key] for key in ('origin_angstrom','lengths_angstrom','periodic')},
        producer_engine_identity={'synthetic_manifest_fixture':True},engine_version_banner='NOT_AN_EXECUTABLE_CHECKPOINT')
    checkpoint['checkpoint_canonical_sha256']=_digest_record(checkpoint)
    return checkpoint


class NeighborCapacityPolicyTests(unittest.TestCase):
    def test_size_derived_storage_is_finite_bounded_and_page_is_sufficient(self):
        expected={513:(2000,100000),4704:(9408,100000),8192:(16384,163840)}
        for count,(one,page) in expected.items():
            with self.subTest(count=count):
                commands=_neighbor_capacity_commands(count)
                self.assertEqual(commands,(f'neigh_modify one {one} page {page}',))
                self.assertGreaterEqual(one,2*count)
                self.assertGreaterEqual(page,10*one)
                self.assertLessEqual(one,16384)
                self.assertLessEqual(page,163840)

    def test_no_commands_for_legacy_size_and_invalid_counts_rejected(self):
        for count in (1,128,512):self.assertEqual(_neighbor_capacity_commands(count),())
        for count in (True,False,0,-1,513.,8193,float('nan'),'4704',10**100):
            with self.subTest(count=str(count)),self.assertRaises(ValueError):
                _neighbor_capacity_commands(count)

    def test_existing_small_legacy_setup_and_data_hashes_remain_exact(self):
        cell=json.loads((ROOT/'evidence/geometry_v03/cell.json').read_text())['payload']
        legacy=assign_periodic_alkanes(cell);compact=assign_compact_periodic_alkanes(cell)
        setup,data=build_periodic_alkane_setup(legacy)
        self.assertEqual(hashlib.sha256(setup.encode()).hexdigest(),
                         '3c564cd58067e11357db530d398ae1a5ea809a130a3b5f5546020ddd30043b2b')
        self.assertEqual(hashlib.sha256(data.encode()).hexdigest(),
                         '977b46d38455bb11197dbd2b8c876fb95724e5d0c877711e7c8c94be70fae876')
        self.assertEqual(build_periodic_alkane_setup(compact),(setup,data))
        self.assertNotIn('neigh_modify one',setup)

    def test_existing_small_native_continuation_script_remains_exact(self):
        cell=json.loads((ROOT/'evidence/geometry_v03/cell.json').read_text())['payload']
        assigned=assign_periodic_alkanes(cell)
        checkpoint=json.loads((ROOT/'evidence/restart_v06/preparation/protocol.json').read_text())['checkpoint']
        plan=plan_periodic_alkane_continuation(assigned,checkpoint,PeriodicAlkaneContinuationSpec(steps=40),
            checkpoint_directory='not-read',run_directory='not-created')
        old=(ROOT/'evidence/restart_v06/continuous_40/input.lammps').read_text()
        self.assertEqual(plan['input_script'],old)

    def test_large_thermo_preparation_and_restart_plans_share_explicit_storage(self):
        cell=pack_cell([analytic_ethane()]*65,box_lengths_angstrom=[80.,80.,80.],
            seed=84021,min_distance_angstrom=1.5).to_dict()
        assigned=assign_compact_periodic_alkanes(cell)
        self.assertEqual(len(assigned.to_dict()['atoms']),520)
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            with patch('polychemmd.engines.periodic_thermo.run_lammps') as run:
                thermo=plan_periodic_alkane_thermo_singlepoint(assigned,run_directory=root/'thermo')
                preparation=plan_periodic_alkane_protocol(assigned,PeriodicAlkaneProtocolSpec(),run_directory=root/'prep')
                resume=plan_periodic_alkane_continuation(assigned,synthetic_checkpoint(assigned),
                    PeriodicAlkaneContinuationSpec(),checkpoint_directory=root/'absent',run_directory=root/'resume')
                run.assert_not_called()
            self.assertFalse(list(root.iterdir()))
        command='neigh_modify one 2000 page 100000'
        for plan in (thermo,preparation,resume):
            script=plan['input_script']
            self.assertEqual(script.count(command),1)
            self.assertLess(script.index('neigh_modify every 1 delay 0 check yes'),script.index(command))
            self.assertIn('neighbor 2.0 bin',script)
        self.assertIn('read_restart source.restart',resume['input_script'])
        self.assertLess(resume['input_script'].index(command),resume['input_script'].index('run 0'))


if __name__=='__main__':
    unittest.main()

"""Native restart schema and failure fixtures; no engine is invoked here."""

from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from polychemmd.forcefields import assign_periodic_alkanes
from polychemmd.protocols.alkane import PeriodicAlkaneProtocolSpec, parse_periodic_final_dump
from polychemmd.protocols.restart import (
    PeriodicAlkaneContinuationSpec, _audit_data_identity, _digest_record, _state_equal,
    checkpoint_manifest, plan_periodic_alkane_continuation, run_periodic_alkane_continuation,
)

ROOT=Path(__file__).resolve().parents[1]


def assigned_fixture():
    record=json.loads((ROOT/'evidence/geometry_v03/cell.json').read_text())
    return assign_periodic_alkanes(record.get('payload',record))


def write_checkpoint_fixture(directory,assigned):
    # Snapshot values are from an old actual engine artifact, but the binary and
    # producer below deliberately represent a non-executable schema fixture.
    lines=(ROOT/'evidence/periodic_v05/preparation_pilot/final.dump').read_text().splitlines()
    lines[8]+=' vx vy vz'
    lines[9:]=[line+' 0.001 -0.002 0.003' for line in lines[9:]]
    (directory/'checkpoint.dump').write_text('\n'.join(lines)+'\n')
    (directory/'checkpoint.restart').write_bytes(b'NON_ENGINE_NATIVE_CHECKPOINT_SCHEMA_FIXTURE')
    return checkpoint_manifest(assigned.to_dict(),PeriodicAlkaneProtocolSpec(),
        {'engine_identity':{'fixture':'no_actual_engine'},'engine_version_banner':'FIXTURE_NOT_LAMMPS'},
        directory,step=200)


class ContinuationSpecTests(unittest.TestCase):
    def test_roundtrip_and_finite_resource_bounds(self):
        spec=PeriodicAlkaneContinuationSpec(steps=40,thermo_stride_steps=5)
        self.assertEqual(PeriodicAlkaneContinuationSpec.from_mapping(spec.to_dict()),spec)
        for overrides in ({'steps':True},{'steps':0},{'steps':10001},{'steps':1.2},
                {'timeout_seconds':False},{'timeout_seconds':float('nan')},
                {'timeout_seconds':121},{'thermo_stride_steps':0},{'thermo_stride_steps':True}):
            with self.subTest(overrides=overrides),self.assertRaises(ValueError):
                PeriodicAlkaneContinuationSpec(**overrides)
        with self.assertRaises(ValueError):
            PeriodicAlkaneContinuationSpec.from_mapping({'temperature_k':350})


class CheckpointTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.assigned=assigned_fixture()

    def test_plan_reads_no_artifacts_and_preserves_state_controls(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            checkpoint=write_checkpoint_fixture(root,self.assigned)
            with patch('polychemmd.protocols.restart.run_lammps') as run:
                plan=plan_periodic_alkane_continuation(self.assigned,checkpoint,
                    PeriodicAlkaneContinuationSpec(steps=40),checkpoint_directory=root/'absent',
                    run_directory=root/'new',executable='definitely_absent')
                run.assert_not_called()
            self.assertFalse((root/'new').exists())
            script=plan['input_script']
            self.assertIn('read_restart source.restart',script)
            self.assertIn('fix polychemmd_npt all npt temp 300 300 100.0',script)
            self.assertNotIn('reset_timestep',script)
            self.assertNotIn('velocity all create',script)
            self.assertNotIn('minimize ',script)
            self.assertEqual(plan['start_step'],200)
            self.assertEqual(plan['end_step_requested'],240)
            self.assertEqual(plan['duration_fs_requested'],10)
            self.assertEqual(plan['checkpoint_artifact_validation'],'deferred_until_execution')
            self.assertFalse(plan['production_ready'])

    def test_checkpoint_tampering_and_cross_assignment_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);checkpoint=write_checkpoint_fixture(root,self.assigned)
            changes=[('step',201),('units_style','metal'),('model_id','UNSUPPORTED'),
                     ('assignment_canonical_sha256','0'*64),('production_ready',True),
                     ('equilibrium_status','equilibrated'),('atom_count',1),
                     ('producer_engine_identity',{}),('engine_version_banner',None)]
            for field,value in changes:
                bad=deepcopy(checkpoint);bad[field]=value
                # Recomputed manifest SHA must not bypass semantic identity gates.
                bad['checkpoint_canonical_sha256']=_digest_record(bad)
                with self.subTest(field=field),self.assertRaises(ValueError):
                    plan_periodic_alkane_continuation(self.assigned,bad,PeriodicAlkaneContinuationSpec(),
                        checkpoint_directory=root,run_directory=root/'new')
            bad=deepcopy(checkpoint);bad['protocol_spec']['temperature_k']=350
            # Any unsigned edit without updating its digest also fails.
            with self.assertRaisesRegex(ValueError,'checksum'):
                plan_periodic_alkane_continuation(self.assigned,bad,PeriodicAlkaneContinuationSpec(),
                    checkpoint_directory=root,run_directory=root/'new')

    def test_binary_and_snapshot_hashes_fail_before_engine_or_output_creation(self):
        for artifact in ('checkpoint.restart','checkpoint.dump'):
            with self.subTest(artifact=artifact),tempfile.TemporaryDirectory() as temp:
                root=Path(temp);checkpoint=write_checkpoint_fixture(root,self.assigned)
                (root/artifact).write_bytes((root/artifact).read_bytes()+b'CHANGED')
                with patch('polychemmd.protocols.restart.run_lammps') as run,self.assertRaisesRegex(ValueError,'checksum'):
                    run_periodic_alkane_continuation(self.assigned,checkpoint,PeriodicAlkaneContinuationSpec(),
                        checkpoint_directory=root,run_directory=root/'new')
                    run.assert_not_called()
                self.assertFalse((root/'new').exists())

    def test_unsafe_symlink_checkpoint_fails_before_execution(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);checkpoint=write_checkpoint_fixture(root,self.assigned)
            (root/'saved.restart').write_bytes((root/'checkpoint.restart').read_bytes())
            (root/'checkpoint.restart').unlink();(root/'checkpoint.restart').symlink_to(root/'saved.restart')
            with self.assertRaisesRegex(ValueError,'unsafe'):
                run_periodic_alkane_continuation(self.assigned,checkpoint,PeriodicAlkaneContinuationSpec(),
                    checkpoint_directory=root,run_directory=root/'new')

    def test_malformed_box_values_fail_with_value_error(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);checkpoint=write_checkpoint_fixture(root,self.assigned)
            for field,value in (('lengths_angstrom',[10**1000,24,24]),
                    ('origin_angstrom',[True,0,0]),('periodic',[1,1,1]),
                    ('lengths_angstrom','invalid')):
                bad=deepcopy(checkpoint);bad['box'][field]=value
                bad['checkpoint_canonical_sha256']=_digest_record(bad)
                with self.subTest(field=field),self.assertRaises(ValueError):
                    plan_periodic_alkane_continuation(self.assigned,bad,PeriodicAlkaneContinuationSpec(),
                        checkpoint_directory=root,run_directory=root/'new')

    def test_artifact_size_budget_checked_before_reading(self):
        from polychemmd.protocols.restart import _artifact
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);artifact=root/'checkpoint.restart'
            with artifact.open('wb') as stream:stream.truncate(64*1024*1024+1)
            with patch.object(Path,'read_bytes',side_effect=AssertionError('must not read oversized input')):
                with self.assertRaisesRegex(ValueError,'64 MiB'):_artifact(root,'checkpoint.restart')

    def test_velocity_image_and_box_restore_checks_are_independent(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);write_checkpoint_fixture(root,self.assigned)
            original=parse_periodic_final_dump((root/'checkpoint.dump').read_text(),self.assigned.to_dict(),
                                              200,include_velocities=True)
            _state_equal(original,deepcopy(original))
            for field,value in (('velocities_angstrom_fs',[0,0,0]),('image_flags',[99,0,0]),
                                ('coordinates_angstrom',[0,0,0]),('forces_kj_mol_angstrom',[999,0,0])):
                changed=deepcopy(original);changed['atoms'][0][field]=value
                with self.subTest(field=field),self.assertRaises(ValueError):
                    _state_equal(original,changed)
            changed=deepcopy(original);changed['box']['origin_angstrom'][0]+=.1
            with self.assertRaisesRegex(ValueError,'box'):_state_equal(original,changed)

    def test_native_text_readback_obeys_pinned_coefficient_serialization(self):
        text=(ROOT/'evidence/periodic_v05/preparation_pilot/prepared.data').read_text()
        # OPLS source CCC coefficient is 58.350023900573618 kcal/mol; the native
        # binary retains it, while AngleHarmonic::write_data serializes %g=58.35.
        _audit_data_identity(text,self.assigned,PeriodicAlkaneProtocolSpec())
        for before,after in (('1 58.35 112.7','1 58.3501 112.7'),('1 12.011','1 12.111'),
                             ('1 1.3 -0.05 0.2 0','1 1.3 -0.05 0.3 0')):
            self.assertIn(before,text)
            with self.subTest(change=after),self.assertRaises(ValueError):
                _audit_data_identity(text.replace(before,after),self.assigned,PeriodicAlkaneProtocolSpec())

    def test_interaction_row_reordering_and_full_reversal_preserve_topology(self):
        text=(ROOT/'evidence/periodic_v05/preparation_pilot/prepared.data').read_text()
        lines=text.splitlines();current=None;counter=0
        for index,line in enumerate(lines):
            clean=line.split('#',1)[0].strip()
            if clean in {'Bonds','Angles','Dihedrals'}:current=clean;counter=0;continue
            if not clean:continue
            if current is not None:
                if not clean[0].isdigit():current=None;continue
                tokens=clean.split();counter+=1
                # Unique replacement row IDs and full reverse do not change the
                # type/endpoint multiset. Physical atom IDs remain exact.
                lines[index]=' '.join([str(10000-counter),tokens[1],*tokens[2:][::-1]])
        _audit_data_identity('\n'.join(lines),self.assigned,PeriodicAlkaneProtocolSpec())
        original=text.splitlines();index=original.index('Bonds')+2
        tokens=original[index].split();tokens[1]='99';original[index]=' '.join(tokens)
        with self.assertRaisesRegex(ValueError,'Bonds'):
            _audit_data_identity('\n'.join(original),self.assigned,PeriodicAlkaneProtocolSpec())


if __name__=='__main__':unittest.main()

"""Process-contract fixtures, not scientific LAMMPS integration results.

Temporary Python executables below deliberately emulate process outcomes only.
Actual LAMMPS execution must be recorded separately with its binary identity.
"""

from pathlib import Path
import hashlib
import json
import os
import sys
import tempfile
import unittest

from polychemmd.engines import (
    LammpsError, LammpsRunSpec, discover_lammps, parse_thermo_log,
    plan_lammps, run_lammps,
)


GOOD_LOG = '''LAMMPS process fixture; not an actual solver
Step Temp PotEng TotEng Press
0 1 -6 -4 0.5
20 0.98 -5.98 -4.01 0.49
Loop time of 0.001 on 1 procs for 20 steps with 108 atoms
POLYCHEMMD_RUN_COMPLETE
'''


@unittest.skipIf(os.name != "posix", "temporary executable process fixtures use POSIX shebangs")
class AdapterContractTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def fake(self, mode="success", log=GOOD_LOG):
        executable = self.root / f"fake solver {mode}.py"
        program = f'''#!{sys.executable}
import os, sys, time
from pathlib import Path
mode={mode!r}
if '-h' in sys.argv:
    if mode=='help-failure':
        print('loader error: deliberately emulated',file=sys.stderr)
        sys.exit(127)
    print('LAMMPS PROCESS-CONTRACT TEST FIXTURE, NOT SCIENTIFIC EVIDENCE')
    sys.exit(0)
if mode=='timeout':
    time.sleep(2)
if mode!='missing':
    path=Path(sys.argv[sys.argv.index('-log')+1])
    path.write_text('' if mode=='empty' else {log!r})
    if mode=='stale':
        os.utime(path,(1,1))
Path('smoke.data').write_text('fixture artifact only')
if mode=='generic-error-message':
    print('host diagnostic ERROR word is not a solver result',file=sys.stderr)
if mode=='nonzero':
    sys.exit(3)
'''
        executable.write_text(program)
        executable.chmod(0o755)
        return executable

    def spec(self, executable, **kwargs):
        return LammpsRunSpec("# trusted process fixture\nrun 20\n", self.root / "new run",
                             executable=executable, **kwargs)

    def test_success_checks_artifacts_hashes_and_parses_thermo(self):
        executable = self.fake()
        spec = self.spec(executable, expected_artifacts=("log.lammps", "smoke.data"),
                         input_files=(("seed.data", "input fixture"),))
        result = run_lammps(spec).to_dict()
        self.assertEqual(result['status'], 'completed')
        self.assertEqual(result['input_sha256'], hashlib.sha256(spec.input_script.encode()).hexdigest())
        self.assertEqual(result['executable_sha256'], hashlib.sha256(executable.read_bytes()).hexdigest())
        self.assertEqual(result['staged_input_sha256']['seed.data'], hashlib.sha256(b'input fixture').hexdigest())
        self.assertEqual(result['thermo_blocks'][0]['rows'][-1][0], 20)
        self.assertEqual(result['material_validation'], 'not_performed')
        self.assertEqual(result['equilibration_validation'], 'not_performed')
        self.assertEqual(result, json.loads((spec.run_directory/'run_record.json').read_text()))
        self.assertEqual(result['command'][0], str(executable.resolve()))

    def test_zero_exit_alone_is_insufficient(self):
        for mode in ('missing', 'empty', 'stale'):
            with self.subTest(mode=mode):
                spec = LammpsRunSpec('run 20', self.root/mode, executable=self.fake(mode))
                with self.assertRaises(LammpsError) as caught:
                    run_lammps(spec)
                record = caught.exception.record
                self.assertEqual(record['status'], 'failed')
                self.assertEqual(record['failure_phase'], 'artifact_validation')
                self.assertEqual(record['returncode'], 0)

    def test_nonzero_records_failure_despite_complete_artifacts(self):
        spec = self.spec(self.fake('nonzero'))
        with self.assertRaises(LammpsError) as caught:
            run_lammps(spec)
        self.assertEqual(caught.exception.record['returncode'], 3)
        self.assertEqual(caught.exception.record['failure_phase'], 'engine_run')

    def test_timeout_records_failure(self):
        spec = self.spec(self.fake('timeout'), timeout_seconds=.25)
        with self.assertRaises(LammpsError) as caught:
            run_lammps(spec)
        self.assertEqual(caught.exception.record['status'], 'failed')
        self.assertEqual(caught.exception.record['failure_phase'], 'engine_run')
        self.assertTrue((spec.run_directory/'run_record.json').is_file())

    def test_help_must_actually_execute(self):
        spec = self.spec(self.fake('help-failure'))
        with self.assertRaises(LammpsError) as caught:
            run_lammps(spec)
        self.assertFalse(caught.exception.record['simulation_performed'])
        self.assertEqual(caught.exception.record['failure_phase'], 'engine_identification')

    def test_existing_directory_cannot_reuse_old_result(self):
        spec = self.spec(self.fake())
        spec.run_directory.mkdir()
        previous = spec.run_directory/'log.lammps'
        previous.write_text(GOOD_LOG)
        with self.assertRaisesRegex(LammpsError, 'already exists'):
            run_lammps(spec)
        self.assertEqual(previous.read_text(), GOOD_LOG)

    def test_missing_engine_never_creates_directory(self):
        spec = self.spec(self.root/'missing-executable')
        with self.assertRaisesRegex(LammpsError, 'unavailable'):
            run_lammps(spec)
        self.assertFalse(spec.run_directory.exists())

    def test_marker_must_be_output_not_echoed_command(self):
        log = GOOD_LOG.replace('POLYCHEMMD_RUN_COMPLETE', 'print "POLYCHEMMD_RUN_COMPLETE"')
        with self.assertRaisesRegex(LammpsError, 'completion marker'):
            run_lammps(self.spec(self.fake(log=log)))

    def test_marker_must_follow_last_loop(self):
        log = 'POLYCHEMMD_RUN_COMPLETE\n' + GOOD_LOG.replace('POLYCHEMMD_RUN_COMPLETE', '')
        with self.assertRaisesRegex(LammpsError, 'completion marker'):
            run_lammps(self.spec(self.fake(log=log)))

    def test_nonfinite_thermo_is_not_completed_evidence(self):
        log = GOOD_LOG.replace('20 0.98', '20 NaN')
        with self.assertRaisesRegex(LammpsError, 'nonfinite'):
            run_lammps(self.spec(self.fake(log=log)))

    def test_generic_stderr_diagnostic_does_not_replace_process_gates(self):
        result = run_lammps(self.spec(self.fake('generic-error-message'))).to_dict()
        self.assertEqual(result['status'], 'completed')

    def test_explicit_lammps_error_is_rejected(self):
        with self.assertRaisesRegex(LammpsError, 'reported an ERROR'):
            run_lammps(self.spec(self.fake(log=GOOD_LOG+'ERROR: deliberate process fixture\n')))

    def test_units_and_version_do_not_mislabel_diagnostics(self):
        executable = self.fake()
        # The fake help text is deliberately not a real solver version.
        spec = LammpsRunSpec('units lj\nrun 20', self.root/'reduced run', executable=executable)
        result = run_lammps(spec).to_dict()
        self.assertEqual(result['input_units_style'], 'lj')
        self.assertIn('not_performed', result['thermo_unit_conversion'])
        self.assertIsNone(result['engine_version_banner'])

    def test_plan_runs_nothing_and_creates_no_result(self):
        spec = self.spec(self.root/'not-installed')
        result = plan_lammps(spec)
        self.assertEqual(result['status'], 'planned')
        self.assertFalse(result['simulation_performed'])
        self.assertNotIn('thermo_blocks', result)
        self.assertFalse(spec.run_directory.exists())


class ThermoAndSpecTests(unittest.TestCase):
    def test_clear_resets_units_to_default_lj(self):
        with tempfile.TemporaryDirectory() as directory:
            spec = LammpsRunSpec('units real\nclear\nrun 0', directory)
            result = plan_lammps(spec)
            self.assertEqual(result['input_units_style'], 'lj')
            self.assertEqual(result['units_style_confidence'], 'recipe_literal_not_engine_verified')
            with self.assertRaisesRegex(ValueError, 'conflicts'):
                LammpsRunSpec('units real\nclear\nrun 0', directory, units_style='real')

    def test_units_after_clear_and_earlier_runs_are_distinguished(self):
        with tempfile.TemporaryDirectory() as directory:
            self.assertEqual(plan_lammps(LammpsRunSpec('clear\nunits real\nrun 0', directory))['input_units_style'], 'real')
            # A reset after a completed run does not retroactively change its thermo units.
            self.assertEqual(plan_lammps(LammpsRunSpec('units real\nrun 0\nclear', directory))['input_units_style'], 'real')
            multiple = 'units real\nrun 0\nclear\nrun 0'
            result = plan_lammps(LammpsRunSpec(multiple, directory))
            self.assertIsNone(result['input_units_style'])
            self.assertEqual(result['recipe_run_units_styles'], ['real', 'lj'])
            with self.assertRaisesRegex(ValueError, 'conflicts'):
                LammpsRunSpec(multiple, directory, units_style='real')

    def test_dynamic_recipe_units_stay_unresolved_despite_declaration(self):
        with tempfile.TemporaryDirectory() as directory:
            for command in ('include setup.in', 'jump SELF loop', 'variable u string real',
                            'if "1" then "units lj"', 'python callback invoke',
                            'units ${u}', 'run 10 every 5 "clear"', 'units real &\n'):
                with self.subTest(command=command):
                    spec = LammpsRunSpec('units real\n'+command+'\nrun 0', directory, units_style='real')
                    result = plan_lammps(spec)
                    self.assertIsNone(result['input_units_style'])
                    self.assertEqual(result['declared_units_style'], 'real')
                    self.assertEqual(result['units_style_confidence'], 'unresolved')

    def test_multiple_completed_blocks_and_warning(self):
        log = GOOD_LOG.replace('POLYCHEMMD_RUN_COMPLETE', '') + GOOD_LOG.replace('0 1 -6', 'WARNING: test line\n0 1 -6')
        blocks = parse_thermo_log(log)
        self.assertEqual(len(blocks), 2)
        self.assertEqual(blocks[0].columns[0], 'Step')

    def test_truncated_malformed_and_missing_required_column(self):
        variants = [GOOD_LOG.split('Loop time')[0],
                    GOOD_LOG.replace('20 0.98 -5.98 -4.01 0.49', '20 0.98'),
                    GOOD_LOG.replace('Temp', 'Temperature'),
                    GOOD_LOG.replace('0 1 -6 -4 0.5', 'not a numeric data line'),
                    GOOD_LOG.replace('0 1 -6 -4 0.5', '30 1 -6 -4 0.5')]
        for log in variants:
            with self.subTest(log=log), self.assertRaises(ValueError):
                parse_thermo_log(log)

    def test_empty_or_no_rows_is_invalid(self):
        for log in ('', 'Step Temp\nLoop time of 0 seconds\n'):
            with self.assertRaises(ValueError):
                parse_thermo_log(log)

    def test_path_traversal_collisions_and_invalid_timeout(self):
        with tempfile.TemporaryDirectory() as directory:
            for value in (True, 0, -1, float('nan'), 10**400, 3601):
                with self.subTest(timeout=value), self.assertRaises(ValueError):
                    LammpsRunSpec('run 0', directory, timeout_seconds=value)
            for name in ('../old.data', '/tmp/old.data', 'a/../old.data', 'a\\old.data'):
                with self.subTest(name=name), self.assertRaises(ValueError):
                    LammpsRunSpec('run 0', directory, expected_artifacts=(name,))
            with self.assertRaises(ValueError):
                LammpsRunSpec('run 0', directory, input_files=(("log.lammps", "stale"),))
            with self.assertRaises(ValueError):
                LammpsRunSpec('run 0', directory, expected_artifacts='log.lammps')
            with self.assertRaisesRegex(ValueError, 'conflicts'):
                LammpsRunSpec('units lj\nrun 0', directory, units_style='real')


if __name__ == '__main__':
    unittest.main()

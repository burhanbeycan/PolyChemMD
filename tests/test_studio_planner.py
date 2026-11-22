"""Check the local Ollama wire contract, whitelist, output size and proposal budget."""

import http.client
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
from threading import Thread
import unittest
from unittest.mock import patch

from polychemmd.studio.planner import PAIRS, suggest


class TestPlanner(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        class Stub(BaseHTTPRequestHandler):
            def do_POST(self):
                req = json.loads(self.rfile.read(int(self.headers["Content-Length"])))
                self.server.captured.append(
                    (self.path, self.headers["Content-Type"], req)
                )
                body = self.server.reply
                self.send_response(self.server.code)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)

            def log_message(self, *args):
                pass

        cls.server = ThreadingHTTPServer(("127.0.0.1", 0), Stub)
        cls.server.captured = []
        cls.thread = Thread(target=cls.server.serve_forever, daemon=True)
        cls.thread.start()

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.thread.join()
        cls.server.server_close()

    def setUp(self):
        # Test over real loopback HTTP on an ephemeral port without touching a
        # researcher's Ollama daemon. The requested *production* URL is asserted.
        real_connection = http.client.HTTPConnection

        def fixture_connection(host, port, timeout):
            self.assertEqual((host, port), ("127.0.0.1", 11434))
            return real_connection(host, self.server.server_port, timeout=timeout)

        transport = patch(
            "polychemmd.studio.planner.http.client.HTTPConnection",
            side_effect=fixture_connection,
        )
        self.client = transport.start()
        self.addCleanup(transport.stop)
        self.server.code = 200
        self.server.captured.clear()

    def reply(self, envelope):
        self.server.reply = json.dumps(envelope).encode()
        return suggest("ollama", "smollm2:135m", "termination")

    def test_http_schema_and_job_budget(self):
        result = self.reply(
            {"done": True, "response": json.dumps({"pair": "termination"})}
        )
        path, content_type, body = self.server.captured[0]
        self.assertEqual((path, content_type), ("/api/generate", "application/json"))
        self.assertFalse(body["stream"])
        self.assertEqual(body["model"], "smollm2:135m")
        self.assertEqual(body["options"]["num_predict"], 64)
        schema = body["format"]
        self.assertEqual(set(schema["properties"]["pair"]["enum"]), set(PAIRS))
        self.assertFalse(schema["additionalProperties"])
        self.assertEqual(len(result["jobs"]), 2)
        self.assertFalse(result["submitted"])
        self.assertEqual([j["config"]["seed"] for j in result["jobs"]], [17, 17])

    def test_invalid_local_responses(self):
        for env in (
            None,
            [],
            {"done": False},
            {"done": True},
            {"done": True, "response": {}},
            {"done": True, "response": "not JSON"},
            {"done": True, "response": '{"pair":[]}'},
            {"done": True, "response": '{"pair":"termination","shell":"whoami"}'},
            {"done": True, "response": '{"pair":"other"}'},
        ):
            with self.subTest(envelope=env), self.assertRaises(ValueError):
                self.reply(env)
        self.server.code = 404
        with self.assertRaisesRegex(ValueError, "HTTP 404"):
            self.reply({})
        self.server.code = 200
        with self.assertRaises(ValueError):
            self.reply({"done": True, "response": "x" * 65537})

    def test_rules_independent_and_immutable(self):
        r = suggest("rules", objective="hydrolysis")
        r["jobs"][0]["config"]["seed"] = 123
        self.assertEqual(
            suggest("rules", objective="hydrolysis")["jobs"][0]["config"]["seed"], 17
        )
        self.client.assert_not_called()

    def test_bad_inputs_make_no_model_request(self):
        for mode, model, objective in (
            ("ollama", "x y", "model"),
            ("no", "ok", "model"),
            ("ollama", "ok", "typo"),
            ("rules", None, ["termination"]),
        ):
            with (
                self.subTest(mode=mode, model=model, objective=objective),
                self.assertRaises(ValueError),
            ):
                suggest(mode, model, objective)
        self.client.assert_not_called()

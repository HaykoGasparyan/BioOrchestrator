"""Safety and scientific labeling checks for the workbench's data presentation."""

import csv
import importlib.util
import io
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from biocentric.ui import artifact_bytes, candidate_csv, is_fixture, pipeline_summary, safe_data
from biocentric.viewer import resolve_scene, script_json, viewer_html


class WorkbenchPresentationTests(unittest.TestCase):
    def test_fixture_labels_follow_persisted_mode_and_not_project_name(self):
        self.assertTrue(is_fixture({"project": {"spec": {"mode": "offline_fixture"}}}))
        self.assertTrue(is_fixture({"data": {"fixture": "public-ui-example"}}))
        self.assertFalse(is_fixture({"project": {"spec": {"name": "Real example"}}}))

    def test_blocked_and_skipped_stages_are_not_reported_as_complete(self):
        summary = pipeline_summary({"steps": [
            {"id": "dossier", "status": "completed"},
            {"id": "benchmark", "status": "blocked"},
            {"id": "screen", "status": "skipped"},
        ], "jobs": [{"status": "completed"}]})
        self.assertEqual(summary["complete"], 1)
        self.assertEqual(summary["current"]["id"], "benchmark")
        self.assertFalse(summary["active"])

    def test_only_contained_regular_artifacts_can_be_downloaded(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "project"
            root.mkdir()
            owned = root / "result.txt"
            owned.write_text("scientific record")
            outside = Path(temporary) / "secret.txt"
            outside.write_text("private")
            (root / "link.txt").symlink_to(outside)
            self.assertEqual(artifact_bytes({"path": str(owned)}, root), b"scientific record")
            self.assertIsNone(artifact_bytes({"path": str(owned), "sha256": "wrong-hash"}, root))
            self.assertIsNone(artifact_bytes({"path": "../secret.txt"}, root))
            self.assertIsNone(artifact_bytes({"path": str(root / "link.txt")}, root))

    def test_csv_preserves_numeric_scores_but_neutralizes_formula_text(self):
        output = candidate_csv([{"compound_id": '=HYPERLINK("https://example.invalid")', "score": -7.3,
                                 "flags": ["unknown ADMET"], "affinity": None}])
        row = next(csv.DictReader(io.StringIO(output)))
        self.assertTrue(row["compound_id"].startswith("'="))
        self.assertEqual(row["score"], "-7.3")
        self.assertEqual(row["affinity"], "")
        self.assertIn("unknown ADMET", row["flags"])

    def test_call_inspection_redacts_nested_credentials_but_keeps_usage(self):
        actual = safe_data({"request": {"Authorization": "Bearer private", "api_key": "private"},
                            "usage": {"input_tokens": 12}, "observations": [{"access_token": "private"}]})
        self.assertEqual(actual["request"]["Authorization"], "[redacted]")
        self.assertEqual(actual["observations"][0]["access_token"], "[redacted]")
        self.assertEqual(actual["usage"]["input_tokens"], 12)


class LocalGeometryViewerTests(unittest.TestCase):
    def test_untrusted_structure_text_cannot_close_script_or_fetch_external_data(self):
        malicious = "</script><script>fetch('https://example.invalid/private')</script><img src=x>\u2028"
        encoded = script_json({"text": malicious})
        self.assertNotIn("<", encoded)
        self.assertEqual(json.loads(encoded)["text"], malicious)
        document = viewer_html({"models": [{"text": malicious}], "coordinate_frame": malicious})
        self.assertNotIn(malicious, document)
        self.assertIn("connect-src 'none'", document)
        self.assertIn("worker-src 'none'", document)
        self.assertNotIn('<script src=', document)

    def test_registered_scene_requires_frame_and_verified_candidate_linkage(self):
        from biocentric.service import Service

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            service = Service(root)
            project_id = service.create_example()
            snapshot = service.snapshot(project_id)
            receptor = next(a for a in snapshot["artifacts"] if a["kind"] == "structure")
            docking = next(a for a in snapshot["artifacts"] if a["kind"] == "receptor")
            reference = next(a for a in snapshot["artifacts"] if a["kind"] == "reference")
            with self.assertRaisesRegex(ValueError, "Coordinate frame"):
                resolve_scene(snapshot, root, receptor["id"])
            with self.assertRaisesRegex(ValueError, "Confirm"):
                resolve_scene(snapshot, root, receptor["id"], reference["id"], declared_frame="original input")
            scene = resolve_scene(snapshot, root, receptor["id"], declared_frame="original input")
            self.assertEqual(scene["models"][0]["sha256"], receptor["sha256"])
            self.assertEqual(scene["frame_source"], "user-supplied display annotation")
            snapshot["data"]["protocol"] = {"input_hashes": {"structure": receptor["sha256"], "receptor": docking["sha256"]},
                                                "box": {"coordinate_frame": "untransformed original receptor"}}
            pose = {"artifact_ids": [reference["id"]], "pose_sha256": reference["sha256"], "receptor_sha256": docking["sha256"],
                    "pose_path": "/outside/never/read.sdf"}
            scene = resolve_scene(snapshot, root, receptor["id"], reference["id"], pose)
            self.assertEqual(len(scene["models"]), 2)
            self.assertEqual(scene["frame_source"], "recorded protocol")
            with self.assertRaisesRegex(ValueError, "artifact IDs"):
                resolve_scene(snapshot, root, receptor["id"], reference["id"], {**pose, "artifact_ids": []})
            with self.assertRaisesRegex(ValueError, "pose hash"):
                resolve_scene(snapshot, root, receptor["id"], reference["id"], {**pose, "pose_sha256": "wrong"})
            with self.assertRaisesRegex(ValueError, "approved protocol"):
                resolve_scene(snapshot, root, receptor["id"], reference["id"], {**pose, "receptor_sha256": "wrong"})


@unittest.skipUnless(importlib.util.find_spec("streamlit"), "Streamlit is not installed")
class WorkbenchInteractionTests(unittest.TestCase):
    def setUp(self):
        self.classic_mode = patch.dict(os.environ, {"BIOCENTRIC_UI_MODE": "classic"})
        self.classic_mode.start()
        self.addCleanup(self.classic_mode.stop)

    def test_intake_persists_typed_inputs_without_a_model_or_scientific_call(self):
        import streamlit as st
        from streamlit.testing.v1 import AppTest
        from biocentric.service import Service

        with tempfile.TemporaryDirectory() as temporary, patch.dict(os.environ, {"BIOCENTRIC_HOME": temporary}):
            st.cache_resource.clear()
            app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / "app.py"), default_timeout=20).run()
            next(item for item in app.text_input if item.label == "Project name").set_value("Private research")
            next(item for item in app.text_area if item.label == "Protein, accession, or sequence").set_value("P00519")
            cpu_input = next(item for item in app.number_input if item.label == "CPU cores")
            self.assertEqual(cpu_input.max, 16)
            cpu_input.set_value(16)
            next(button for button in app.button if button.label == "Create project").click().run()
            self.assertFalse(app.exception)
            self.assertFalse(app.error)
            project = Service(temporary).list_projects()[0]
            self.assertEqual(project["spec"]["target"], "P00519")
            self.assertEqual(project["spec"]["cpu"], 16)
            self.assertTrue(project["spec"]["local_only"])
            self.assertFalse(project["spec"]["allow_model"])
            self.assertIsNone(project["spec"]["intended_effect"])
            self.assertEqual(Service(temporary).snapshot(project["id"])["jobs"], [])
            st.cache_resource.clear()

    def test_executable_example_registers_real_inputs_without_starting_a_job(self):
        import streamlit as st
        from streamlit.testing.v1 import AppTest
        from biocentric.service import Service

        with tempfile.TemporaryDirectory() as temporary, patch.dict(os.environ, {"BIOCENTRIC_HOME": temporary}):
            st.cache_resource.clear()
            app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / "app.py"), default_timeout=20).run()
            next(button for button in app.button if button.label == "Load executable smoke example").click().run()
            self.assertFalse(app.exception)
            self.assertFalse(any("OFFLINE DEMONSTRATION" in item.value for item in app.markdown))
            self.assertFalse(next(button for button in app.button if button.label == "▶ Run pipeline").disabled)
            service = Service(temporary)
            snapshot = service.snapshot(service.list_projects()[0]["id"])
            self.assertEqual(snapshot["project"]["spec"]["mode"], "live")
            self.assertEqual(snapshot["project"]["spec"]["organism"], "Mus musculus")
            self.assertEqual(snapshot["jobs"], [])
            self.assertIn("receptor", {artifact["kind"] for artifact in snapshot["artifacts"]})
            app.checkbox(key="file_geometry_enabled").check().run()
            self.assertFalse(app.exception)
            self.assertTrue(any("Coordinate frame" in error.value for error in app.error))
            app.text_input(key="file_geometry_frame").set_value("Original uploaded PDB coordinates").run()
            self.assertFalse(app.exception)
            self.assertFalse(app.error)
            self.assertEqual(len(app.get("iframe")), 1)
            st.cache_resource.clear()

    def test_offline_demo_discloses_fixture_and_exposes_all_ten_stages(self):
        import streamlit as st
        from streamlit.testing.v1 import AppTest

        with tempfile.TemporaryDirectory() as temporary, patch.dict(os.environ, {"BIOCENTRIC_HOME": temporary}):
            st.cache_resource.clear()
            app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / "app.py"), default_timeout=20).run()
            self.assertFalse(app.exception)
            next(button for button in app.button if button.label == "Load offline demonstration").click().run()
            self.assertFalse(app.exception)
            self.assertTrue(any("OFFLINE DEMONSTRATION" in item.value for item in app.markdown))
            self.assertEqual(len([button for button in app.button if button.label == "Inspect stage"]), 10)
            self.assertTrue(next(button for button in app.button if button.label == "▶ Run pipeline").disabled)
            self.assertEqual([tab.label for tab in app.tabs], ["Workflow", "Model activity", "Evidence", "Plan & Runs", "Candidates", "Files & handoff"])
            st.cache_resource.clear()

    def test_model_status_and_opt_in_example_do_not_make_api_calls(self):
        import streamlit as st
        from streamlit.testing.v1 import AppTest
        from biocentric.service import Service

        with tempfile.TemporaryDirectory() as temporary, patch.dict(os.environ, {"BIOCENTRIC_HOME": temporary, "OPENAI_API_KEY": "synthetic-test-key"}), patch("biocentric.agent.advise") as advise:
            service = Service(temporary)
            local = service.create_example()
            st.cache_resource.clear()
            app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / "app.py"), default_timeout=20).run()
            self.assertFalse(app.exception)
            self.assertTrue(any("LLM disabled" in item.value and "0 API requests" in item.value for item in app.warning))
            next(button for button in app.button if button.label == "Load GPT-enabled public example").click().run()
            self.assertFalse(app.exception)
            project = service.list_projects()[0]
            self.assertNotEqual(project['id'], local)
            self.assertTrue(project['spec']['allow_model'])
            self.assertFalse(project['spec']['local_only'])
            self.assertEqual(service.snapshot(project['id'])['jobs'], [])
            self.assertFalse(service.snapshot(local)['project']['spec']['allow_model'])
            self.assertTrue(any("LLM waiting" in item.value for item in app.info))
            advise.assert_not_called()
            st.cache_resource.clear()

    def test_saved_model_advice_is_readable_without_expanding_json(self):
        import streamlit as st
        from streamlit.testing.v1 import AppTest
        from biocentric.service import Service

        with tempfile.TemporaryDirectory() as temporary, patch.dict(os.environ, {"BIOCENTRIC_HOME": temporary}):
            service = Service(temporary)
            pid = service.create_example(allow_model=True)
            advice = {'summary': 'A matching receptor is available; review the observed binding site.',
                      'proposed_strategy': 'structure_based', 'claims': [
                          {'text': 'The supplied structure is the 1IEP tutorial input.',
                           'category': 'retrieved_observation', 'evidence_ids': ['example-1iep']}],
                      'limitations': ['Predictive performance has not been validated.'], 'unresolved_questions': []}
            service.store.update(pid, data={'model_advice': advice})
            st.cache_resource.clear()
            app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / "app.py"), default_timeout=20).run()
            self.assertFalse(app.exception)
            visible = '\n'.join(item.value for item in app.markdown)
            self.assertIn(advice['summary'], visible)
            self.assertIn('GPT explanation / advisory proposal', visible)
            self.assertIn('Decision summary', visible)
            st.cache_resource.clear()

    def test_local_identity_review_is_persisted_from_the_interface(self):
        import streamlit as st
        from streamlit.testing.v1 import AppTest
        from biocentric.service import Service

        with tempfile.TemporaryDirectory() as temporary, patch.dict(os.environ, {"BIOCENTRIC_HOME": temporary}):
            service = Service(temporary)
            project_id = service.create_project({"name": "Local target", "target": "P00519", "organism": "Homo sapiens"})
            service.store.update(project_id, state="blocked", data={"gate": {"kind": "identity"}})
            st.cache_resource.clear()
            app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / "app.py"), default_timeout=20).run()
            self.assertFalse(app.exception)
            identity = {"accession": "P00519", "name": "ABL1", "organism": "Homo sapiens"}
            next(item for item in app.text_area if item.label == "TargetIdentity JSON").set_value(json.dumps(identity))
            next(item for item in app.text_input if item.label == "Identity provenance and rationale").set_value("Reviewed laboratory construct record")
            next(button for button in app.button if button.label == "Record local identity").click().run()
            self.assertFalse(app.exception)
            decision = service.snapshot(project_id)["decisions"][-1]
            self.assertEqual(decision["kind"], "identity")
            self.assertEqual(decision["value"]["accession"], "P00519")
            self.assertEqual(decision["reason"], "Reviewed laboratory construct record")
            st.cache_resource.clear()


if __name__ == "__main__":
    unittest.main()

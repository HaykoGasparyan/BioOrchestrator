"""Offline 3Dmol presentation of registered, integrity-checked molecular geometry."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

VIEWER_VERSION = "2.5.5"
VIEWER_SHA256 = "f7cc78921ae72e7623e89cdd111434f58c2efddd2ffda1cd212644b406fb8016"


def script_json(value: Any) -> str:
    """JSON in an HTML script element must not contain a literal closing tag."""
    return (json.dumps(value, ensure_ascii=True, allow_nan=False)
            .replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026"))


def geometry_artifacts(snapshot: dict) -> tuple[list[dict], list[dict]]:
    """Only the selected project's registered structure and ligand records qualify."""
    project_id = snapshot.get("project", {}).get("id")
    artifacts = [a for a in snapshot.get("artifacts", []) if a.get("project_id") == project_id
                 and a.get("id") and a.get("sha256")]
    receptors = [a for a in artifacts if a.get("kind") == "structure"
                 and Path(a.get("name", "")).suffix.lower() in {".pdb", ".cif", ".mmcif"}]
    ligands = [a for a in artifacts if a.get("kind") in {"reference", "scientific_output"}
               and Path(a.get("name", "")).suffix.lower() == ".sdf"]
    return receptors, ligands


def resolve_scene(snapshot: dict, root: Path, receptor_id: str, ligand_id: str | None = None,
                  pose_record: dict | None = None, declared_frame: str = "",
                  manual_overlay_confirmed: bool = False) -> dict:
    """Never dereference a pose_path or accept a structure outside this project."""
    from .ui import artifact_bytes

    receptors, ligands = geometry_artifacts(snapshot)
    receptor = next((a for a in receptors if a["id"] == receptor_id), None)
    ligand = next((a for a in ligands if a["id"] == ligand_id), None) if ligand_id else None
    if receptor is None or (ligand_id and ligand is None):
        raise ValueError("Choose registered structure and SDF artifact IDs from this project.")
    protocol = snapshot.get("data", {}).get("protocol") or {}
    hashes = protocol.get("input_hashes") or {}
    recorded_frame = (protocol.get("box") or {}).get("coordinate_frame")
    frame = recorded_frame or declared_frame.strip()
    if not frame:
        raise ValueError("Coordinate frame is not recorded. Supply the frame/source before displaying geometry.")
    if len(frame) > 500:
        raise ValueError("Coordinate-frame description must be at most 500 characters.")
    if pose_record is not None:
        if not ligand or ligand["id"] not in pose_record.get("artifact_ids", []):
            raise ValueError("The selected pose must be referenced by this candidate's registered artifact IDs.")
        if ligand["sha256"] != pose_record.get("pose_sha256"):
            raise ValueError("Candidate pose hash does not match the registered SDF.")
        if not hashes.get("structure") or receptor["sha256"] != hashes["structure"]:
            raise ValueError("The receptor structure does not match this pose's approved protocol.")
        if not hashes.get("receptor") or pose_record.get("receptor_sha256") != hashes["receptor"]:
            raise ValueError("The pose's docking receptor does not match the approved protocol.")
    elif ligand and not manual_overlay_confirmed:
        raise ValueError("Confirm that the selected structure and SDF share the same coordinate frame before overlaying them.")

    models = []
    for artifact, role in [(receptor, "receptor"), (ligand, "ligand")]:
        if artifact is None:
            continue
        # Restrict to the selected project's immutable artifact directory, then verify SHA.
        path = Path(artifact.get("path", "")).resolve()
        project_root = (root / "projects" / snapshot["project"]["id"]).resolve()
        if not path.is_relative_to(project_root):
            raise ValueError("Registered geometry must remain inside its project artifact directory.")
        content = artifact_bytes(artifact, root)
        if content is None:
            raise ValueError("Geometry artifact is missing, changed, or outside the allowed artifact store.")
        if len(content) > 10 * 1024 * 1024:
            raise ValueError("3D preview is limited to 10 MB per geometry artifact.")
        text = content.decode("utf-8-sig", errors="strict")
        suffix = Path(artifact["name"]).suffix.lower()
        format_name = "sdf" if role == "ligand" else "pdb" if suffix == ".pdb" else "cif"
        if format_name == "pdb" and not any(line.startswith(("ATOM  ", "HETATM")) for line in text.splitlines()):
            raise ValueError("PDB artifact contains no recognizable atom records.")
        if format_name == "cif" and "_atom_site." not in text:
            raise ValueError("mmCIF artifact contains no atom-site records.")
        if format_name == "sdf" and ("M  END" not in text or "$$$$" not in text):
            raise ValueError("SDF artifact is incomplete or unsupported for this viewer.")
        models.append({"role": role, "format": format_name, "text": text,
                       "name": artifact["name"], "artifact_id": artifact["id"], "sha256": artifact["sha256"]})
    return {"models": models, "coordinate_frame": frame,
            "frame_source": "recorded protocol" if recorded_frame else "user-supplied display annotation",
            "interpretation": "Geometry inspection only; a plausible pose does not confirm binding."}


def viewer_html(scene: dict) -> str:
    library = (Path(__file__).parent / "assets" / f"3Dmol-{VIEWER_VERSION}-min.js").read_bytes()
    if hashlib.sha256(library).hexdigest() != VIEWER_SHA256:
        raise ValueError("The bundled 3D viewer failed its version integrity check.")
    # A fixed vendored library runs inside a frame whose CSP disables all network access.
    source = library.decode("utf-8").replace("</script", "<\\/script")
    return """<!doctype html><html><head><meta charset="utf-8">
<meta http-equiv="Content-Security-Policy" content="default-src 'none'; script-src 'unsafe-inline'; style-src 'unsafe-inline'; img-src data: blob:; connect-src 'none'; worker-src 'none'; object-src 'none'; base-uri 'none'; form-action 'none'">
<style>body{margin:0;font:13px system-ui;color:#294552;background:#f5f9fa}.toolbar{display:flex;align-items:center;gap:8px;padding:10px;border-bottom:1px solid #dbe7e9}button{padding:7px 12px;background:white;border:1px solid #b8cfd2;border-radius:6px;color:#176773;cursor:pointer}#viewer{height:430px;position:relative}#status{padding:10px 12px;font-size:12px;line-height:1.5}#error{color:#a4412f;padding:15px;display:none}.hint{margin-left:auto;color:#587780}</style></head><body>
<div class="toolbar"><button id="all">Fit structure</button><button id="site">Focus ligand</button><button id="spin">Rotate</button><span class="hint">Drag to rotate · wheel / pinch to zoom</span></div>
<div id="viewer"></div><div id="error" role="alert"></div><div id="status"></div>
<script>""" + source + """</script><script type="application/json" id="scene-data">""" + script_json(scene) + """</script>
<script>(function(){'use strict';try{
const scene=JSON.parse(document.getElementById('scene-data').textContent);
const viewer=$3Dmol.createViewer(document.getElementById('viewer'),{backgroundColor:'#f5f9fa',antialias:true});
let ligand=null, rotating=false;
for(const model of scene.models){const added=viewer.addModel(model.text,model.format,{keepH:true});
if(added.selectedAtoms({}).length===0)throw new Error('No atoms could be displayed from '+model.name);
if(model.role==='receptor'){added.setStyle({},{cartoon:{color:'#69aab0',opacity:0.8}});added.setStyle({hetflag:true},{stick:{radius:0.16,colorscheme:'grayCarbon'}});}
else{ligand=added;added.setStyle({},{stick:{radius:0.22,colorscheme:model.role==='reference'?'cyanCarbon':'orangeCarbon'}});}}
viewer.zoomTo();viewer.render();
document.getElementById('all').onclick=()=>{viewer.zoomTo();viewer.render();};
document.getElementById('site').disabled=!ligand;
document.getElementById('site').onclick=()=>{viewer.zoomTo({model:ligand});viewer.zoom(0.7);viewer.render();};
document.getElementById('spin').onclick=()=>{rotating=!rotating;viewer.spin(rotating?'y':false);document.getElementById('spin').textContent=rotating?'Stop rotation':'Rotate';};
document.getElementById('status').textContent=scene.interpretation+' Frame: '+scene.coordinate_frame+' ('+scene.frame_source+'). 3Dmol.js 2.5.5 · local assets only.';
window.addEventListener('resize',()=>{viewer.resize();viewer.render();});
}catch(error){document.getElementById('error').style.display='block';document.getElementById('error').textContent='3D preview unavailable: '+error.message;}})();</script></body></html>"""

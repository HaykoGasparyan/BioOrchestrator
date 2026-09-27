"""Download a pinned official Vina binary for this platform into .tools/.

Run explicitly: python examples/install_vina.py
This installer does not overwrite an existing executable.
"""
import hashlib
import json
import os
import platform
import subprocess
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VERSION = "1.2.7"
# Verified Apple Silicon asset; remaining platforms retain computed hashes in manifest.
MAC_ARM64_SHA256 = json.loads((ROOT / "examples/1iep/sources.json").read_text())["vina_binary"]["sha256"]


def main():
    system = {"Darwin": "mac", "Linux": "linux", "Windows": "win"}.get(platform.system())
    machine = {"arm64": "aarch64", "aarch64": "aarch64", "x86_64": "x86_64", "AMD64": "x86_64"}.get(platform.machine())
    if not system or not machine:
        raise SystemExit("Unsupported platform; install official Vina separately and set BIOCENTRIC_VINA.")
    name = f"vina_{VERSION}_win.exe" if system == "win" else f"vina_{VERSION}_{system}_{machine}"
    url = f"https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v{VERSION}/{name}"
    destination = ROOT / ".tools" / ("vina.exe" if system == "win" else "vina")
    if destination.exists():
        raise SystemExit(f"Already present: {destination}; left unchanged.")
    destination.parent.mkdir(exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "BioCentricGPT-install/0.1"})
    with urllib.request.urlopen(request, timeout=60) as response:
        content = response.read(20 * 1024 * 1024)
    checksum = hashlib.sha256(content).hexdigest()
    if system == "mac" and machine == "aarch64" and checksum != MAC_ARM64_SHA256:
        raise SystemExit("Official release bytes differ from verified pinned SHA-256; download rejected.")
    with destination.open("xb") as handle:
        handle.write(content)
    destination.chmod(0o755)
    version = subprocess.run([str(destination), "--version"], capture_output=True, text=True, timeout=10, check=True).stdout.strip()
    (destination.parent / "vina-install.json").write_text(json.dumps({"url": url, "sha256": checksum, "version": version, "platform": platform.platform()}, indent=2))
    print(f"Installed {version}: {destination}")


if __name__ == "__main__":
    main()

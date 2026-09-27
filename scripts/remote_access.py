"""Password-protected, temporary HTTPS access to the local Bio Orchestrator app."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import secrets
import subprocess
import time
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
PRIVATE = ROOT / ".biocentric/remote-access"
GATEWAY = "bio-coordinator-gateway"
TUNNEL = "bio-coordinator-tunnel"


def command(*args: str, check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(args, capture_output=True, text=True, check=check)


def active(name: str) -> bool:
    return command("systemctl", "--user", "is-active", name, check=False).stdout.strip() == "active"


def launch(name: str, *args: str) -> None:
    if active(name):
        return
    command("systemd-run", "--user", "--collect", f"--unit={name}",
            "--property=Restart=on-failure", "--property=RestartSec=5",
            f"--working-directory={ROOT}", *args)


def status() -> str | None:
    if not active(TUNNEL):
        print("Remote access is stopped.")
        return None
    # Restrict logs to this service invocation; an old hostname may no longer work.
    invocation = command("systemctl", "--user", "show", TUNNEL,
                         "--property=InvocationID", "--value").stdout.strip()
    logs = command("journalctl", "--user", f"_SYSTEMD_INVOCATION_ID={invocation}",
                   "--no-pager", "-o", "cat").stdout
    urls = re.findall(r"https://[a-z0-9-]+\.trycloudflare\.com", logs)
    if urls:
        url = urls[-1]
        (PRIVATE / "url.txt").write_text(url + "\n")
        print(url)
        print(f"Login details: {PRIVATE / 'credentials.json'}")
        return url
    print("Tunnel is starting; run status again in a few seconds.")
    return None


def start() -> None:
    PRIVATE.mkdir(parents=True, exist_ok=True, mode=0o700)
    PRIVATE.chmod(0o700)
    caddy = ROOT / ".tools/linux/caddy"
    cloudflared = ROOT / ".tools/linux/cloudflared"
    if not caddy.is_file() or not cloudflared.is_file():
        raise SystemExit("Install the official Caddy and cloudflared binaries in .tools/linux first.")
    credentials = PRIVATE / "credentials.json"
    if not credentials.exists():
        with open(credentials, "x", opener=lambda path, flags: os.open(path, flags, 0o600)) as out:
            json.dump({"username": "biocoordinator", "password": secrets.token_urlsafe(24)}, out, indent=2)
    credentials.chmod(0o600)
    login = json.loads(credentials.read_text())
    hashed = subprocess.run([str(caddy), "hash-password"], input=login["password"] + "\n",
                            capture_output=True, text=True, check=True).stdout.strip()
    config = PRIVATE / "Caddyfile"
    config.write_text("{\n    admin off\n    auto_https off\n    persist_config off\n}\n"
                      "http://:8592 {\n    bind 127.0.0.1\n"
                      "    basic_auth {\n        " + login["username"] + " " + hashed + "\n    }\n"
                      "    header {\n        Cache-Control no-store\n"
                      "        X-Robots-Tag \"noindex, nofollow\"\n        -Server\n    }\n"
                      "    reverse_proxy 127.0.0.1:8502 {\n        header_up -Authorization\n    }\n}\n")
    config.chmod(0o600)
    command(str(caddy), "validate", "--config", str(config))
    try:
        urllib.request.urlopen("http://127.0.0.1:8502/_stcore/health", timeout=3).close()
    except OSError:
        launch("bio-coordinator-app", str(ROOT / "scripts/launch_ubuntu.sh"))
    launch(GATEWAY, str(caddy), "run", "--config", str(config))
    # Verify fail-closed authentication before exposing the gateway.
    for _ in range(20):
        try:
            urllib.request.urlopen("http://127.0.0.1:8592/", timeout=2).close()
            raise SystemExit("Gateway did not require authentication; tunnel was not started.")
        except urllib.error.HTTPError as error:
            if error.code == 401:
                break
            raise
        except urllib.error.URLError:
            time.sleep(0.25)
    else:
        raise SystemExit("Gateway did not become ready; tunnel was not started.")
    launch(TUNNEL, str(cloudflared), "tunnel", "--no-autoupdate", "--protocol", "http2",
           "--url", "http://127.0.0.1:8592")
    time.sleep(5)
    status()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["start", "status", "stop"])
    args = parser.parse_args()
    if args.action == "start":
        start()
    elif args.action == "stop":
        command("systemctl", "--user", "stop", TUNNEL, GATEWAY, check=False)
        print("Remote access stopped. Local app and experiment workers are unchanged.")
    else:
        status()

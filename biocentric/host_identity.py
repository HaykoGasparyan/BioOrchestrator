"""Worker identity is scoped to a machine, boot, and process incarnation."""
from __future__ import annotations

import hashlib
import os
import platform
import subprocess
import uuid
from pathlib import Path


_PROCESS_TOKEN = uuid.uuid4().hex


def host_identity() -> tuple[str, str]:
    try:
        machine = Path("/etc/machine-id").read_text().strip()
    except OSError:
        machine = f"{platform.system()}:{platform.node()}:{uuid.getnode()}"
    try:
        boot = Path("/proc/sys/kernel/random/boot_id").read_text().strip()
    except OSError:
        boot = ""
        if platform.system() == "Darwin":
            try:
                boot = subprocess.check_output(["/usr/sbin/sysctl", "-n", "kern.boottime"], text=True, timeout=2).strip()
            except (OSError, subprocess.SubprocessError):
                pass
    # An unverified boot can never authorize a PID lookup in another process.
    return hashlib.sha256(machine.encode()).hexdigest(), boot or _PROCESS_TOKEN


def process_identity(pid: int | None) -> str | None:
    if not isinstance(pid, int) or isinstance(pid, bool) or pid <= 0:
        return None
    try:
        # comm may contain spaces and parentheses; starttime is field 22.
        fields = Path(f"/proc/{pid}/stat").read_text().rsplit(")", 1)[1].split()
        if fields[0] == "Z":
            return None
        return fields[19]
    except (OSError, IndexError):
        # Safe portable fallback only recognizes this exact process.
        return _PROCESS_TOKEN if pid == os.getpid() else None

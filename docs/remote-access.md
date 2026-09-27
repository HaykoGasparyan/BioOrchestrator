# Remote access on this machine

Bio Orchestrator runs locally on port 8502. A Caddy gateway bound only to
127.0.0.1:8592 requires a password for every request, including WebSockets and
downloads. Cloudflare Tunnel publishes that gateway over HTTPS without opening
router ports. The application and experiment files remain on this machine;
remote traffic passes through Cloudflare.

From the project directory:

```bash
.venv-linux/bin/python scripts/remote_access.py start
.venv-linux/bin/python scripts/remote_access.py status
.venv-linux/bin/python scripts/remote_access.py stop
```

`status` prints the current URL. Login details are in
`.biocentric/remote-access/credentials.json`, readable only by your Linux user.
The login grants access to the same shared experiment library and controls as
the local app. Use the HTTPS URL when connecting remotely.

Gateway and tunnel run as background systemd user services with restart on
failure. Closing a browser or this chat does not stop them. They are temporary
services: after reboot, run `start` again. Stopping remote access leaves the local
app and experiment workers running. Keep this machine awake and online.

This uses a Cloudflare Quick Tunnel. Its URL can change when the tunnel process
restarts; run `status` for the current address. A fixed hostname requires a named
tunnel configured with a Cloudflare account/domain. Quick Tunnels have no uptime
guarantee; see [Cloudflare documentation](https://developers.cloudflare.com/cloudflare-one/networks/connectors/cloudflare-tunnel/do-more-with-tunnels/trycloudflare/).

Installed binaries are `.tools/linux/cloudflared` (2026.9.3) and
`.tools/linux/caddy` (2.11.4), downloaded from their official GitHub releases and
verified against the published SHA-256 asset digests. Neither secrets nor runtime
configuration belongs in source control; `.biocentric/` and `.tools/` are ignored.

# Preserving the copied Mac evidence on Ubuntu

The copied store contains seven projects and 1,508 registered immutable artifacts. Both requested completed projects are present: `6a2b1de5ee184b1b8c03b35249ed2057` and `726fbd0fea4e4b7496088bc3e458fd5b`. The supplied `runs/diverse-docking-verification/` and `runs/interactive-workspace-verification/` are historical Mac evidence, not new Ubuntu calculations.

Before any database migration, `runs/ubuntu-history-backup-20260926/` preserved a consistent SQLite snapshot, original artifact/job records, an archive of every store evidence file, individual SHA256 inventory, and archive checksums. A second backup, `runs/ubuntu-preserved-reference-runs-20260926/`, additionally includes both historical validation run directories. Neither archive is needed during normal operation; keep a copy off this disk for independent disaster recovery.

The explicit mapping in `.biocentric/relocation.json` resolves `/Users/hayk/Documents/BioGPT` to `/home/sargis/BioGPT/BioGPT`. Configuration checked all 1,508 artifact hashes and project content-addressed paths before activation. Original SQLite artifact paths, project JSON, events, manifests, scientific checkpoints and their recorded hashes were not rewritten. `Store.artifacts()` returns a local `path` and preserves `original_path`; each download still authenticates its bytes. The snapshot also exposes the root mapping. Historical JSON downloads deliberately retain their original paths and remain verifiable evidence.

For another move, configure a reviewed mapping on the untouched copied store before starting the app:

```sh
python3 -m biocentric.relocation --store .biocentric \
  --old-root /Users/hayk/Documents/BioGPT \
  --new-root /absolute/new/project \
  --backup /absolute/new/backup-directory
```

The command refuses an existing mapping, broad filesystem roots, traversal, escaping symlinks, missing or mismatched artifacts. Failure leaves the backup and never activates a partial mapping. `RelocationResolver.verified(original_path, recorded_sha256)` can locate historical scientific files without modifying their contents. Store artifact lookups additionally require the selected project and its registered content address.

Job claims now record machine identity, boot identity and process start identity. Legacy or foreign-host PIDs and future leases are not authoritative; no signal, including a PID-existence signal, is sent to them. Claims from the wrong host/boot cannot heartbeat or finish. A reused local PID cannot retain an expired job. Cancelling remains a durable database request observed by the owning worker. Existing database schema upgrades first copy SQLite and artifact metadata into `.biocentric/metadata-backups/`.

A relocated store creates new calculations in `projects/PROJECT/work/host-HASH/`, leaving historical work and checkpoint directories untouched. Completed projects remain completed. New Ubuntu validation uses separately named `runs/ubuntu-*` directories and current executable/source fingerprints.

Verification evidence is `runs/ubuntu-relocation-verification/verification.json`, including all artifact checks, unchanged historical database tables, and real CSV/SDF exports for both requested projects. The focused relocation/store/review/download suite passed 57 tests on Linux. It also exposed a real porting issue: immediate equal-length file edits may retain both nanosecond mtime and ctime on this filesystem. Downloads now rehash bytes on every read rather than trusting cached metadata; chemistry parsing remains cached by content.

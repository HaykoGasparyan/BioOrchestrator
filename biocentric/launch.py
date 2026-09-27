"""Local launcher with a private, memory-only OpenAI credential prompt."""
from __future__ import annotations

import argparse
import getpass
import os
import socket
import subprocess
import sys
import warnings
from pathlib import Path


class LaunchError(RuntimeError):
    pass


def read_key(environment: dict[str, str], *, interactive: bool, prompt=None) -> str:
    key = environment.get('OPENAI_API_KEY', '').strip()
    if not key:
        if not interactive:
            raise LaunchError('Run this command in your own interactive Terminal to enter the key privately.')
        # Never let getpass fall back to echoing a credential on an unsupported terminal.
        with warnings.catch_warnings():
            warnings.simplefilter('error', getpass.GetPassWarning)
            try:
                key = (prompt or getpass.getpass)('OpenAI API key (hidden; not saved): ').strip()
            except (getpass.GetPassWarning, EOFError) as exc:
                raise LaunchError('A terminal with hidden password entry is required.') from None
    if not key or any(character.isspace() for character in key):
        raise LaunchError('A nonempty API key without whitespace is required.')
    return key


def verify_connection(key: str, model: str, client_factory=None) -> dict:
    if client_factory is None:
        from openai import OpenAI
        client_factory = OpenAI
    try:
        with client_factory(api_key=key, base_url='https://api.openai.com/v1', timeout=30, max_retries=0) as client:
            response = client.responses.create(model=model, input='Reply with exactly CONNECTED.',
                                               max_output_tokens=16, store=False)
        if response.status != 'completed' or response.output_text.strip() != 'CONNECTED':
            raise LaunchError('The connection returned an incomplete or unexpected test response. No example was started.')
        return {'status': 'connected', 'model': model}
    except LaunchError:
        raise
    except Exception as exc:
        # Provider messages/bodies can contain request data. Never echo them or a traceback.
        status = getattr(exc, 'status_code', None)
        suffix = f', HTTP {status}' if isinstance(status, int) else ''
        raise LaunchError(f'OpenAI check failed ({type(exc).__name__}{suffix}). Check your API key, model access, API billing and internet connection.') from None


def child_environment(key: str, model: str, root: Path) -> dict[str, str]:
    environment = os.environ.copy()
    environment.update(OPENAI_API_KEY=key, OPENAI_MODEL=model,
                       OPENAI_BASE_URL='https://api.openai.com/v1',
                       PYTHONPATH=str(root), OPENAI_LOG='error')
    return environment


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--with-openai', action='store_true', required=True,
                        help='Ask privately for a key and make one small paid connection request.')
    parser.add_argument('--example', action='store_true',
                        help='Start a new public 1IEP research example with at most 3 advisory model calls; scientific review gates remain in place.')
    parser.add_argument('--model', default=os.getenv('OPENAI_MODEL', 'gpt-4.1-mini'))
    parser.add_argument('--port', type=int, default=8502,
                        help='Separate port to preserve the already running local-only app (default 8502).')
    args = parser.parse_args(argv)
    if not 1024 <= args.port <= 65535:
        parser.error('Choose a port between 1024 and 65535.')
    root = Path(__file__).resolve().parents[1]
    try:
        try:
            with socket.create_connection(('127.0.0.1', args.port), timeout=0.3):
                raise LaunchError(f'Port {args.port} is already in use; choose a different --port.')
        except ConnectionRefusedError:
            pass
        print(f'Model: {args.model}. One short API connection check will be billed to your API account.')
        if args.example:
            print('The public 1IEP example will retrieve public evidence and allow up to 3 model calls, each capped at 1,500 output tokens. Scientific approval gates remain enabled.')
        key = read_key(os.environ, interactive=sys.stdin.isatty())
        verify_connection(key, args.model)
        print('OpenAI connection verified. Credential held only in process memory.')
        environment = child_environment(key, args.model, root)
        del key
        project = None
        if args.example:
            from .service import Service
            service = Service(environment.get('BIOCENTRIC_HOME', str(root / '.biocentric')))
            project = service.create_example(allow_model=True, model=args.model)
            environment['BIOCENTRIC_HOME'] = str(service.root)
            print(f'Public model-enabled project created: {project}')
        command = [sys.executable, '-m', 'streamlit', 'run', str(root / 'app.py'),
                   '--server.address', '127.0.0.1', '--server.port', str(args.port),
                   '--server.headless', 'true', '--browser.gatherUsageStats', 'false']
        server = subprocess.Popen(command, cwd=root, env=environment)
        try:
            if project:
                # The short-lived CLI child launches the durable worker with the same private environment.
                subprocess.run([sys.executable, '-m', 'biocentric', '--root', environment['BIOCENTRIC_HOME'],
                                'run', project], cwd=root, env=environment, check=True)
            print(f'Open http://127.0.0.1:{args.port} . Keep this terminal running. Ctrl+C stops the web interface; use Cancel in the app to stop a worker.')
            return server.wait()
        except KeyboardInterrupt:
            server.terminate()
            server.wait(timeout=10)
            return 0
        except Exception:
            server.terminate()
            raise
    except KeyboardInterrupt:
        return 0
    except LaunchError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    except Exception as exc:
        print(f'Launcher stopped ({type(exc).__name__}); no credential or response body is printed.', file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())

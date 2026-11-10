#!/usr/bin/env python
"""Runs a Python script exactly as `python <script> [arguments ...]` would, and writes every setting it runs with.

Usage:
  python run_with_settings.py --settings <file.json> [--trace <module>:<function> ...] -- <script> [arguments ...]

<file.json> holds:
  command            this process's command line, the script with its sha256 and arguments, folder, host, CPUs
  options            the script's options as its own argparse parser resolved them, defaults included;
                     options_from_defaults lists those not given on the command line
  script_constants   the script's module-level UPPER_CASE constants, read from its source
  calls              every call of a --trace function: its arguments bound to the function's signature, with the
                     library's defaults filled in for those not passed (named in library_defaults). Data are described
                     by type and shape. Each call is recorded just before it runs.
  subprocesses       every command the script starts with subprocess.run. For pySCENIC's arboreto_with_multiprocessing.py
                     (GRNBoost2) also its options as its own parser resolves them, the call it makes for each target
                     gene, and every parameter of the regressor it fits.
  inputs             size, modification time and (up to 256 MB) sha256 of every existing file named in the arguments;
                     the entries of every folder named
  python, packages, all_packages, environment
                     interpreter; versions of the main packages and of every installed package; thread-count, cache,
                     temporary-folder and SLURM variables
  started, finished, seconds, exit_code, status ("running", "succeeded" or "failed"), errors

The script runs in this process (runpy, as __main__, with its own folder first on sys.path). Each --trace function is
replaced, in every module that refers to it, by a wrapper that records the call and then calls the function with the
same arguments, returning its result unchanged; worker processes the script starts call through without recording.
The file is rewritten (atomically) at every event, so a run that is killed keeps what it had reached. A problem while
recording a detail goes to "errors" and does not stop the script; a --trace that cannot be installed stops the run
before the script starts.

The file is written as <file.json>.partial.<pid> and renamed to <file.json> when the script exits with status 0. An
existing <file.json> is never overwritten: the .partial file is then kept and a message goes to stderr. Nothing else
is printed.
"""
import argparse
import ast
import datetime
import functools
import hashlib
import importlib
import importlib.metadata
import importlib.util
import inspect
import json
import numbers
import os
import platform
import runpy
import socket
import subprocess
import sys
import time
import traceback

SHA256_LIMIT = 256 * 2 ** 20  # files above this size are not hashed
MAX_CALLS = 20  # calls listed per traced function; later calls are only counted
KEY_PACKAGES = ("pyscenic", "arboreto", "ctxcore", "scikit-learn", "numpy", "scipy", "pandas", "anndata", "scanpy",
                "loompy", "h5py", "geosketch", "dask", "distributed", "numba", "matplotlib")
ENV_KEYS = ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS", "NUMBA_NUM_THREADS",
            "PYTHONHASHSEED", "HDF5_USE_FILE_LOCKING", "TMPDIR", "MPLCONFIGDIR", "NUMBA_CACHE_DIR", "XDG_CACHE_HOME",
            "SLURM_JOB_ID", "SLURM_ARRAY_JOB_ID", "SLURM_ARRAY_TASK_ID", "SLURM_JOB_NAME", "SLURM_JOB_PARTITION",
            "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE", "SLURM_JOB_NODELIST", "SLURM_SUBMIT_DIR")
ARBORETO = "arboreto_with_multiprocessing.py"


def now():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 24), b""):
            h.update(block)
    return h.hexdigest()


def type_name(value):
    t = type(value)
    return t.__qualname__ if t.__module__ == "builtins" else f"{t.__module__}.{t.__qualname__}"


def summarize(value, depth=0):
    """A JSON-ready description: scalars as they are, short containers item by item, data by type and shape."""
    if value is None or isinstance(value, (bool, str)):
        return value
    if isinstance(value, numbers.Integral):
        return int(value)
    if isinstance(value, numbers.Real):
        return float(value)
    if getattr(value, "shape", None) == () and hasattr(value, "item"):  # a numpy scalar
        return summarize(value.item(), depth)
    if isinstance(value, (bytes, bytearray)):
        return {"type": type_name(value), "bytes": len(value)}
    if inspect.isfunction(value) or inspect.isbuiltin(value) or inspect.ismethod(value) or inspect.isclass(value):
        return {"type": type_name(value),
                "name": f"{getattr(value, '__module__', '')}.{getattr(value, '__qualname__', repr(value))}"}
    if hasattr(value, "read") and hasattr(value, "name") and hasattr(value, "mode"):  # an open file (argparse.FileType)
        return {"type": "file", "name": str(value.name), "mode": str(value.mode)}
    if isinstance(value, argparse.Namespace):
        value = vars(value)
    shape = getattr(value, "shape", None)
    if shape is not None and not isinstance(value, (dict, list, tuple)):
        out = {"type": type_name(value)}
        try:
            out["shape"] = [int(s) for s in shape]
        except TypeError:
            out["shape"] = str(shape)
        dtype = getattr(value, "dtype", None)
        if dtype is not None and not hasattr(dtype, "index"):
            out["dtype"] = str(dtype)
        return out
    if isinstance(value, dict):
        if len(value) <= 50 and depth < 4:
            return {str(k): summarize(v, depth + 1) for k, v in value.items()}
        return {"type": type_name(value), "len": len(value)}
    if isinstance(value, (list, tuple, set, frozenset)):
        items = sorted(value, key=repr) if isinstance(value, (set, frozenset)) else list(value)
        simple = all(x is None or isinstance(x, (bool, str, numbers.Number)) for x in items)
        if depth < 4 and ((simple and len(items) <= 50) or len(items) <= 10):
            return [summarize(x, depth + 1) for x in items]
        out = {"type": type_name(value), "len": len(items)}
        if items:
            out["first"] = summarize(items[0], depth + 1)
        return out
    text = repr(value)
    return {"type": type_name(value), "repr": text if len(text) <= 300 else text[:300] + "..."}


def bind(signature, args, kwargs):
    """(arguments with the defaults filled in, names of those not passed) of one call."""
    bound = signature.bind(*args, **kwargs)
    passed = set(bound.arguments)
    bound.apply_defaults()
    return ({k: summarize(v) for k, v in bound.arguments.items()},
            [k for k in bound.arguments if k not in passed])


def script_constants(path):
    """Module-level UPPER_CASE assignments of a script: literal values, or the expression as written."""
    source = open(path, encoding="utf-8").read()
    out = {}
    for node in ast.parse(source).body:
        if isinstance(node, ast.Assign):
            targets, value = node.targets, node.value
        elif isinstance(node, ast.AnnAssign) and node.value is not None:
            targets, value = [node.target], node.value
        else:
            continue
        for target in targets:
            if isinstance(target, ast.Name) and target.id.upper() == target.id and target.id[0].isalpha():
                try:
                    out[target.id] = summarize(ast.literal_eval(value))
                except (ValueError, TypeError, SyntaxError, MemoryError, RecursionError):
                    out[target.id] = {"expression": (ast.get_source_segment(source, value) or "")[:300]}
    return out


def call_keywords(source, function, namespace):
    """Keyword arguments of the first call of `function` written in `source`, resolved where possible."""
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Call) and getattr(node.func, "id", getattr(node.func, "attr", None)) == function:
            out = {}
            for kw in node.keywords:
                text = ast.get_source_segment(source, kw.value) or ""
                try:
                    out[kw.arg] = summarize(ast.literal_eval(kw.value))
                except (ValueError, TypeError, SyntaxError):
                    out[kw.arg] = summarize(namespace[text]) if text in namespace else {"expression": text}
            return out
    return None


class Recorder:
    def __init__(self, path, script, script_args, traces):
        self.path = os.path.abspath(path)
        self.pid = os.getpid()
        self.partial = f"{self.path}.partial.{self.pid}"
        self.t0 = time.time()
        self.internal = False  # set while the recorder parses arguments itself
        self.counts = {}
        try:
            cpus = len(os.sched_getaffinity(0))
        except AttributeError:
            cpus = None
        self.rec = {
            "status": "running",
            "started": now(),
            "command": list(getattr(sys, "orig_argv", None) or [sys.executable] + sys.argv),
            "script": {"path": os.path.abspath(script), "arguments": list(script_args)},
            "working_folder": os.getcwd(),
            "host": socket.gethostname(),
            "cpus": {"available_to_this_process": cpus, "on_the_node": os.cpu_count()},
            "traced": {spec: None for spec in traces},
            "options": None,
            "calls": [],
            "subprocesses": [],
            "errors": [],
        }

    def write(self):
        tmp = f"{self.partial}.tmp"
        with open(tmp, "w") as fh:
            json.dump(self.rec, fh, indent=1, default=str)
            fh.write("\n")
        os.replace(tmp, self.partial)

    def error(self, where):
        self.rec["errors"].append({"where": where, "error": traceback.format_exc(limit=4)[-2000:]})
        try:
            self.write()
        except Exception:
            pass

    def start(self):
        """Everything known before the script runs; written at once, so a settings file that cannot be written stops
        the run before anything is computed."""
        self.write()
        script = self.rec["script"]["path"]
        for key, fill in (("sha256", lambda: sha256(script)),
                          ("script_constants", lambda: script_constants(script)),
                          ("inputs", lambda: self.inputs(self.rec["script"]["arguments"])),
                          ("python", lambda: {"executable": sys.executable, "version": sys.version,
                                              "platform": platform.platform()}),
                          ("packages", self.packages),
                          ("environment", lambda: {k: os.environ[k] for k in ENV_KEYS if k in os.environ})):
            try:
                value = fill()
            except Exception:
                self.error(f"recording {key}")
                continue
            if key == "sha256":
                self.rec["script"]["sha256"] = value
            elif key == "packages":
                self.rec["packages"], self.rec["all_packages"] = value
            else:
                self.rec[key] = value
        self.write()

    @staticmethod
    def inputs(arguments):
        out, seen = [], set()
        for arg in arguments:
            for cand in (arg, arg.split("=", 1)[1] if arg.startswith("-") and "=" in arg else None):
                if not cand or cand in seen or not os.path.exists(cand):
                    continue
                seen.add(cand)
                path = os.path.abspath(cand)
                if os.path.isdir(path):
                    entries = sorted(os.listdir(path))
                    out.append({"path": path, "type": "folder", "entries": len(entries), "names": entries[:50]})
                    continue
                st = os.stat(path)
                entry = {"path": path, "bytes": st.st_size,
                         "modified": datetime.datetime.fromtimestamp(st.st_mtime).astimezone().isoformat(timespec="seconds")}
                if st.st_size <= SHA256_LIMIT:
                    entry["sha256"] = sha256(path)
                out.append(entry)
        return out

    @staticmethod
    def packages():
        every = {}
        for dist in importlib.metadata.distributions():
            name = dist.metadata["Name"]
            if name:
                every.setdefault(name, dist.version)
        key = {}
        for name in KEY_PACKAGES:
            try:
                key[name] = importlib.metadata.version(name)
            except importlib.metadata.PackageNotFoundError:
                pass
        return key, dict(sorted(every.items(), key=lambda kv: kv[0].lower()))

    def options(self, parser, namespace, args):
        argv = [str(a) for a in (sys.argv[1:] if args is None else args)]
        values = vars(namespace)
        if self.rec["options"] is not None:  # a later parse: kept apart from the script's own options
            self.rec.setdefault("later_parses", []).append(
                {"parser": parser.prog, "options": {k: summarize(v, 1) for k, v in values.items()}})
            self.write()
            return
        actions = list(parser._actions)
        for action in parser._actions:
            if isinstance(action, argparse._SubParsersAction):
                for name, sub in action.choices.items():
                    if name in argv:
                        self.rec["subcommand"] = name
                        actions += sub._actions
        default = {a.dest: a.default for a in actions}
        given = set()
        for a in actions:
            if not a.option_strings:
                given.add(a.dest)
            elif any(t == s or (s.startswith("--") and t.startswith(s + "=")) for t in argv for s in a.option_strings):
                given.add(a.dest)
        self.rec["options_parser"] = parser.prog
        self.rec["options"] = {k: summarize(v, 1) for k, v in values.items()}
        self.rec["options_from_defaults"] = sorted(
            k for k, v in values.items()
            if k not in given and (k not in default or summarize(v, 1) == summarize(default[k], 1)))
        self.write()

    def call(self, label, signature, args, kwargs):
        n = self.counts[label] = self.counts.get(label, 0) + 1
        if n > MAX_CALLS:
            return
        entry = {"function": label, "call": n, "time": now()}
        try:
            if signature is None:
                raise TypeError("no signature")
            entry["arguments"], entry["library_defaults"] = bind(signature, args, kwargs)
        except TypeError as e:  # the function itself will raise on such a call
            entry["not_bound"] = str(e)
            entry["args"] = [summarize(a) for a in args]
            entry["kwargs"] = {k: summarize(v) for k, v in kwargs.items()}
        self.rec["calls"].append(entry)
        self.write()

    def subprocess_run(self, popenargs, kwargs):
        cmd = popenargs[0] if popenargs else kwargs.get("args")
        entry = {"time": now(), "command": [str(x) for x in cmd] if isinstance(cmd, (list, tuple)) else str(cmd)}
        if isinstance(cmd, (list, tuple)):
            hits = [i for i, x in enumerate(cmd) if os.path.basename(str(x)) == ARBORETO]
            if hits:
                try:
                    entry["grnboost2"] = self.arboreto(str(cmd[hits[0]]), [str(x) for x in cmd[hits[0] + 1:]])
                except Exception:
                    self.error(f"recording the settings of {ARBORETO}")
        self.rec["subprocesses"].append(entry)
        self.write()

    def arboreto(self, script, argv):
        """Settings of an arboreto_with_multiprocessing.py command, from the script itself: its parser, its constants,
        the keywords of its infer_partial_network() call, and the regressor arboreto.core.fit_model builds from them.
        Importing the script runs only its definitions (its main() runs only as __main__)."""
        spec = importlib.util.spec_from_file_location("_arboreto_settings", script)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        parser = mod.create_argument_parser()
        self.internal = True
        try:
            ns = parser.parse_args(argv)
        except SystemExit as e:
            raise RuntimeError(f"its parser rejected the arguments (exit {e.code})") from None
        finally:
            self.internal = False
        import arboreto.core as core
        regressor_type, kwargs_name = {"grnboost2": ("GBM", "SGBM_KWARGS"), "genie3": ("RF", "RF_KWARGS")}[ns.method]
        regressor_kwargs = dict(getattr(mod, kwargs_name))
        cls = core.SKLEARN_REGRESSOR_FACTORY[regressor_type]
        parameters, library_defaults = bind(inspect.signature(cls), (), dict(random_state=ns.seed, **regressor_kwargs))
        defaults = {a.dest: a.default for a in parser._actions}
        values = vars(ns)
        source = open(script, encoding="utf-8").read()
        return {
            "script": script,
            "script_sha256": sha256(script),
            "options": {k: summarize(v, 1) for k, v in values.items()},
            "options_from_defaults": sorted(k for k, v in values.items() if k in defaults and v == defaults[k]),
            "per_target_gene": {
                "function": "arboreto.core.infer_partial_network, once for every gene of the loom (Pool.imap)",
                "regressor_type": regressor_type,
                "regressor_kwargs": summarize(regressor_kwargs),
                "seed": ns.seed,
                "keywords_as_written": call_keywords(source, "infer_partial_network", vars(mod)),
            },
            "regressor": {
                "class": f"{cls.__module__}.{cls.__qualname__}",
                "built_as": "SKLEARN_REGRESSOR_FACTORY[regressor_type](random_state=seed, **regressor_kwargs)"
                            " in arboreto.core.fit_model",
                "parameters": parameters,
                "library_defaults": library_defaults,
            },
            "early_stopping": {
                "used": bool(core.is_oob_heuristic_supported(regressor_type, regressor_kwargs)),
                "window_length": summarize(mod.EARLY_STOP_WINDOW_LENGTH),
                "monitor": "arboreto.core.EarlyStopMonitor, passed to the regressor's fit()",
            },
            "arboreto_core": {"file": core.__file__, "sha256": sha256(core.__file__)},
        }

    def finish(self, code):
        exit_code = 0 if code is None else code if isinstance(code, int) else 1
        self.rec.update(finished=now(), seconds=round(time.time() - self.t0, 1), exit_code=exit_code,
                        status="succeeded" if exit_code == 0 else "failed")
        if any(n > MAX_CALLS for n in self.counts.values()):
            self.rec["calls_per_function"] = dict(self.counts)
        self.write()
        if exit_code == 0:
            if os.path.exists(self.path):
                print(f"run_with_settings: {self.path} exists, so the settings stay in {self.partial}", file=sys.stderr)
            else:
                os.replace(self.partial, self.path)


def install_trace(recorder, spec):
    """Replaces <module>:<function> in every loaded module that refers to it by a recording wrapper."""
    module_name, _, attr = spec.partition(":")
    if not module_name or not attr:
        raise ValueError(f"--trace {spec}: expected <module>:<function>")
    module = importlib.import_module(module_name)
    original = getattr(module, attr)
    label = f"{module_name}.{attr}"
    try:
        signature = inspect.signature(original)
    except (TypeError, ValueError):
        signature = None

    @functools.wraps(original)
    def traced(*args, **kwargs):
        if os.getpid() == recorder.pid and not recorder.internal:
            try:
                recorder.call(label, signature, args, kwargs)
            except Exception:
                recorder.error(f"recording a call of {label}")
        return original(*args, **kwargs)

    replaced = 0
    for mod in list(sys.modules.values()):
        try:
            items = list(vars(mod).items())
        except TypeError:
            continue
        for name, value in items:
            if value is original:
                setattr(mod, name, traced)
                replaced += 1
    recorder.rec["traced"][spec] = {"defined_in": f"{original.__module__}.{original.__qualname__}",
                                    "signature": str(signature), "references_replaced": replaced}


def main():
    argv = sys.argv[1:]
    if "--" not in argv or argv.index("--") == len(argv) - 1:
        sys.exit("usage: run_with_settings.py --settings <file.json> [--trace <module>:<function> ...] "
                 "-- <script> [arguments ...]")
    cut = argv.index("--")
    ap = argparse.ArgumentParser(prog="run_with_settings.py")
    ap.add_argument("--settings", required=True, help="the JSON file to write")
    ap.add_argument("--trace", action="append", default=[], metavar="MODULE:FUNCTION",
                    help="record every call of this function (repeatable)")
    own = ap.parse_args(argv[:cut])
    script, script_args = argv[cut + 1], argv[cut + 2:]

    recorder = Recorder(own.settings, script, script_args, own.trace)
    recorder.start()
    for spec in own.trace:
        try:
            install_trace(recorder, spec)
        except Exception:
            recorder.error(f"installing --trace {spec}")
            recorder.finish(2)
            print(f"run_with_settings: cannot trace {spec} (see {recorder.partial}); the script was not run",
                  file=sys.stderr)
            sys.exit(2)
    recorder.write()

    parse_args = argparse.ArgumentParser.parse_args

    def recording_parse_args(self, args=None, namespace=None):
        ns = parse_args(self, args, namespace)
        if os.getpid() == recorder.pid and not recorder.internal:
            try:
                recorder.options(self, ns, args)
            except Exception:
                recorder.error("recording the options")
        return ns

    run = subprocess.run

    def recording_run(*popenargs, **kwargs):
        if os.getpid() == recorder.pid:
            try:
                recorder.subprocess_run(popenargs, kwargs)
            except Exception:
                recorder.error("recording a subprocess")
        return run(*popenargs, **kwargs)

    argparse.ArgumentParser.parse_args = recording_parse_args
    subprocess.run = recording_run

    # As `python <script>`: the script's folder first on sys.path; runpy sets sys.argv[0] to the script.
    sys.argv = [script] + script_args
    sys.path[0] = os.path.dirname(os.path.abspath(script))
    code = 0
    try:
        runpy.run_path(script, run_name="__main__")
    except SystemExit as e:
        code = e.code
        raise
    except BaseException:
        code = 1
        recorder.rec["exception"] = traceback.format_exc(limit=6)[-4000:]
        raise
    finally:
        recorder.finish(code)


if __name__ == "__main__":
    main()

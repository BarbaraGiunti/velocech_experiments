"""
measure_run.py
==============

Runs ONE command with an optional time limit and measures, from outside the
program, its wall-clock time and peak memory. Appends one JSON line with the
result to a record file. Used by run_all_experiments.sh for every run.

    python measure_run.py --record runs.jsonl [--timeout SECONDS]
                          [--stdout FILE] [--stderr FILE] [--set key=value ...]
                          -- COMMAND ARG ...

Record (one JSON line): the --set values, plus
    status          ok | timeout | error
    exit_code       the command's exit code (None if it was stopped)
    wall_time_s     from starting the command to its end (includes starting
                    the process, reading inputs and writing outputs)
    peak_memory_mb  peak memory of the command (see below)
    memory_method   how the peak memory was obtained
    timeout_s, started

How it is measured (nothing runs inside the measured program, and this
helper does not poll it -- it only waits, so it does not disturb the timing):
  * Linux / macOS: the command runs in its own process group; the helper
    blocks in wait4(), which returns the operating system's resource usage
    of the finished command: ru_maxrss = peak RESIDENT memory of the largest
    process among the command and the child processes it waited for (e.g.
    alpha_at_half_diameter.py and the alphadual it starts). The time limit is
    an interval timer: at the limit the whole process group gets SIGTERM, and
    SIGKILL 5 s later if it is still running.
  * Windows: the command is placed in a Windows "job object"; after it ends,
    the job reports PeakProcessMemoryUsed = peak COMMITTED memory of the
    largest process in the job (the command and all processes it started).
    At the time limit the whole job is terminated. If the job cannot be
    created, the command's own peak (PeakPagefileUsage) is used and the time
    limit uses "taskkill /T /F".
The exit code of this helper is 0 whatever happened to the command (the
outcome is in the record), so a calling loop simply continues.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import signal
import subprocess
import sys
import time


def run_posix(cmd, timeout, out, err):
    state = {"timed_out": False}
    t0 = time.perf_counter()
    p = subprocess.Popen(cmd, stdout=out, stderr=err, start_new_session=True)

    def on_alarm(signum, frame):
        try:
            if not state["timed_out"]:
                state["timed_out"] = True
                os.killpg(p.pid, signal.SIGTERM)
                signal.setitimer(signal.ITIMER_REAL, 5.0)     # SIGKILL 5 s later
            else:
                os.killpg(p.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass

    if timeout and timeout > 0:
        signal.signal(signal.SIGALRM, on_alarm)
        signal.setitimer(signal.ITIMER_REAL, timeout)
    _, status, ru = os.wait4(p.pid, 0)          # blocks until the command ends
    wall = time.perf_counter() - t0
    signal.setitimer(signal.ITIMER_REAL, 0)
    p.returncode = os.waitstatus_to_exitcode(status) if hasattr(os, "waitstatus_to_exitcode") \
        else (os.WEXITSTATUS(status) if os.WIFEXITED(status) else -os.WTERMSIG(status))
    if state["timed_out"]:
        try:
            os.killpg(p.pid, signal.SIGKILL)    # leftovers of the process group
        except ProcessLookupError:
            pass
    scale = 1048576.0 if sys.platform == "darwin" else 1024.0   # bytes on macOS, KB on Linux
    return (state["timed_out"], p.returncode, wall, ru.ru_maxrss / scale,
            "peak resident memory (wait4 ru_maxrss), largest process")


def run_windows(cmd, timeout, out, err):
    import ctypes
    from ctypes import wintypes

    k32 = ctypes.WinDLL("kernel32", use_last_error=True)

    class IO_COUNTERS(ctypes.Structure):
        _fields_ = [(n, ctypes.c_ulonglong) for n in (
            "ReadOperationCount", "WriteOperationCount", "OtherOperationCount",
            "ReadTransferCount", "WriteTransferCount", "OtherTransferCount")]

    class BASIC_LIMIT(ctypes.Structure):
        _fields_ = [("PerProcessUserTimeLimit", ctypes.c_int64),
                    ("PerJobUserTimeLimit", ctypes.c_int64),
                    ("LimitFlags", wintypes.DWORD),
                    ("MinimumWorkingSetSize", ctypes.c_size_t),
                    ("MaximumWorkingSetSize", ctypes.c_size_t),
                    ("ActiveProcessLimit", wintypes.DWORD),
                    ("Affinity", ctypes.c_size_t),
                    ("PriorityClass", wintypes.DWORD),
                    ("SchedulingClass", wintypes.DWORD)]

    class EXTENDED_LIMIT(ctypes.Structure):
        _fields_ = [("BasicLimitInformation", BASIC_LIMIT),
                    ("IoInfo", IO_COUNTERS),
                    ("ProcessMemoryLimit", ctypes.c_size_t),
                    ("JobMemoryLimit", ctypes.c_size_t),
                    ("PeakProcessMemoryUsed", ctypes.c_size_t),
                    ("PeakJobMemoryUsed", ctypes.c_size_t)]

    class PROCESS_MEMORY_COUNTERS(ctypes.Structure):
        _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                    ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
                    ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                    ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t)]

    JobObjectExtendedLimitInformation = 9
    JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE = 0x2000
    k32.CreateJobObjectW.restype = wintypes.HANDLE
    k32.CreateJobObjectW.argtypes = [ctypes.c_void_p, wintypes.LPCWSTR]
    k32.SetInformationJobObject.argtypes = [wintypes.HANDLE, ctypes.c_int, ctypes.c_void_p, wintypes.DWORD]
    k32.QueryInformationJobObject.argtypes = [wintypes.HANDLE, ctypes.c_int, ctypes.c_void_p,
                                              wintypes.DWORD, ctypes.c_void_p]
    k32.AssignProcessToJobObject.argtypes = [wintypes.HANDLE, wintypes.HANDLE]
    k32.TerminateJobObject.argtypes = [wintypes.HANDLE, wintypes.UINT]
    k32.CloseHandle.argtypes = [wintypes.HANDLE]
    k32.K32GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.c_void_p, wintypes.DWORD]

    job = k32.CreateJobObjectW(None, None)
    if job:
        lim = EXTENDED_LIMIT()
        lim.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE
        k32.SetInformationJobObject(job, JobObjectExtendedLimitInformation,
                                    ctypes.byref(lim), ctypes.sizeof(lim))
    t0 = time.perf_counter()
    p = subprocess.Popen(cmd, stdout=out, stderr=err)
    in_job = bool(job) and bool(k32.AssignProcessToJobObject(job, wintypes.HANDLE(int(p._handle))))
    timed_out = False
    try:
        p.wait(timeout=timeout if timeout and timeout > 0 else None)   # blocks (no polling)
    except subprocess.TimeoutExpired:
        timed_out = True
        if in_job:
            k32.TerminateJobObject(job, 1)
        else:
            subprocess.run(["taskkill", "/T", "/F", "/PID", str(p.pid)], capture_output=True)
        p.wait()
    wall = time.perf_counter() - t0
    peak, method = None, ""
    if in_job:
        q = EXTENDED_LIMIT()
        if k32.QueryInformationJobObject(job, JobObjectExtendedLimitInformation,
                                         ctypes.byref(q), ctypes.sizeof(q), None):
            peak = q.PeakProcessMemoryUsed / 1048576.0
            method = "peak committed memory (job object), largest process"
    if peak is None:
        pmc = PROCESS_MEMORY_COUNTERS()
        pmc.cb = ctypes.sizeof(pmc)
        if k32.K32GetProcessMemoryInfo(wintypes.HANDLE(int(p._handle)), ctypes.byref(pmc), pmc.cb):
            peak = pmc.PeakPagefileUsage / 1048576.0
            method = "peak committed memory of the command only (no job object)"
    if job:
        k32.CloseHandle(job)
    return timed_out, (None if timed_out else p.returncode), wall, peak, method


def main():
    if "--" not in sys.argv:
        sys.exit("usage: measure_run.py --record FILE [--timeout S] [--stdout F] [--stderr F] "
                 "[--set key=value ...] -- COMMAND ...")
    k = sys.argv.index("--")
    ap = argparse.ArgumentParser()
    ap.add_argument("--record", required=True)
    ap.add_argument("--timeout", type=float, default=0.0)
    ap.add_argument("--stdout")
    ap.add_argument("--stderr")
    ap.add_argument("--set", action="append", default=[])
    args = ap.parse_args(sys.argv[1:k])
    cmd = sys.argv[k + 1:]
    if not cmd:
        sys.exit("measure_run.py: no command given after --")

    rec = {}
    for kv in args.set:
        key, _, val = kv.partition("=")
        rec[key] = val
    rec["started"] = dt.datetime.now().isoformat(timespec="seconds")
    rec["timeout_s"] = args.timeout or None
    out = open(args.stdout, "wb") if args.stdout else subprocess.DEVNULL
    err = open(args.stderr, "wb") if args.stderr else subprocess.DEVNULL
    try:
        if os.name == "nt":
            timed_out, code, wall, peak, method = run_windows(cmd, args.timeout, out, err)
        else:
            timed_out, code, wall, peak, method = run_posix(cmd, args.timeout, out, err)
        status = "timeout" if timed_out else ("ok" if code == 0 else "error")
    except FileNotFoundError as e:
        timed_out, code, wall, peak, method, status = False, None, 0.0, None, "", "error"
        rec["message"] = f"command not found: {e}"
    finally:
        for f in (out, err):
            if f is not subprocess.DEVNULL:
                f.close()
    rec.update(status=status, exit_code=code, wall_time_s=wall, peak_memory_mb=peak,
               memory_method=method)
    with open(args.record, "a") as fh:            # one line, forced to disk
        fh.write(json.dumps(rec) + "\n")
        fh.flush()
        os.fsync(fh.fileno())
    sys.exit(0)


if __name__ == "__main__":
    main()

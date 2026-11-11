"""Conservative CPU and memory limits for independent simulation workers."""

import os
from pathlib import Path


def available_memory_bytes() -> int:
    """Read available host memory, bounded by any readable cgroup limit."""
    memory = int(os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_AVPHYS_PAGES"))
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable:"):
                memory = int(line.split()[1]) * 1024
    except OSError:
        pass
    for name in ["/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes"]:
        try:
            limit = int(Path(name).read_text())
            if limit > 0:
                memory = min(memory, limit)
        except (OSError, ValueError):
            pass
    return memory


def worker_limit(requested: int, jobs: int) -> int:
    """Reserve half the CPUs and approximately 512 MiB per worker."""
    if requested < 1 or jobs < 1:
        raise ValueError("Worker count and trajectory count must be positive")
    cores = os.cpu_count() or 2
    if hasattr(os, "sched_getaffinity"):
        cores = min(cores, len(os.sched_getaffinity(0)))
    memory_workers = max(1, (available_memory_bytes() - 512 * 1024**2) // (512 * 1024**2))
    return min(requested, max(1, cores // 2), memory_workers, jobs)

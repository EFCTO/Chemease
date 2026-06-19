import multiprocessing as mp
from typing import Any, Callable, Iterable, List

from chemease.logger import logger


def parallelize(func: Callable, data: Iterable, n_jobs: int = -1) -> List[Any]:
    if n_jobs == -1:
        n_jobs = mp.cpu_count()

    logger.info(f"Running {func.__name__} in parallel with {n_jobs} cores.")

    with mp.Pool(processes=n_jobs) as pool:
        results = pool.map(func, data)

    return results

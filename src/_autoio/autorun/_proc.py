""" Process stuff
"""

import os
import random
import errno
import datetime
import signal
import multiprocessing
from functools import wraps
import numpy


def utc_time():
    """ get the current UTC time

    :returns: the utc time
    :rtype: datetime
    """
    return datetime.datetime.utcnow()


def set_nprocs(nobjs, nprocs='auto'):
    """ Set the number of processors to use for some task
    """
    if nprocs is None:
        _nprocs = 1
    elif nprocs == 'auto':
        _nprocs = len(os.sched_getaffinity(0)) - 1
    elif isinstance(nprocs, int):
        _nprocs = nprocs
    else:
        raise NotImplementedError

    # Set number of processors equal to obj number if more available
    _nprocs = min(_nprocs, nobjs)

    return _nprocs


def execute_function_in_parallel(fxn, objs, args, nprocs='auto'):
    """ MultiProcessing wrapper function that can execute some
        function using a set of arguments across multiple processors.
        Function will largely loop over a list of objects

        Note that the function to parallelize must take
        arguments in a specific manner.

        fxn(args, objs, output_queue)

        where:
        objs is a list of objects that the task will execute over
        output_queue is a variable for a multiprocessing.Queue()
    """

    # Allot number of objects for each processor
    num_obj = len(objs)
    rand_objs = numpy.empty((num_obj,), dtype=object)

    nprocs = set_nprocs(num_obj, nprocs=nprocs)
    # Randomize objects the objects (for distributing workload?)
    # This function may it very hard to debug - removed for now
    # obj_per_proc = math.ceiling(num_obj / nprocs)
    if nprocs > 1:
        rand_objs[:] = random.sample(objs, num_obj)
        print('Begin parallel job array on {:g} processors'.format(nprocs))
    else:
        # rand_objs[:] = objs
        rand_objs = objs
    rand_objs = numpy.array(rand_objs, dtype=object)
    rand_objs_splt = numpy.array_split(rand_objs, nprocs)

    # Loop over each processor and launch the process
    output_queue = multiprocessing.Queue()
    procs = []
    for proc_n in range(nprocs):

        # Generate list of objects to work with for each processor
        # obj_start = proc_n * obj_per_proc
        # obj_end = (proc_n+1) * obj_per_proc
        # if proc_n != nprocs-1 else num_obj
        # obj_lst = rand_objs[obj_start:obj_end]
        obj_lst = rand_objs_splt[proc_n]
        # Create full args to pass to function including
        fxn_args = tuple(args) + (obj_lst, output_queue)

        # Create a process object for each processor; Launch associated process
        proc = multiprocessing.Process(target=fxn, args=fxn_args)
        procs.append(proc)
        proc.start()

    # Collect the results of each process together
    output_lst = ()
    for _ in procs:
        output_lst += output_queue.get()

    # Wait for all processes to finish
    for proc in procs:
        proc.join()

    return output_lst


def timeout(seconds=10, error_message=os.strerror(errno.ETIME)):
    """ check if process has died
    """

    def decorator(func):
        def _handle_timeout():
            print(error_message)
            raise TimeoutError(error_message)

        def wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result

        return wraps(func)(wrapper)

    return decorator

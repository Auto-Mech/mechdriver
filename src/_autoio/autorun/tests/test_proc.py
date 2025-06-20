""" test the autorun._proc
"""

import autorun


def test__():
    """ test autorun._proc.multiproc_task
    """

    def _adder(num1, num2, num3, num_lst, output_queue=None):
        """ Test adder function that adds three nums in some combination
            to numbers from a list
        """
        added_nums = ()
        for num in num_lst:
            added_nums += (num + num1 + 2*num2 + 3*num3,)

        # Either set the return for direct fxn call or
        # Add results to multiprocessing.Queue obj passed in
        if output_queue is None:
            ret = added_nums
        else:
            ret = None
            output_queue.put(added_nums)

        return ret

    adder_input = (1, 5, 10)
    adder_args = (2, 4, 6)

    results1 = autorun.execute_function_in_parallel(
        _adder, adder_input, adder_args, nprocs='auto')
    results2 = autorun.execute_function_in_parallel(
        _adder, adder_input, adder_args, nprocs=1)
    assert set(results1) == set(results2) == {29, 33, 38}


if __name__ == '__main__':
    test__()

"""Defines abstract base class for automatically handling model units."""

import abc
import functools
from collections.abc import Sequence
from typing import ClassVar

import numpy
from numpy.typing import NDArray

from ..util.type_ import Frozen
from . import dim
from .dim import Dimension
from .system import UNITS, Units, UnitsData


class UnitManager(Frozen, abc.ABC):
    _dimensions: ClassVar[dict[str, Dimension]]

    def __init__(self, units: UnitsData | None = None, **kwargs):
        if units is not None:
            units0 = Units.model_validate(units)
            for key, dim_ in self.__class__._dimensions.items():
                val0 = kwargs.get(key)
                kwargs[key] = dim.convert(units0, UNITS, dim_, val0, **kwargs)

        super().__init__(**kwargs)


def manage_units(arg_dims: Sequence[Dimension], ret_dim: Dimension | None = None):
    """Transform function into a unit managing function.

    TODO: Fix handling of args vs kwargs for input. This is currently broken when the
    user passes in e.g. P=1 for the first argument (this gets passed in kwargs, not
    args). To resolve this, have `arg_dims` be a dict of {name: dimension}

    Converts arguments to internal units, calls the function, then converts the return
    value back to desired units.

    :param arg_dims: Argument dimensions
    :param ret_dim: Return dimensions
    :return: Function decorator
    """

    def manage_units_(func0):
        @functools.wraps(func0)
        def func(
            self, *args, units: UnitsData | None = None, **kwargs
        ) -> float | NDArray[numpy.float64]:
            # If no units were specified, return as-is
            if units is None:
                return func0(self, *args, units=units, **kwargs)

            # Process units
            units0 = Units.model_validate(units)

            # Convert arguments
            args_ = tuple(
                dim.convert(units0, UNITS, dim_, arg, **self.model_dump())
                for dim_, arg in zip(arg_dims, args, strict=True)
            )

            # Call function
            ret = func0(self, *args_, units=units, **kwargs)

            # Convert return
            if ret_dim is not None:
                ret = dim.convert(UNITS, units0, ret_dim, ret, **self.model_dump())

            return ret

        return func

    return manage_units_

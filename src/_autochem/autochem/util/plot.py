"""Plotting helpers."""

import itertools
from collections.abc import Sequence
from typing import Literal, Protocol

import altair
import numpy
import pandas
from numpy.typing import ArrayLike, NDArray

from .. import unit_
from ..unit_ import UNITS, Units, UnitsData


class RateFunction(Protocol):
    """Protocol for callable rate function."""

    def __call__(
        self,
        T: ArrayLike,  # noqa: N803
        P: ArrayLike = 1,  # noqa: N803
        units: UnitsData | None = None,
    ) -> NDArray[numpy.float128]: ...


class Color:
    # Line colors:
    blue = "#0066ff"
    red = "#ff0000"
    green = "#1ab73a"
    orange = "#ef7810"
    purple = "#8533ff"
    pink = "#d0009a"
    yellow = "#ffcd00"
    # Point colors:
    black = "#000000"
    gray = "#808080ff"
    light_gray = "#bfbfbfff"
    brown = "#916e6e"


LINE_COLOR_CYCLE = [
    Color.blue,
    Color.red,
    Color.green,
    Color.purple,
    Color.pink,
    Color.yellow,
    Color.orange,
]


POINT_COLOR_CYCLE = [
    Color.black,
    Color.gray,
    Color.light_gray,
    Color.brown,
]


class Mark:
    point = "point"
    line = "line"


MARKS = (Mark.point, Mark.line)


def arrhenius(
    ks: ArrayLike,
    T: Sequence[float],  # noqa: N803
    order: int = 1,
    units: UnitsData | None = None,
    labels: Sequence[str] | None = None,
    colors: Sequence[str] | None = None,
    x_label: str = "1000/𝑇",  # noqa: RUF001
    y_label: str = "𝑘",
    mark: str = Mark.line,
    domain: tuple[float, float] | None = None,
) -> altair.Chart:
    """Display as an Arrhenius plot.

    :param others: Other rate constants
    :param others_labels: Labels for other rate constants
    :param T_range: Temperature range
    :param P: Pressure
    :param units: Units
    :param x_label: X-axis label
    :param y_label: Y-axis label
    :param point: Whether to mark with points instead of a line
    :return: Chart
    """
    assert mark in MARKS, f"{mark} not in {MARKS}"
    color_cycle = LINE_COLOR_CYCLE if mark == Mark.line else POINT_COLOR_CYCLE

    nk, nT = numpy.shape(ks)
    colors = colors or list(itertools.islice(itertools.cycle(color_cycle), nk))
    keep_legend = labels is not None
    labels = labels or [f"k{i+1}" for i in range(nk)]
    assert len(T) == nT, f"{T} !~ {ks}"
    assert len(labels) == nk, f"{labels} !~ {ks}"

    # Process units
    units = UNITS if units is None else Units.model_validate(units)
    x_unit = unit_.pretty_string(units.temperature**-1)
    y_unit = unit_.pretty_string(units.rate_constant(order))

    # Add units to labels
    x_label = f"{x_label} ({x_unit})"
    y_label = f"{y_label} ({y_unit})"

    # Gather data from functons
    data_dct = dict(zip(labels, ks, strict=True))
    data = pandas.DataFrame({"x": numpy.divide(1000, T), **data_dct})

    # Determine exponent range
    if domain is None:
        vals_arr = numpy.array(list(data_dct.values()))
        is_nan = numpy.isnan(vals_arr)
        exp_arr = numpy.log10(vals_arr, where=~is_nan)
        exp_arr[is_nan] = 0.0
        exp_arr = numpy.rint(exp_arr).astype(int)
        exp_max = numpy.max(exp_arr).item()
        exp_min = numpy.min(exp_arr).item()
        y_vals = [10**x for x in range(exp_min, exp_max + 2)]
        domain = altair.Undefined
    else:
        y_vals = altair.Undefined

    # Prepare encoding parameters
    x = altair.X("x", title=x_label, scale=altair.Scale(zero=False))
    y = (
        altair.Y("value:Q", title=y_label)
        .scale(type="log", domain=domain)
        .axis(format=".1e")
        .axis(format=".1e", values=y_vals)
    )
    color = (
        altair.Color("key:N", scale=altair.Scale(domain=labels, range=colors))
        if keep_legend
        else altair.value(colors[0])
    )

    chart = altair.Chart(data)
    chart = (
        chart.mark_point(filled=True, opacity=1)
        if mark == Mark.point
        else chart.mark_line()
    )

    # Create chart
    return chart.transform_fold(fold=list(data_dct.keys())).encode(
        x=x, y=y, color=color
    )

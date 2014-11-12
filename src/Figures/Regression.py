__author__ = 'agross'

import pandas as pd
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from Figures.Pandas import series_scatter
from Figures.FigureHelpers import init_ax, prettify_ax
from Helpers.Pandas import match_series


def linear_regression(a, b):
    a, b = match_series(a, b)
    res = sp.stats.linregress(a, b)
    return pd.Series({'slope': res[0], 'intercept': res[1], 'r-value': res[2],
                      'p-value': res[3], 'stderr': res[4]})


def remove_leading_zero(f, places=2):
    f = round(f, places)
    if abs(f - 1) < .01:
        return '1.0'
    elif abs(f) < .01:
        return ''
    elif abs(f) > 1:
        f = str(f)
    elif f > 0:
        f = str(f)[1:]
    else:
        f = '-' + str(f)[2:]
    return f


def regression_string(reg):
    r_value = round(reg['r-value'], 2)
    r_value = str(r_value)[1:] if r_value != 1 else '1.0'
    r_value = '$r={}$'.format(r_value)

    slope = remove_leading_zero(reg['slope'])
    intercept = remove_leading_zero(reg['intercept'])
    if (len(intercept) != 0) and (intercept[0] != '-'):
        intercept = '+' + intercept
    line = '$y={}x {}$'.format(slope, intercept)
    return '\n'.join([r_value, line])


def line_me(slope, intercept, start=0, end=100, ax=None,
            **plot_args):
    if ax is None:
        ax = plt.gca()
    line = lambda x: (slope * x) + intercept
    ax.plot([start, end], [line(start), line(end)],
            **plot_args)


def plot_regression_plain(x, y, ax=None):
    x, y = match_series(x, y)
    fig, ax = init_ax(ax, figsize=(5, 5))
    series_scatter(x, y, ax=ax, ann=None)
    reg = linear_regression(x, y)
    ax.annotate(regression_string(reg), (.5, .05),
                xycoords='axes fraction', size=18)
    line_me(reg['slope'], reg['intercept'], start=x.min(),
            end=x.max(), ax=ax)
    line_me(1, 0, start=x.min(),
            end=x.max(), ax=ax)

    xy = x.append(y)
    ax.set_xbound(xy.min() - 3, xy.max() + 3)
    ax.set_ybound(xy.min() - 3, xy.max() + 3)
    prettify_ax(ax)


def plot_regression_density(x, y, rad=3, ax=None):
    """
    Color density modified from Gordon Bean's Matlab code.
    https://github.com/brazilbean/bean-matlab-toolkit/blob/master/denscat.m
    """
    fig, ax = init_ax(ax, figsize=(5, 5))

    x, y = match_series(x, y)
    pts = pd.concat([x, y], axis=1)
    d = sp.spatial.distance_matrix(pts, pts)
    d = pd.DataFrame(d, pts.index, pts.index)
    area = np.pi * rad ** 2
    dens = 1. * (d < rad).sum() / area
    idx = dens.order().index
    series_scatter(x.ix[idx], y.ix[idx], c=list(dens.ix[idx]), alpha=1, lw=0,
                   ann=None, ax=ax, cmap=cm.jet)

    reg = linear_regression(x, y)
    ax.annotate(regression_string(reg), (.7,.05),
                xycoords='axes fraction', size=18)
    line_me(1, 0, start=x.min(), end=x.max(), ax=ax,
            lw=4, ls='--', color='grey', dash_capstyle='round',
            alpha=.75)

    line_me(reg['slope'], reg['intercept'], start=x.min(),
            end=x.max(), ax=ax, lw=4, ls='-', color='grey',
            dash_capstyle='round', alpha=.75)

    xy = x.append(y)
    ax.set_xbound(xy.min() - 3, xy.max() + 3)
    ax.set_ybound(xy.min() - 3, xy.max() + 3)
    prettify_ax(ax)


def plot_regression(x, y, density=False, rad=3, ax=None):
    if density is True:
        plot_regression_density(x, y, rad, ax)
    else:
        plot_regression_plain(x, y, ax)
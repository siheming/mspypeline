import difflib
from collections.abc import Sized
from difflib import SequenceMatcher
from itertools import combinations


def get_number_rows_cols_for_fig(obj):
    if isinstance(obj, Sized):
        obj = len(obj)
    n_rows, n_cols = 0, 0
    while n_rows * n_cols < obj:
        if n_cols <= n_rows:
            n_cols += 1
        else:
            n_rows += 1
    return n_rows, n_cols


def barplot_annotate_brackets(ax, num1, num2, data, center, height, yerr=None, dh=.05, barh=.05, fs=None,
                              maxasterix=None):
    """
    From: https://stackoverflow.com/questions/11517986/indicating-the-statistically-significant-difference-in-bar-graph
    Annotate barplot with p-values.

    :param ax: axis of plot to put the annotaion brackets
    :param num1: number of left bar to put bracket over
    :param num2: number of right bar to put bracket over
    :param data: string to write or number for generating asterixes
    :param center: centers of all bars (like plt.bar() input)
    :param height: heights of all bars (like plt.bar() input)
    :param yerr: yerrs of all bars (like plt.bar() input)
    :param dh: height offset over bar / bar + yerr in axes coordinates (0 to 1)
    :param barh: bar height in axes coordinates (0 to 1)
    :param fs: font size
    :param maxasterix: maximum number of asterixes to write (for very small p-values)
    """

    if type(data) is str:
        text = data
    else:
        # * is p < 0.05
        # ** is p < 0.005
        # *** is p < 0.0005
        # etc.
        text = ''
        p = .05

        while data < p:
            text += '*'
            p /= 10.

            if maxasterix and len(text) == maxasterix:
                break

        if len(text) == 0:
            text = 'n. s.'

    lx, ly = center[num1], max(height)
    rx, ry = center[num2], max(height)

    if yerr:
        ly += yerr[num1]
        ry += yerr[num2]

    ax_y0, ax_y1 = ax.get_ylim()
    log = False
    if log:
        dh *= 10 ** (ax_y1 - ax_y0)
        barh *= 10 ** (ax_y1 - ax_y0)
    else:
        dh *= (ax_y1 - ax_y0)
        barh *= (ax_y1 - ax_y0)

    y = max(ly, ry) + dh

    barx = [lx, lx, rx, rx]
    bary = [y, y + barh, y + barh, y]
    mid = ((lx + rx) / 2, y + barh)

    ax.plot(barx, bary, c='black')

    kwargs = dict(ha='center', va='bottom')
    if fs is not None:
        kwargs['fontsize'] = fs

    ax.text(*mid, text, **kwargs)


def venn_names(named_sets):
    names = set(named_sets)
    for i in range(1, len(named_sets) + 1):
        for to_intersect in combinations(sorted(named_sets), i):
            others = names.difference(to_intersect)
            intersected = set.intersection(*(named_sets[k] for k in to_intersect))
            unioned = set.union(*(named_sets[k] for k in others)) if others else set()
            yield to_intersect, others, intersected - unioned


def install_r_dependencies(r_package_names, r_bioconducter_package_names):
    from rpy2.robjects.packages import importr
    import rpy2.robjects.packages as rpackages

    r_packages_uninstalled = [x for x in r_package_names if not rpackages.isinstalled(x)]
    r_bioconducter_packages_uninstalled = [x for x in r_bioconducter_package_names if not rpackages.isinstalled(x)]
    if r_packages_uninstalled:
        utils = importr('utils')
        utils.chooseCRANmirror(ind=1)
        for p in r_packages_uninstalled:
            utils.install_packages(p)

    if r_bioconducter_packages_uninstalled:
        biocm = importr("BiocManager")
        for p in r_bioconducter_packages_uninstalled:
            biocm.install(p)


def get_number_of_non_na_values(x):
    import numpy as np
    percentage = 1 / (1 + np.exp(0.5 * x - 3.5)) * 0.5 + 0.5
    return int(np.round(percentage * x))


def string_similarity_ratio(a, b):
    return SequenceMatcher(None, a, b).ratio()


def get_overlap(s1, s2):
    s = difflib.SequenceMatcher(None, s1, s2)
    pos_a, pos_b, size = s.find_longest_match(0, len(s1), 0, len(s2))
    return s1[pos_a:pos_a + size]

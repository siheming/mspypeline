from collections.abc import Sized


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

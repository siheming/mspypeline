from abc import abstractmethod, ABC
from collections import defaultdict as ddict
from typing import Type, Callable, Optional
import pandas as pd
import numpy as np
import logging
import warnings

from mspypeline.helpers import get_logger


def interpolate_data(data: pd.DataFrame) -> pd.DataFrame:
    """
    Performs interpolation on the data by sampling from the same distribution as the input distribution.
    Adopted from https://github.com/bmbolstad/preprocessCore,
    more specifically: https://github.com/bmbolstad/preprocessCore/blob/master/src/qnorm.c
    Parameters
    ----------
    data
        A DataFrame with columns being the experiments and rows and being the features

    Returns
    -------
    A DataFrame where all values have been replaced by interpolating from the old data column wise.
    For the non missing entries the new values are very close to the old,
    while the for the missing entries a sampled value is assigned
    """
    data_arg_sort = np.argsort(data.values, axis=0)
    data_sorted = np.take_along_axis(data.values, data_arg_sort, axis=0)
    data_index = data.index[data_arg_sort]

    rows, cols = data_sorted.shape

    num_na = data.isna().sum(axis=0)
    not_na = rows - num_na

    # create a linspace for each column, each with #rows entries that span from 0 to #non-missing values in that column
    float_index = (np.tile(np.linspace(0, 1, rows), (cols, 1)) * (not_na.values[:, np.newaxis] - 1)).T
    index = np.floor(float_index).astype(int)
    decimal = float_index % 1

    index_min = np.minimum(index + 1, not_na.values[np.newaxis, :] - 1)
    # take the weighted average of the two neighboring values, based on the decimal
    # slightly simplified expression (1 - dec) * sorted_data[index] + decimal * sorted_data[index + 1]
    new_values = (1 - decimal) * np.take_along_axis(data_sorted, index, axis=0) +\
                 decimal * np.take_along_axis(data_sorted, index_min, axis=0)

    # now the index of the new values needs to be reconstructed, since each column was sorted differently
    result = []
    # reconstruct one column at a time
    for column_index, column_name in enumerate(data.columns):
        series_index = np.empty((rows,))
        series_index[:] = np.nan

        index_mapping = ddict(lambda: (-1, np.inf))
        for i, (index_val, float_val) in enumerate(zip(index[:, column_index], float_index[:, column_index])):
            decimal = float_val - index_val
            if index_mapping[index_val][1] > decimal:
                index_mapping[index_val] = i, decimal

            if index_mapping[index_val + 1][1] > 1 - decimal:
                index_mapping[index_val + 1] = i, 1 - decimal

        # the index that was added last needs to be removed
        index_mapping.pop(index_val + 1)

        for key, (idx, _) in index_mapping.items():
            series_index[idx] = key

        series_index[np.isnan(series_index)] = np.arange(np.nanmax(series_index) + 1, rows, 1)
        s_ind = data_index[:, column_index][series_index.astype(int)]
        result.append(pd.Series(new_values[:, column_index], index=s_ind, name=column_name))

    return pd.concat(result, axis=1, sort=False)


def median_polish(data: pd.DataFrame, max_iter: int = 100, tol: float = 0.001):
    overall = np.nanmedian(data.values)
    row_effect = pd.Series([0] * data.shape[0], index=data.index)
    column_effect = pd.Series([0] * data.shape[1], index=data.columns)
    residuals = data - overall

    for i in range(max_iter):
        # row collapse
        row_medians = residuals.median(axis=1)
        overall += column_effect.median()
        row_effect += row_medians
        residuals = residuals.sub(row_medians, axis=0)
        column_effect -= column_effect.median()
        # column collapse
        column_medians = residuals.median(axis=0)
        overall += row_effect.median()
        column_effect += column_medians
        residuals = residuals.sub(column_medians, axis=1)
        row_effect -= row_effect.median()
        # check stop condition
        # could be adapted to check difference between last two updates
        if column_medians.abs().sum() + row_medians.abs().sum() <= tol:
            break
    else:
        warnings.warn("stopping because max iter was reached")
    return {"ave": overall, "row_effect": row_effect, "col_effect": column_effect, "residual": residuals}


class BaseNormalizer(ABC):
    def __init__(self, input_scale: str = "log2",
                 output_scale: str = "normal",
                 col_name_prefix: Optional[str] = None,
                 loglevel: int = logging.DEBUG,
                 **kwargs):
        """
        Abstract base class for Normalizers

        Parameters
        ----------
        input_scale
            Scale of the input data. Either normal or log2
        output_scale
            Scale of the output data. Either normal or log2
        col_name_prefix
            If not None the prefix is added to each column name
        loglevel
            loglevel of the logger
        kwargs
        """
        self.logger = get_logger(self.__class__.__name__, loglevel)
        allowed_scales = ("log2", "normal")
        if input_scale not in allowed_scales:
            raise ValueError("input_scale should be one of: " + ", ".join(allowed_scales))
        if output_scale not in allowed_scales:
            raise ValueError("out_scale should be one of: " + ", ".join(allowed_scales))
        self.input_scale = input_scale
        self.output_scale = output_scale
        self.col_name_prefix = col_name_prefix

    @abstractmethod
    def fit(self, data: pd.DataFrame):
        raise NotImplementedError

    @abstractmethod
    def transform(self, data: pd.DataFrame):
        raise NotImplementedError

    def fit_transform(self, data: pd.DataFrame):
        return self.fit(data).transform(data)


class MedianNormalizer(BaseNormalizer):
    def __init__(self, input_scale: str = "log2",
                 output_scale: str = "normal",
                 col_name_prefix: Optional[str] = None,
                 loglevel: int = logging.DEBUG,
                 **kwargs):
        super().__init__(input_scale, output_scale, col_name_prefix, loglevel, **kwargs)
        self.factors = None

    def fit(self, data: pd.DataFrame):
        if self.input_scale == "normal":
            data = np.log2(data)
        medians = data.median()
        self.factors = medians - medians.mean()
        return self

    def transform(self, data: pd.DataFrame):
        if self.factors is None:
            raise ValueError("Please call fit first or use fit_transform")
        if self.input_scale == "normal":
            data = np.log2(data)
        result = data - self.factors
        if self.output_scale == "normal":
            result = np.exp2(result)
        if self.col_name_prefix is not None:
            result.rename(lambda x: f"{self.col_name_prefix} {x}", axis=1, inplace=True)
        return result


class QuantileNormalizer(BaseNormalizer):
    """
    Quantile Normalizer as described on wikipedia
    https://en.wikipedia.org/wiki/Quantile_normalization
    """
    def __init__(self, missing_value_handler: Optional[Callable] = interpolate_data,
                 input_scale: str = "log2",
                 output_scale: str = "normal",
                 col_name_prefix: Optional[str] = None,
                 loglevel: int = logging.DEBUG,
                 **kwargs):
        super().__init__(input_scale, output_scale, col_name_prefix, loglevel, **kwargs)
        self.missing_value_handler = missing_value_handler
        self.rank_replace = {}

    def fit(self, data: pd.DataFrame):
        if self.input_scale == "normal":
            data = np.log2(data)
        if self.missing_value_handler is not None:
            data = self.missing_value_handler(data)
        sorted_values = pd.DataFrame(np.sort(data.values, axis=0))
        sorted_mean = sorted_values.mean(axis=1)
        self.rank_replace.update({rank + 1.: intensity for rank, intensity in sorted_mean.iteritems()})
        # also include all half ranks
        half_steps = pd.concat(
            (sorted_mean, sorted_mean[1:].append(pd.Series([np.nan]), ignore_index=True)),
            axis=1).mean(axis=1)
        self.rank_replace.update({rank + 1.5: intensity for rank, intensity in half_steps.iteritems()})
        return self

    def transform(self, data: pd.DataFrame):
        if not self.rank_replace:
            raise ValueError("Please call fit first or use fit_transform")
        if self.missing_value_handler is not None:
            na_mask = ~data.isna()
            data = self.missing_value_handler(data)
        result = data.rank()
        if self.missing_value_handler is not None:
            result = result[na_mask]
        result = result.apply(pd.Series.map, arg=self.rank_replace)
        if self.output_scale == "normal":
            result = np.exp2(result)
        if self.col_name_prefix is not None:
            result.rename(lambda x: f"{self.col_name_prefix} {x}", axis=1, inplace=True)
        return result


class TailRobustNormalizer(BaseNormalizer):
    """
    https://www.biorxiv.org/content/10.1101/2020.04.17.046227v1.full
    """
    def __init__(self, normalizer: Type[BaseNormalizer] = QuantileNormalizer,
                 missing_value_handler: Optional[Callable] = interpolate_data,
                 input_scale: str = "log2",
                 output_scale: str = "normal",
                 col_name_prefix: Optional[str] = None,
                 loglevel: int = logging.DEBUG,
                 **kwargs):
        super().__init__(input_scale, output_scale, col_name_prefix, loglevel, **kwargs)
        self.offset_factor = None
        self.normalizer = normalizer
        self.missing_value_handler = missing_value_handler

    def fit(self, data: pd.DataFrame):
        if self.input_scale == "normal":
            data = np.log2(data)
        self.offset_factor = data.mean(axis=1)
        return self

    def transform(self, data: pd.DataFrame):
        if self.offset_factor is None:
            raise ValueError("Please call fit first or use fit_transform")
        if self.input_scale == "normal":
            data = np.log2(data)
        result = data.subtract(self.offset_factor, axis=0)
        result = self.normalizer(missing_value_handler=self.missing_value_handler, input_scale="log2",
                                 output_scale="log2", col_name_prefix=None, loglevel=self.logger.getEffectiveLevel()
                                 ).fit_transform(result)
        result = result.add(self.offset_factor, axis=0)
        if self.output_scale == "normal":
            result = np.exp2(result)
        if self.col_name_prefix is not None:
            result.rename(lambda x: f"{self.col_name_prefix} {x}", axis=1, inplace=True)
        return result


default_normalizers = {
    "median_norm": MedianNormalizer(),
    "quantile_norm": QuantileNormalizer(missing_value_handler=None),
    "quantile_norm_missing_handled": QuantileNormalizer(),
    "trqn": TailRobustNormalizer(missing_value_handler=None),
    "trqn_missing_handled": TailRobustNormalizer(),
    "trmn": TailRobustNormalizer(normalizer=MedianNormalizer, missing_value_handler=None),
}
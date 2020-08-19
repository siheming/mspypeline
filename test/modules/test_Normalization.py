import pandas as pd
import numpy as np
import pytest
from sklearn.exceptions import ConvergenceWarning


def test_interpolate_data():
    from mspypeline.modules.Normalization import interpolate_data
    assert interpolate_data(pd.DataFrame()).equals(pd.DataFrame())
    data = pd.DataFrame(np.random.random((100, 100)))
    data[np.random.random((100, 100)) > 0.5] = np.nan
    assert interpolate_data(data).isna().sum().sum() == 0


def test_median_polish():
    from mspypeline.modules.Normalization import median_polish
    with pytest.warns(RuntimeWarning) as record:
        median_polish(pd.DataFrame())
        # check that only one warning was raised
        assert len(record) == 1
        # check that the message matches
        assert record[0].message.args[0] == "Mean of empty slice"
    with pytest.warns(ConvergenceWarning) as record:
        median_polish(pd.DataFrame(np.random.random((10, 10))), max_iter=1)
        assert len(record) == 1
    # TODO testcase with known data and result


def test_base_normalizer():
    from mspypeline.modules.Normalization import BaseNormalizer

    class NormTest(BaseNormalizer):
        def fit(self, data):
            super().fit(data)

        def transform(self, data):
            super().transform(data)

    nt = NormTest()
    with pytest.raises(NotImplementedError):
        nt.fit(pd.DataFrame())
    with pytest.raises(NotImplementedError):
        nt.transform(pd.DataFrame())

    with pytest.raises(ValueError):
        NormTest(input_scale="test")
    with pytest.raises(ValueError):
        NormTest(output_scale="test")


def test_default_normalizers():
    from mspypeline.modules.Normalization import default_normalizers
    data = pd.DataFrame(np.random.random((100, 100)))
    data_copy = data.copy()
    for norm_name, norm in default_normalizers.copy().items():
        setattr(norm, "input_scale", "normal")
        setattr(norm, "col_name_prefix", norm_name)
        with pytest.raises(ValueError):
            norm.transform(data)
        norm.fit_transform(data)
        assert data.equals(data_copy)


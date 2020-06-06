def test_logger():
    import pytest
    from mspypeline.helpers import get_logger
    from logging import Logger
    assert isinstance(get_logger("test", "warning"), Logger)
    assert isinstance(get_logger(loglevel=10), Logger)

    with pytest.warns(RuntimeWarning) as record:
        assert isinstance(get_logger(loglevel="asd"), Logger)
        assert len(record) == 1

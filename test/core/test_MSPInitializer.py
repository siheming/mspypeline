import os
import shutil


def test_initializer():
    from mspypeline import MSPInitializer
    script_loc = os.path.dirname(os.path.realpath(__file__))
    dir_name = "initializer_test"
    dir_path = os.path.join(script_loc, dir_name)
    shutil.rmtree(dir_path, ignore_errors=True)
    os.makedirs(os.path.join(dir_path, "txt"), exist_ok=True)
    ini = MSPInitializer(dir_path)
    ini.start_dir = os.path.join(dir_path, "txt")
    assert os.path.split(ini.start_dir)[1] != "txt"
    ini.init_config()
    ini.update_config_file()
    ini.configs["test"] = True
    ini.update_config_file()
    # TODO set some pathways
    ini.init_interest_from_txt()
    del ini
    ini = MSPInitializer(dir_path)
    assert ini.configs["test"] is True
    del ini
    with open(os.path.join(dir_path, "config", MSPInitializer.yml_file_name), "a") as f:
        f.write("\n")
        f.write("test1: true")
    ini = MSPInitializer(dir_path)
    assert ini.configs["test1"] is True
    shutil.rmtree(dir_path)

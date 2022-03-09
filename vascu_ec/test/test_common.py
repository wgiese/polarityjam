import logging
import os
import tempfile
import unittest
from pyunpack import Archive
from pathlib import Path

import yaml

from vascu_ec.logging import LOGGER_NAME
from vascu_ec.utils.io import list_files_recursively


class TestCommon(unittest.TestCase):
    """Base class all tests inherit from. Responsible for temporary directory, cleanup, parameters."""
    def setUp(self) -> None:
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.current_path = Path(os.path.dirname(os.path.realpath(__file__)))
        self.parameters, self.param_base_file = self.load_parameters()

    def tearDown(self) -> None:
        self.tmp_dir.cleanup()
        logging.getLogger(LOGGER_NAME).handlers.clear()

    def extract_test_data(self):
        Archive(str(self.current_path.joinpath("resources", "data.zip"))).extractall(str(self.tmp_dir.name))

    def get_test_image_path(self, image_name):
        data_path = Path(self.tmp_dir.name).joinpath("data")
        if not data_path.exists():
            self.extract_test_data()

        files_list = list_files_recursively(data_path)
        test_image_path = files_list[[i.stem + i.suffix for i in files_list].index(image_name)]

        return test_image_path

    def get_test_image_folder(self, mode):
        # "gn" = golgi_nuclei, "g" = "no golgi", "n" = "no nuclei"

        data_path = Path(self.tmp_dir.name).joinpath("data")
        if not data_path.exists():
            self.extract_test_data()

        if mode == "gn":
            return data_path.joinpath("golgi_nuclei")
        elif mode == "n":
            return data_path.joinpath("no_nuclei")
        elif mode == "g":
            return data_path.joinpath("no_golgi")
        else:
            return None

    def get_test_parameter_file(self, name):
        return self.current_path.joinpath("resources", name)

    def get_test_key_file(self):
        return self.current_path.joinpath("resources", "test_key_file.csv")

    def load_parameters(self):
        param_base_file = Path(self.current_path).joinpath("..", "utils", "base", "parameters.yml")

        with open(param_base_file, 'r') as yml_f:
            parameters = yaml.safe_load(yml_f)

        parameters["channel_junction"] = 3
        parameters["channel_nucleus"] = 2
        parameters["channel_golgi"] = 0

        return parameters, param_base_file

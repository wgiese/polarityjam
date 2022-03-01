import os
import tempfile
import unittest
from pathlib import Path

import yaml


class TestCommon(unittest.TestCase):

    def setUp(self) -> None:
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.current_path = Path(os.path.dirname(os.path.realpath(__file__)))
        self.parameters, self.param_base_file = self.load_parameters()

    def tearDown(self) -> None:
        self.tmp_dir.cleanup()

    def get_test_image_path(self, image_name):
        test_image_path = self.current_path.joinpath("resources", image_name)

        # 0 golgi
        # 1 not needed
        # 2 nuclei
        # 3 junctions

        return test_image_path

    def load_parameters(self):
        param_base_file = Path(self.current_path).joinpath("..", "utils", "base", "parameters.yml")

        with open(param_base_file, 'r') as yml_f:
            parameters = yaml.safe_load(yml_f)

        parameters["channel_junction"] = 3
        parameters["channel_nucleus"] = 2
        parameters["channel_golgi"] = 0

        return parameters, param_base_file

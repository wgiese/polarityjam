from pathlib import Path

import numpy as np
import yaml

from vascu_ec.feature_extraction import get_image_for_segmentation
from vascu_ec.test.test_common import TestCommon
from vascu_ec.utils.io import read_parameters, read_image


class TestFunctions(TestCommon):
    """Small tests to ensure basic IO functionality of the tool."""

    def setUp(self) -> None:
        super().setUp()

    def test_read_parameters(self):
        # prepare
        local_parameter_file = {
            "channel_junction": 0,
            "channel_nucleus": 1,
            "channel_golgi": 2
        }
        local_parameter_file_path = Path(self.tmp_dir.name).joinpath("local_parameter_file.yml")

        with open(local_parameter_file_path, "w+") as file:
            file.write(yaml.dump(local_parameter_file, Dumper=yaml.Dumper))

        self.parameters["channel_junction"] = 0
        self.parameters["channel_nucleus"] = 1
        self.parameters["channel_golgi"] = 2

        # call
        r = read_parameters(local_parameter_file_path)

        # assert
        self.assertDictEqual(self.parameters, r)

    def test_read_img(self):
        # prepare
        test_image_path = self.get_test_image_path("060721_EGM2_18dyn_01.tif")

        # call
        r = read_image(test_image_path)

        # assert
        self.assertIsNotNone(r)

        # assert correct type
        self.assertEqual(type(np.array(0)), type(r))

    def test_get_image_for_segmentation_case_junctions(self):
        # prepare
        img = read_image(self.get_test_image_path("060721_EGM2_18dyn_01.tif"))

        # call
        r = get_image_for_segmentation(self.parameters, img)

        # assert
        self.assertEqual((2, 1024, 1024), r.shape)

    def test_get_image_for_segmentation_case_segmentation(self):
        self.parameters["channel_nucleus"] = -1
        img = read_image(self.get_test_image_path("060721_EGM2_18dyn_01.tif"))

        # call
        r = get_image_for_segmentation(self.parameters, img)

        # expect junctions channel only
        self.assertEqual((1024, 1024), r.shape)


class EmptyTestClass:
    pass

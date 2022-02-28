import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch, MagicMock

import numpy as np
import yaml
from vascu_ec.feature_extraction import get_image_for_segmentation,     get_features_from_cellpose_seg_multi_channel
from vascu_ec.utils.io import read_parameters, read_image
from vascu_ec.utils.seg import get_cellpose_segmentation


class TestFunctions(unittest.TestCase):

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
        param_base_file = Path(self.current_path).joinpath("..", "base", "parameters.yml")

        with open(param_base_file, 'r') as yml_f:
            parameters = yaml.safe_load(yml_f)

        parameters["channel_junction"] = 3
        parameters["channel_nucleus"] = 2
        parameters["channel_golgi"] = 0

        return parameters, param_base_file

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
        test_image_path = self.get_test_image_path("multichannel.tif")

        # call
        r = read_image(test_image_path)

        # assert
        self.assertIsNotNone(r)

        # assert correct type
        self.assertEqual(type(np.array(0)), type(r))

        # assert shape is correct
        # todo

    def test_get_image_for_segmentation_case_junctions(self):
        # prepare
        img = read_image(self.get_test_image_path("multichannel.tif"))

        # call
        r = get_image_for_segmentation(self.parameters, img)

        # assert
        self.assertEqual((2, 1024, 1024), r.shape)

    def test_get_image_for_segmentation_case_segmentation(self):
        self.parameters["channel_nucleus"] = -1
        img = read_image(self.get_test_image_path("multichannel.tif"))

        # call
        r = get_image_for_segmentation(self.parameters, img)

        # expect junctions channel only
        self.assertEqual((1024, 1024), r.shape)

    @unittest.skip
    @patch('vascu_ec.functions.skimage.segmentation.clear_border')
    @patch("vascu_ec.functions.get_cellpose_model")
    def test_get_cellpose_segmentation(self, clear_border_mock, get_cellpose_model_mock):
        # mock
        eval_mock = MagicMock(return_value=("mask", "flows", "styles", "diams"))

        fake_model = EmptyTestClass()
        fake_model.eval = eval_mock

        get_cellpose_model_mock.return_value = fake_model

        # call
        r = get_cellpose_segmentation(self.parameters, "myImgSeg")

        # assert
        self.assertEqual("mask", r)
        eval_mock.assert_called_once_with("myImgSeg", diameter=100, channels=[1, 2])
        clear_border_mock.assert_called_once()

    def test_get_features_from_cellpose_seg(self):
        img = read_image(self.get_test_image_path("multichannel.tif"))
        mask = np.load('/home/jpa/PycharmProjects/vascu_ec-app/vascu_ec/test/resources/multichannel_seg.npy',
                       allow_pickle=True)

        get_features_from_cellpose_seg_multi_channel(
            self.parameters,
            img,
            mask,
            "bla",
            '/home/jpa/PycharmProjects/vascu_ec-app/vascu_ec/test/resources/'
        )


class EmptyTestClass:
    pass

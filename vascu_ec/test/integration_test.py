from pathlib import Path

import sys

from vascu_ec.argument_parsing import startup
from vascu_ec.test.test_common import TestCommon
from vascu_ec.utils.io import write_dict_to_yml


class TestIntegration(TestCommon):

    def setUp(self) -> None:
        super().setUp()
        self.output_path = Path(self.tmp_dir.name).joinpath("output")

    def test_run(self):
        # prepare
        param_file = str(Path(self.tmp_dir.name).joinpath("param_file_test_run.yml"))
        out = str(self.output_path)
        in_file = str(self.get_test_image_path("multichannel.tif"))

        write_dict_to_yml(param_file, self.parameters)

        # build arguments
        sys.argv = sys.argv + ['run', param_file, out, in_file, '--filename_prefix=myfile']

        # call
        startup()

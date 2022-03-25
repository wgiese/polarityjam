import sys
from pathlib import Path

from vascu_ec.argument_parsing import startup
from vascu_ec.test.test_common import TestCommon


class TestIntegration(TestCommon):

    def setUp(self) -> None:
        super().setUp()
        self.output_path = Path(self.tmp_dir.name).joinpath("output")

    def test_run(self):
        # prepare
        in_file = str(self.get_test_image_path("060721_EGM2_18dyn_01.tif"))
        param_file = str(self.get_test_parameter_file("parameters_golgi_nuclei.yml"))
        out_path = str(self.output_path)

        # build arguments
        sys.argv = [sys.argv[0]] + ['run', param_file, out_path, in_file, '--filename_prefix=myfile']

        # call
        startup()

        # todo: compare results!

    def test_run_stack(self):
        in_path = str(self.get_test_image_folder("gn").joinpath("set_2"))
        param_file = str(self.get_test_parameter_file("parameters_golgi_nuclei.yml"))
        out_path = str(self.output_path)

        # build arguments
        sys.argv = [sys.argv[0]] + ['run-stack', param_file, out_path, in_path]

        # call
        startup()

    def test_run_stack_no_golgi(self):
        in_path = str(self.get_test_image_folder("g"))
        param_file = str(self.get_test_parameter_file("parameters_no_golgi.yml"))
        out_path = str(self.output_path)

        # build arguments
        sys.argv = [sys.argv[0]] + ['run-stack', param_file, out_path, in_path]

        # call
        startup()

    def test_run_stack_no_nuclei(self):
        in_path = str(self.get_test_image_folder("n"))
        param_file = str(self.get_test_parameter_file("parameters_no_nuclei.yml"))
        out_path = str(self.output_path)

        # build arguments
        sys.argv = [sys.argv[0]] + ['run-stack', param_file, out_path, in_path]

        # call
        startup()

    def test_run_key(self):
        in_path = str(self.get_test_image_folder("gn"))
        in_key = str(self.get_test_key_file())
        param_file = str(self.get_test_parameter_file("parameters_golgi_nuclei.yml"))
        out_path = str(self.output_path)

        # build arguments
        sys.argv = [sys.argv[0]] + ['run-key', param_file, out_path, in_path, in_key]

        # call
        startup()

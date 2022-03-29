import time
import unittest

from vascu_ec.vascu_ec_logging import get_logger
from vascu_ec.test import test_integration, test_feature_extraction


def start_tests():
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    # unittests
    suite.addTests(loader.loadTestsFromModule(test_feature_extraction))

    # integration tests
    suite.addTests(loader.loadTestsFromModule(test_integration))

    runner = unittest.TextTestRunner(verbosity=3)
    result = runner.run(suite)
    if result.wasSuccessful():
        time.sleep(5)
        get_logger().info("Success")
        exit(0)
    else:
        get_logger().warning("Failed")
        exit(1)


if __name__ == "__main__":
    start_tests()

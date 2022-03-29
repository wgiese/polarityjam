import argparse
import sys

import vascu_ec
from vascu_ec.commandline import run, run_stack, run_key
from vascu_ec.vascu_ec_logging import configure_logger, get_logger


def startup():
    """Entry points of `vascu-ec`."""
    configure_logger('INFO')

    parser = create_parser()
    args = parser.parse_args()

    __run_subcommand(args, parser)


def __run_subcommand(args, parser):
    """Calls a specific album subcommand."""
    command = ""
    try:
        command = sys.argv[1]  # command always expected at second position
    except IndexError:
        parser.error("Please provide a valid action!")
    get_logger().debug("Running %s subcommand..." % command)

    get_logger().info("VascuShare Version %s" %
                      (vascu_ec.__version__))

    args.func(args)  # execute entry point function


def create_parser():
    parser = VascuParser()

    # run action
    p = parser.create_file_command_parser('run', run, 'Feature extraction from a single tiff image.')
    p.add_argument('in_file', type=str, help='Path to the input tif file.')
    p.add_argument('--filename_prefix', type=str, help='prefix for the output file.', required=False, default=None)

    # stack action
    p = parser.create_file_command_parser(
        'run-stack', run_stack, 'Feature extraction from an input folder containing several tiff files.'
    )
    p.add_argument('in_path', type=str, help='Name for the input folder containing tif images.')

    # key action
    p = parser.create_file_command_parser(
        'run-key',
        run_key,
        'Feature extraction from a given list of folders provided within a csv file.'
    )
    p.add_argument('in_path', type=str, help='Path prefix. Note: Base path for all folders listed in the key file.')
    p.add_argument(
        'in_key', type=str, help='Path to the input key file. File must contain column '
                                 '"folder_name" and "condition_name". Note: Folder name relative to the path prefix.'
    )

    return parser.parser


class ArgumentParser(argparse.ArgumentParser):
    """Override default error method of all parsers to show help of sub-command."""

    def error(self, message):
        self.print_help()
        self.exit(2, '%s: error: %s\n' % (self.prog, message))


class VascuParser(ArgumentParser):
    def __init__(self):
        super().__init__()
        self.parent_parser = self.create_parent_parser()
        self.parser = self.create_parser()
        self.subparsers = self.parser.add_subparsers(title='actions', help='sub-command help')

    @staticmethod
    def create_parent_parser():
        """Parent parser for all subparsers to have the same set of arguments."""
        parent_parser = ArgumentParser(add_help=False)
        # parse logging
        parent_parser.add_argument(
            '--log',
            required=False,
            help='Logging level for your command. Choose between %s' %
                 ", ".join(['INFO', 'DEBUG']),
            default='INFO'
        )
        return parent_parser

    def create_parser(self):
        """Creates the main parser for the framework."""
        parser = ArgumentParser(
            add_help=True,
            description='VascuShare - feature extraction pipeline for multichannel images.',
            parents=[self.parent_parser])
        return parser

    def create_command_parser(self, command_name, command_function, command_help):
        """Creates a parser for a VascuShare command, specified by a name, a function and a help description."""
        parser = self.subparsers.add_parser(command_name, help=command_help, parents=[self.parent_parser])
        parser.set_defaults(func=command_function)
        return parser

    def create_file_command_parser(self, command_name, command_function, command_help):
        """Creates a parser for a VascuShare command dealing with a file.

        Parser is specified by a name, a function and a help description.
        """
        parser = self.create_command_parser(command_name, command_function, command_help)
        parser.add_argument('param', type=str, help='Path to the parameter file.')
        parser.add_argument('out_path', type=str, help='Path to the output folder.')
        return parser

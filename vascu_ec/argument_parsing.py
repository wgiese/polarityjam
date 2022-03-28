import argparse
import sys

import vascu_ec
from vascu_ec.commandline import run, run_stack, run_key
from vascu_ec.logging import configure_logger, get_logger

# todo: nice help messages


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

    get_logger().info("VascuShare Version %s | contact via %s " %
                      (vascu_ec.__version__, vascu_ec.__email__))

    args.func(args)  # execute entry point function


def create_parser():
    parser = VascuParser()

    # run action
    p = parser.create_file_command_parser('run', run, 'help')
    p.add_argument('in_file', type=str, help='path to the input tif file.')
    p.add_argument('--filename_prefix', type=str, help='prefix for the output file.', required=False, default=None)

    # stack action
    p = parser.create_file_command_parser('run-stack', run_stack, 'help')
    p.add_argument('in_path', type=str, help='name for the input folder containing tifs.')

    # key action
    p = parser.create_file_command_parser('run-key', run_key, 'help')
    p.add_argument('in_path', type=str, help='name for the input folder containing tifs.')
    p.add_argument('in_key', type=str, help='path to the key file.')

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
        parser.add_argument('param', type=str, help='path to the parameter file.')
        parser.add_argument('out_path', type=str, help='output path.')
        return parser

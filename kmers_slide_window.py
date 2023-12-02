#!/usr/bin/env python
# -*- coding: utf-8 -*-


def parse_arguments():
    """Parse the command line arguments to the script xxxx.
    Sets up the argparse command-line parser and calls it. These args can be accessed
    using args.args.
    """
    parser = argparse.ArgumentParser(
        description=""" """,
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-',
                        '--',
                        metavar='',
                        type=str,
                        required=,
                        dest='',
                        help='')

    return parser.parse_args()


def main():
    # starting count the staring time of the script
    start = time()
    # checking the current directory and printing it
    cwd = os.getcwd()
    print(colored(f'\nThe working directory: {cwd}\n',
                  'green',
                  attrs=['bold']))
    # passing the arguments to the script
    args = parse_arguments()

    end = time()
    # print some info
    print(colored(f'Total time for the script finishes: {round(end - start, 2)}.',
                  'red',
                  attrs=['bold']))
    print(colored('Done!',
                  'green',
                  attrs=['bold']))


if __name__ == "__main__":
    sys.exit(main())

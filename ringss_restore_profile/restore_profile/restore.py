import sys

from ringss_restore_profile.core import read_json
from ringss_restore_profile.core import read_stm


def run_restore_profile():
    """Main entrypoint"""

    print("Hello from entrypoint")
    if len(sys.argv) < 3:
        print("Usage: python readstm.py <par-file.json> <file.stm>")
        sys.exit()

    parameters_filename = sys.argv[1]
    stm_filename = sys.argv[2]
    print('STM file: ' + stm_filename)

    # Ingest all information needed
    parameters = read_json(parameters_filename)
    if parameters is None:
        sys.exit()

    read_stm(parameters, stm_filename)


if __name__ == '__main__':
    run_restore_profile()

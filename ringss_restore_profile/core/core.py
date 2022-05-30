import json


def read_json(filename: str):
    """Utility function to read json files

    Args:
        filename (str): Name of json file to open

    Returns:
        a dictionary; None if fails

    """
    try:
        with open(filename, "r") as fp:
            data = json.load(fp=fp)
            return data
    except FileNotFoundError as err:
        print(err)
        # log.error(str(error)
        return None

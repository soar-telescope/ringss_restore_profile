import json


def read_json(filename):
    """Utility function to read json files

    Args:
        filename (str): Name of json file to open

    Returns:
        a dictionary; None if fails

    """
    try:
        file = open(filename, "r")
        p = json.load(file)
        file.close()
        return p
    except FileNotFoundError as err:
        print(err)
        return None

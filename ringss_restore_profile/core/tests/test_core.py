from unittest import TestCase

from ..core import read_json


class TestReadJsonFile(TestCase):

    def test_read_json__file_does_not_exist(self):
        no_file = 'file_that_does_not_exists.json'
        data = read_json(filename=no_file)
        self.assertIsNone(data)


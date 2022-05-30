import json
import os

from unittest import TestCase

from ..core import read_json


class TestReadJsonFile(TestCase):

    file_name = 'test_json_file.json'

    def setUp(self) -> None:
        self.data = {'name': 'RINGSS',
                     'type': 'Telescope'}

        with open(self.file_name, 'w') as fp:
            json.dump(self.data, fp)

    def tearDown(self) -> None:
        if os.path.exists(self.file_name):
            os.unlink(self.file_name)

    def test_read_json__file_does_not_exist(self):
        no_file = 'file_that_does_not_exists.json'
        data = read_json(filename=no_file)
        self.assertIsNone(data)

    def test_read_json__file_exists(self):
        data = read_json(filename=self.file_name)
        self.assertDictEqual(self.data, data)


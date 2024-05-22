# pylint: disable=all
##### currently using notebook for testing, this is a package level testing framework
import unittest
import sys
import os

# relative and absolute imports for running module and script respectively
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
try:
    from ..configs.config1 import Config
    from ..utils.helper import *
except ImportError:
    from configs.config1 import Config
    from ..utils.helper import *

class TestHelper(unittest.TestCase):
    def setUp(self):
        self.config = Config("GM12878", 100000)
    
import unittest
from tieler.tile import powers2


class TestTile(unittest.TestCase):

    def test_powers(self):
        check = all(sum(2**p for p in powers2(n)) == n for n in (9, 13, 425, 123))

        self.assertTrue(check)

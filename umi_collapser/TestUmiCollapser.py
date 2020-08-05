import unittest
import collapser_functions

class TestUmiCollapser(unittest.TestCase):
    def test_posterior_calculation(self):
        [base_call, score, cigar] = collapser_functions.call_base_posterior(('A', 'A', 'G', 'T'), (40, 40, 5, 5))
        print(base_call, ord(score) - 33)
        self.assertGreater(10, 0)


if __name__ == '__main__':
    unittest.main()
import unittest
import collapser_functions


class TestUmiCollapser(unittest.TestCase):
    def test_quality_increase_with_coverage(self):
        # Check that increasing coverage increases quality
        qual = 20
        [base_call, score, cigar] = collapser_functions.call_base_posterior(('A','A'), (qual, qual))
        new_qual = ord(score) - 33
        self.assertEqual(base_call, 'A')
        self.assertGreater(new_qual, qual)

    def test_quality_decrease_with_mismatch(self):
        qual = 20
        [base_call, score, cigar] = collapser_functions.call_base_posterior(('A', 'A'), (qual, qual))
        [base_call_2, score_2, cigar_2] = collapser_functions.call_base_posterior(('A', 'A', 'T'), (qual, qual, qual))
        new_qual = ord(score) - 33
        new_qual_2 = ord(score_2) - 33
        self.assertGreater(new_qual ,new_qual_2)

    def test_best_base_wins(self):
        [base_call, score, cigar] = collapser_functions.call_base_posterior(('A', 'G',  'C','T'), (10, 10, 10, 30))
        self.assertEqual(base_call, 'T')

    def test_tie(self):
        [base_call, score, cigar] = collapser_functions.call_base_posterior(('A', 'A', 'T', 'T'), (10, 10, 10, 10))
        self.assertEqual(base_call, 'N')

    def test_break_tie(self):
        [base_call, score, cigar] = collapser_functions.call_base_posterior(('A', 'A', 'T', 'T'), (11, 10, 10, 10))
        self.assertEqual(base_call, 'A')


if __name__ == '__main__':
    unittest.main()
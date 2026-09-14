"""Regression tests for the source-binding acceptance boundary."""
import tempfile
import unittest
from pathlib import Path
from hashlib import sha256
import random
from Bio.Seq import Seq
import rebuild_assignment_phylogenies as binding


class SourceSequenceBindingTests(unittest.TestCase):
    def alignment(self, query, source):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'source.fa'
            path.write_text('>source\n'+source+'\n')
            return binding.source_alignment(query, str(path), sha256(path.read_bytes()).hexdigest())

    def setUp(self):
        rng = random.Random(73419)
        self.sequence = ''.join(rng.choice('ACGT') for _ in range(3000))

    def test_identical_complete_sequence(self):
        hits = self.alignment(self.sequence, self.sequence)
        self.assertEqual(sum(hit['exact'] for hit in hits), 1)

    def test_reverse_complement(self):
        hits = self.alignment(self.sequence, str(Seq(self.sequence).reverse_complement()))
        self.assertEqual([hit['strand'] for hit in hits if hit['exact']], ['-'])

    def test_omitted_source_insertion_is_counted(self):
        source = self.sequence[:1500]+'ATGCCGGTAGCACTG'+self.sequence[1500:]
        hits = self.alignment(self.sequence, source)
        self.assertEqual([(hit['omitted_source_bases'], hit['mismatching_bases']) for hit in hits if hit['exact']], [(15, 0)])

    def test_substitution_is_not_an_exact_source_binding(self):
        base = 'A' if self.sequence[1500] != 'A' else 'C'
        source = self.sequence[:1500]+base+self.sequence[1501:]
        hits = self.alignment(self.sequence, source)
        self.assertFalse(any(hit['exact'] for hit in hits))
        self.assertTrue(any(hit['mismatching_bases'] == 1 for hit in hits))

    def test_query_insertion_is_not_silently_omitted(self):
        query = self.sequence[:1500]+'ATGCCGGTAGCACTG'+self.sequence[1500:]
        hits = self.alignment(query, self.sequence)
        self.assertFalse(any(hit['exact'] for hit in hits))
        self.assertTrue(any(hit['query_inserted_bases'] == 15 for hit in hits))

    def test_arbitrary_subsequence_is_not_a_source_alignment(self):
        source = ''.join('CGT'+base for base in self.sequence)
        hits = self.alignment(self.sequence, source)
        self.assertFalse(any(hit['exact'] for hit in hits))


if __name__ == '__main__':
    unittest.main()

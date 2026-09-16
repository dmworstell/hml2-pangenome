"""The same minimum product length applies before and after late frameshifts."""
import unittest
from Bio.Seq import Seq
from orf_analysis import analyze_orf_structural_integrity


class ShortORFLengthRule(unittest.TestCase):
    def call(self, sample, reference, feature='gag', locus_type='TypeII'):
        return analyze_orf_structural_integrity(
            sample, str(Seq(reference).translate()), reference,
            feature_name=feature, locus_type=locus_type,
            gapped_ref_dna=reference,
        )[0]

    def test_short_late_frameshift(self):
        reference = 'ATG' + 'AAA' * 98 + 'TAA'
        sample = reference[:170] + '-' * (len(reference) - 170)
        self.assertEqual(self.call(sample, reference), 'Fragment_Intact')

    def test_long_late_frameshift_preserved(self):
        reference = 'ATG' + 'AAA' * 98 + 'TAA'
        sample = reference[:-1] + '-'
        self.assertEqual(self.call(sample, reference), 'Intact_FS_End')

    def test_intact_and_type1_env_unchanged(self):
        reference = 'ATG' + 'AAA' * 98 + 'TAA'
        for feature, locus_type in [('gag', 'TypeII'), ('env', 'TypeI'), ('rec', 'TypeII')]:
            self.assertEqual(self.call(reference, reference, feature, locus_type), 'Intact')

    def test_early_frameshift_not_promoted(self):
        reference = 'ATG' + 'AAA' * 98 + 'TAA'
        sample = reference[:89] + '-' * (len(reference) - 89)
        self.assertEqual(self.call(sample, reference), 'Deletion')


if __name__ == '__main__':
    unittest.main()

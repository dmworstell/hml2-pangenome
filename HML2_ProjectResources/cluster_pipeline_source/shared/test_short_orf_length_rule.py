"""The same minimum product length applies before and after late frameshifts."""
import unittest
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from orf_analysis import analyze_orf_structural_integrity, REC_EXONS_TYPE2, translate_and_check


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

    def test_rec_reference_includes_complete_stop(self):
        reference_file = Path(__file__).resolve().parents[2] / 'data/ref/type2_KCON.fa'
        sequence = str(SeqIO.read(reference_file, 'fasta').seq)
        reference = ''.join(sequence[a:b] for a, b in REC_EXONS_TYPE2)
        protein, status = translate_and_check(reference)
        self.assertEqual(len(reference), 318)
        self.assertEqual(len(protein.rstrip('*')), 105)
        self.assertTrue(protein.endswith('*'))
        self.assertEqual(status, 'intact')
        result = analyze_orf_structural_integrity(reference, protein, reference,
            feature_name='rec', locus_type='TypeII', gapped_ref_dna=reference)
        self.assertEqual(result[0], 'Intact')
        self.assertEqual(result[3], [])
        self.assertEqual(self.call(reference[:-1] + '-', reference, 'rec'), 'Intact_FS_End')

    def test_rec_length_threshold_retains_existing_stop_convention(self):
        reference = 'ATG' + 'AAA' * 104 + 'TGA'
        # Existing caller counts the reference stop in its 60% denominator.
        sample = reference[:189] + '-' * (len(reference) - 189)
        self.assertEqual(self.call(sample, reference, 'rec'), 'Fragment_Intact')


if __name__ == '__main__':
    unittest.main()

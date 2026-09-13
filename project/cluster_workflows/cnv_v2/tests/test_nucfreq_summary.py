import unittest

from nucfreq_summary import count_region_sites


def site(position, counts):
    return '\t'.join(map(str, ['chrTest', position, position + 1, *counts])) + '\n'


class OriginalSiteRuleTests(unittest.TestCase):
    def test_both_thresholds_are_inclusive(self):
        result = count_region_sites([site(100, [17, 3, 0, 0])], 100, 100)
        self.assertEqual(result['nucfreq_het_sites'], 1)

    def test_secondary_count_below_three_is_not_mixed(self):
        result = count_region_sites([site(100, [2, 2, 0, 0])], 100, 100)
        self.assertEqual(result['nucfreq_total_sites'], 1)
        self.assertEqual(result['nucfreq_het_sites'], 0)

    def test_fraction_below_point15_is_not_mixed(self):
        result = count_region_sites([site(100, [18, 3, 0, 0])], 100, 100)
        self.assertEqual(result['nucfreq_het_sites'], 0)

    def test_denominator_is_top_two_not_all_bases(self):
        result = count_region_sites([site(100, [17, 3, 3, 0])], 100, 100)
        self.assertEqual(result['nucfreq_het_sites'], 1)

    def test_counts_are_ranked_not_reference_based(self):
        result = count_region_sites([site(100, [0, 3, 0, 17])], 100, 100)
        self.assertEqual(result['nucfreq_het_sites'], 1)

    def test_any_positive_top_two_depth_enters_denominator(self):
        result = count_region_sites([site(100, [1, 0, 0, 0])], 100, 100)
        self.assertEqual(result['nucfreq_total_sites'], 1)
        self.assertEqual(result['nucfreq_het_frac'], 0)

    def test_zero_depth_is_not_in_denominator(self):
        result = count_region_sites([site(100, [0, 0, 0, 0])], 100, 100)
        self.assertEqual(result['nucfreq_total_sites'], 0)
        self.assertIsNone(result['nucfreq_het_frac'])

    def test_historical_numeric_bounds_are_inclusive_and_unshifted(self):
        result = count_region_sites([site(p, [3, 3, 0, 0]) for p in range(99, 104)], 100, 102)
        self.assertEqual(result['nucfreq_total_sites'], 3)
        self.assertEqual(result['nucfreq_het_sites'], 3)

    def test_bare_ranked_count_variant(self):
        result = count_region_sites(['chrTest\t100\t17\t3\n'], 100, 100)
        self.assertEqual(result['nucfreq_het_sites'], 1)

    def test_mean_secondary_depth_uses_mixed_sites_only(self):
        result = count_region_sites([site(100, [17, 3, 0, 0]), site(101, [85, 15, 0, 0]), site(102, [18, 2, 0, 0])], 100, 102)
        self.assertEqual(result['nucfreq_total_sites'], 3)
        self.assertEqual(result['nucfreq_mean_alt_depth'], 9)

    def test_one_base_shift_changes_only_boundary_coordinate_sets(self):
        for start in range(1, 20):
            for length in range(1, 20):
                end = start + length - 1
                original = set(range(start, end + 1))
                shifted = set(range(start - 1, end))
                self.assertEqual(original - shifted, {end})
                self.assertEqual(shifted - original, {start - 1})


if __name__ == '__main__':
    unittest.main()

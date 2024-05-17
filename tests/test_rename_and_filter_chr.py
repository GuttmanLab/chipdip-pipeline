import collections
import sys
sys.path.append('../scripts/python')
import rename_and_filter_chr

########################
# Test rename_and_filter_chr.reheader()
########################

old_header = collections.OrderedDict(
    [('HD', {'SO': 'coordinate', 'VN': '1.5'}),
     ('SQ',
      [{'LN': 248956422, 'SN': 'chr1'},
       {'LN': 242193529, 'SN': 'chr2'},
       {'LN': 198295559, 'SN': 'chr3'}]),
    ])

# order of new chromosomes is consistent with existing
chrom_map1 = {
    'chr1': 'chr6',
    'chr2': 'chr7',
    'chr3': 'chr8'
}
old_to_new_refID1 = {0: 0, 1: 1, 2: 2}
retains_sorting1 = True
new_header1 = collections.OrderedDict(
    [('HD', {'SO': 'coordinate', 'VN': '1.5'}),
     ('SQ',
      [{'LN': 248956422, 'SN': 'chr6'},
       {'LN': 242193529, 'SN': 'chr7'},
       {'LN': 198295559, 'SN': 'chr8'}]),
    ])

new_header, old_to_new_refID, retains_sorting = rename_and_filter_chr.reheader(old_header, chrom_map1)
assert new_header == new_header1
assert old_to_new_refID == old_to_new_refID1
assert retains_sorting == retains_sorting1

# order of new chromosomes is different than existing --> reads will need to be resorted
chrom_map2 = {
    'chr2': 'chr7',
    'chr1': 'chr6',
    'chr3': 'chr8'
}
old_to_new_refID2 = {0: 1, 1: 0, 2: 2}
retains_sorting2 = False
new_header2 = collections.OrderedDict(
    [('HD', {'SO': 'unsorted', 'VN': '1.5'}),
     ('SQ',
      [{'LN': 242193529, 'SN': 'chr7'},
       {'LN': 248956422, 'SN': 'chr6'},
       {'LN': 198295559, 'SN': 'chr8'}]),
    ])

new_header, old_to_new_refID, retains_sorting = rename_and_filter_chr.reheader(old_header, chrom_map2)
assert new_header == new_header2
assert old_to_new_refID == old_to_new_refID2
assert retains_sorting == retains_sorting2

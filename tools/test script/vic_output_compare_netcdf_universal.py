#!/usr/bin/env python

""" This script checks any two VIC NetCDF output files against each other, even if one is
 time-major and the other is cell-major (just set the command line parameters accordingly).
  
 For usage information: ./vic_output_compare_netcdf_universal.py --help """
 
import sys
import argparse
import os.path
import h5py
import numpy as np
from collections import defaultdict

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s]n' % message)
        self.print_help()
        sys.exit(2)

parser = MyParser()
parser.add_argument('--testfile', action="store", dest="testfile", type=str, help = 'file name and path of the NetCDF file to be tested')
parser.add_argument('-test-time-major', action="store_true", dest="test_is_time_major", default=True, help = 'file to be tested is time-major (default=True)')
parser.add_argument('-test-cell-major', action="store_false", dest="test_is_time_major", default=True, help = 'file to be tested is cell-major (default=False)')
parser.add_argument('-test-depth-position', action="store", dest="pos_depth_dim_test", type=int, default=1, help = 'set the 0-based position of the depth dimension of 4D variables in the NetCDF file to be tested (default=1, assuming time-major test file option [-test-time-major] is set, i.e. <time, depth, lat, lon>)')
parser.add_argument('--basefile', action="store", dest="basefile", type=str, help = 'file name and path of the base NetCDF file to check against')
parser.add_argument('-base-time-major', action="store_true", dest="base_is_time_major", default=True, help = 'base file is time-major (default=True)')
parser.add_argument('-base-cell-major', action="store_false", dest="base_is_time_major", default=True, help = 'base file is cell-major (default=False)')
parser.add_argument('-base-depth-position', action="store", dest="pos_depth_dim_base", type=int, default=1, help = 'set the 0-based position of the depth dimension of 4D variables in the base NetCDF file (default=1, assuming time-major base file option [-base-time-major] is set, i.e. <time, depth, lat, lon>)')
parser.add_argument('-test-start-rec', action="store", dest="test_start_rec", type=int, default=0, help = 'choose which time record in testfile to start comparison from (default=0)')
parser.add_argument('-test-end-rec', action="store", dest="test_end_rec", type=int, default=0, help = 'choose which time record in testfile to end comparison at (default is last time record)')
parser.add_argument('-base-start-rec', action="store", dest="base_start_rec", type=int, default=0, help = 'choose which time record in basefile to start comparison from (default=0)')
parser.add_argument('-base-end-rec', action="store", dest="base_end_rec", type=int, default=0, help = 'choose which time record in basefile to end comparison at (default is last time record)')
parser.add_argument('-tolerance', action="store", dest="tolerance", type=float, default=0, help = 'set the absolute tolerance that values must be within to be considered in agreement')
parser.add_argument('--csv', action="store_true", dest="csv_out", default=False, help = 'set if you want output of data read from NetCDF to CSV')
parser.add_argument('-csv-diffs-only', action="store_true", dest="csv_diffs_only", default=False, help = 'set if you want to output only the differences to CSV and not the original data')
parser.add_argument('--v', action="store_true", dest="verbose", default=False, help = 'for verbose output')

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

options = parser.parse_args()
testfile = options.testfile
basefile = options.basefile
test_is_time_major = options.test_is_time_major
base_is_time_major = options.base_is_time_major
pos_depth_dim_test = options.pos_depth_dim_test
pos_depth_dim_base = options.pos_depth_dim_base
csv_out = options.csv_out
csv_diffs_only = options.csv_diffs_only
tolerance = options.tolerance
test_start_rec = options.test_start_rec
test_end_rec = options.test_end_rec
base_start_rec = options.base_start_rec
base_end_rec = options.base_end_rec
verbose = options.verbose

print '\n------- vic_output_compare_netcdf_universal -------'
# print interpretation of input
print 'Checking testfile {} against basefile {}'.format(testfile, basefile)

if test_is_time_major:
  print 'testfile declared as time-major with depth dimension in position {}'.format(pos_depth_dim_test)
else:
  print 'testfile declared as cell-major with depth dimension in position {}'.format(pos_depth_dim_test)

if base_is_time_major:
  print 'basefile declared as time-major with depth dimension in position {}'.format(pos_depth_dim_base)
else:
  print 'basefile declared as cell-major with depth dimension in position {}'.format(pos_depth_dim_base)

print 'Absolute tolerance (-tol) for agreement between elements: {}'.format(tolerance)

if csv_out:
  print 'CSV output selected (--csv).'
  if csv_diffs_only:
    print 'Only differences between input files will be saved (-csv_diffs_only).'

# Open test output NetCDF file
#testH5 = h5py.File('../out/results.nc', 'r')
testH5 = h5py.File(testfile, 'r')

# Open out_base NetCDF file
#baseH5 = h5py.File('../out_base/results.nc')
baseH5 = h5py.File(basefile, 'r')

# grab the number of time records
test_time_len = len(testH5['time'])
base_time_len = len(baseH5['time'])

if not test_end_rec:
    test_end_rec = test_time_len
    
if not base_end_rec:
    base_end_rec = base_time_len
    
num_test_recs = test_end_rec - test_start_rec
num_base_recs = base_end_rec - base_start_rec 

print 'Checking testfile time record range {}:{} ({} records total)'.format(test_start_rec, test_end_rec, num_test_recs)
print 'Checking basefile time record range {}:{} ({} records total)'.format(base_start_rec, base_end_rec, num_base_recs)
if num_test_recs != num_base_recs:
    print 'Number of time records selected for comparison between testfile and basefile must be equal.  Exiting.\n'
    sys.exit(0)

# grab NaN fill value of non-initialized NetCDF records from one common attribute (precipitation)
try: # the original forcing file precip variable
    fill_value = baseH5['pr'].attrs['_FillValue']
except: # the new disaggregated forcing / model run output precip variable
    fill_value = baseH5['OUT_PREC'].attrs['_FillValue']

# grab lat and lon dimensions of grid cells
lats = testH5['lat'][:]
lons = testH5['lon'][:]

lat_to_idx = dict([(x[1],x[0]) for x in enumerate(lats.tolist())])
lon_to_idx = dict([(x[1],x[0]) for x in enumerate(lons.tolist())])

# initialize with keys from NetCDF file, and empty values
cell_data_keys = dict.fromkeys(testH5.keys()) 
del cell_data_keys['lat']
del cell_data_keys['lon']
del cell_data_keys['time']
del cell_data_keys['bnds']
del cell_data_keys['depth']
#del cell_data_keys['OUT_WIND'] # uncomment this to test new vs. old disaggregated forcing files

# need this in order to nest defaultdict objects beyond 2 levels
def tree(): return defaultdict(tree)

# create a big nested dictionary with all data for all cells, one for test data and one for out_base data
all_test_data = tree()
all_base_data = tree()

# load up all_test_data
for lat in lats:
    for lon in lons:
        cell_label = '{}_{}'.format(lat, lon)
        print '\nLoading cell {} data...'.format(cell_label)
        for variable in cell_data_keys:
            if verbose:
                print 'loading variable {}'.format(variable)
            if test_is_time_major == True: # test file is time-major format
                if len(testH5[variable].shape) == 4: # 4D variable
                    if pos_depth_dim_test == 1: # <time, depth, lat, lon>
                        for depth in range(0, testH5[variable].shape[pos_depth_dim_test]):
                            all_test_data[cell_label][variable][depth] = testH5[variable][test_start_rec:test_end_rec,depth,lat_to_idx[lat],lon_to_idx[lon]]
                    elif pos_depth_dim_test == 3: # <time, lat, lon, depth>
                        for depth in range(0, testH5[variable].shape[pos_depth_dim_test]):
                        #print 'cell_label: {} variable: {} depth: {}'.format(cell_label, variable, depth)
                            all_test_data[cell_label][variable][depth] = testH5[variable][test_start_rec:test_end_rec,lat_to_idx[lat],lon_to_idx[lon],depth]
                    else:
                        print 'The declared depth dimension position {} for a time-major NetCDF test file is not supported.'.format(pos_depth_dim_test)
                        sys.exit(0)
                elif len(testH5[variable].shape) == 3: # 3D variable
                    all_test_data[cell_label][variable] = testH5[variable][test_start_rec:test_end_rec,lat_to_idx[lat],lon_to_idx[lon]]

            elif test_is_time_major == False: # test file is cell-major format
                if len(testH5[variable].shape) == 4: # 4D variable
                    if pos_depth_dim_test == 0: # <depth, lat, lon, time>
                        for depth in range(0, testH5[variable].shape[pos_depth_dim_test]):
                            all_test_data[cell_label][variable][depth] = testH5[variable][depth,lat_to_idx[lat],lon_to_idx[lon],test_start_rec:test_end_rec]
                    else:
                        print 'The declared depth dimension position {} for a cell-major NetCDF test file is not supported.'.format(pos_depth_dim_test)
                        sys.exit(0)
                elif len(testH5[variable].shape) == 3: # 3D variable
                    all_test_data[cell_label][variable] = testH5[variable][lat_to_idx[lat],lon_to_idx[lon],test_start_rec:test_end_rec]
            # overwrite all testfile NaNs (fill values) with 0          
            if len(testH5[variable].shape) == 3: # 3D variable
                all_test_data[cell_label][variable][all_test_data[cell_label][variable] == fill_value] = 0
            elif len(testH5[variable].shape) == 4: # 4D variable
                for depth in range(0, testH5[variable].shape[pos_depth_dim_test]):
                    all_test_data[cell_label][variable][depth][all_test_data[cell_label][variable][depth] == fill_value] = 0

            if base_is_time_major == True: # base file is time-major format
                if len(baseH5[variable].shape) == 4: # 4D variable
                    if pos_depth_dim_base == 1: # <time, depth, lat, lon>
                        for depth in range(0, baseH5[variable].shape[pos_depth_dim_base]):
                            all_base_data[cell_label][variable][depth] = baseH5[variable][base_start_rec:base_end_rec,depth,lat_to_idx[lat],lon_to_idx[lon]]
                    elif pos_depth_dim_base == 3: # <time, lat, lon, depth>
                        for depth in range(0, baseH5[variable].shape[pos_depth_dim_base]):
                            all_base_data[cell_label][variable][depth] = baseH5[variable][base_start_rec:base_end_rec,lat_to_idx[lat],lon_to_idx[lon],depth]
                    else:
                        print 'The declared depth dimension position {} for a time-major NetCDF base file is not supported.'.format(pos_depth_dim_base)
                        sys.exit(0)
                elif len(baseH5[variable].shape) == 3: # 3D variable
                    all_base_data[cell_label][variable] = baseH5[variable][base_start_rec:base_end_rec,lat_to_idx[lat],lon_to_idx[lon]]

            elif base_is_time_major == False: # base file is cell-major format
                if len(baseH5[variable].shape) == 4: # 4D variable
                    if pos_depth_dim_base == 0: # <depth, lat, lon, time>
                        for depth in range(0, baseH5[variable].shape[pos_depth_dim_base]):
                            all_base_data[cell_label][variable][depth] = baseH5[variable][depth,lat_to_idx[lat],lon_to_idx[lon],base_start_rec:base_end_rec]
                    else:
                        print 'The declared depth dimension position {} for a cell-major NetCDF base file is not supported.'.format(pos_depth_dim_base)
                        sys.exit(0)
                elif len(baseH5[variable].shape) == 3: # 3D variable
                    all_base_data[cell_label][variable] = baseH5[variable][lat_to_idx[lat],lon_to_idx[lon],base_start_rec:base_end_rec]
                              
            # overwrite all basefile NaNs (fill values) with 0          
            if len(baseH5[variable].shape) == 3: # 3D variable
                all_base_data[cell_label][variable][all_base_data[cell_label][variable] == fill_value] = 0
            elif len(baseH5[variable].shape) == 4: # 4D variable
                for depth in range(0, baseH5[variable].shape[pos_depth_dim_base]):
                    all_base_data[cell_label][variable][depth][all_base_data[cell_label][variable][depth] == fill_value] = 0


## compare test data against base data
if verbose:
    print'Checking agreement between test output data with base values:'
# create labels for each cell to use to index into defaultdicts we created
cell_labels = []
for lat_label in lats:
    for lon_label in lons:
        cell_labels.append(repr(lat_label) + '_' + repr(lon_label))
for cell in cell_labels:
    diffs_exist = False
    # if we want to output the test dataset to CSV format for inspection
    if csv_out == True: 
        # Create tables that will be written into CSV files
        # time record index read from file will form the first column of the output tables,
        # relative time record index (0 at test/base_start_rec position) in second column
        general_headers = ['time_rec','relative_time_rec']
        test_table = testH5['time'][test_start_rec:test_end_rec]
        test_table = np.column_stack([test_table, range(0,num_test_recs)])
        base_table = baseH5['time'][base_start_rec:base_end_rec]
        base_table = np.column_stack([base_table, range(0,num_base_recs)])
        band_headers = ['time_rec', 'relative_time_rec'] 
        test_band_table = testH5['time'][test_start_rec:test_end_rec]
        test_band_table = np.column_stack([test_band_table, range(0,num_test_recs)])
        base_band_table = baseH5['time'][base_start_rec:base_end_rec]
        base_band_table = np.column_stack([base_band_table, range(0,num_base_recs)])
        diffs_headers = ['time', 'relative_time_rec']
        diffs_table = baseH5['time'][base_start_rec:base_end_rec]
        diffs_table = np.column_stack([diffs_table, range(0,num_base_recs)])
        diffs_band_headers = ['time', 'relative_time_rec']
        diffs_band_table = baseH5['time'][base_start_rec:base_end_rec]
        diffs_band_table = np.column_stack([diffs_band_table, range(0,num_base_recs)])

    if verbose:    
        print'\n'
        print'Cell ' + str(cell) + ': \n'
        print ''
    for variable in cell_data_keys:
        diffs = []
        diffs_depths = []
        if len(testH5[variable].shape) == 3: # 3D variable
            if tolerance > 0:
                agreement = np.allclose(all_test_data[cell][variable], all_base_data[cell][variable], 0, tolerance)
            else:
                agreement = np.array_equal(all_test_data[cell][variable], all_base_data[cell][variable])
            if csv_out == True:
                column_header = str(testH5[variable].attrs['internal_vic_name'])
                general_headers.append(column_header)
                test_table = np.column_stack([test_table, all_test_data[cell][variable]])
                base_table = np.column_stack([base_table, all_base_data[cell][variable]])
            if agreement == False:
                diffs_exist = True
                diffs = abs(all_test_data[cell][variable] - all_base_data[cell][variable])
                if csv_out == True:
                    diffs_table = np.column_stack([diffs_table, diffs])        
                    diffs_headers.append(column_header) 
                num_diffs = len(diffs[diffs > tolerance])
                max_diff = np.max(diffs)
                sum_diffs = np.sum(diffs)
                if verbose:
                    print '    ' + str(testH5[variable].attrs['internal_vic_name']) + ': ' + str(agreement),
                    print '(Number of different entries: {} '.format(num_diffs),
                    print 'Maximum absolute difference: {} '.format(max_diff),
                    print 'Total sum of differences: {})'.format(sum_diffs) 
            else:
                if verbose:
                    print '    ' + str(testH5[variable].attrs['internal_vic_name']) + ': ' + str(agreement)
        elif len(testH5[variable].shape) == 4: # 4D variable
            max_diff = 0
            num_diffs = 0
            for depth in range(0, testH5[variable].shape[pos_depth_dim_test]):
                if tolerance > 0:
                    agreement = np.allclose(all_test_data[cell][variable][depth], all_base_data[cell][variable][depth], 0, tolerance)
                else:
                    agreement = np.array_equal(all_test_data[cell][variable][depth], all_base_data[cell][variable][depth])
                if csv_out == True:
                    test_band_table = np.column_stack([test_band_table, all_test_data[cell][variable][depth]])
                    base_band_table = np.column_stack([base_band_table, all_base_data[cell][variable][depth]])
                    column_header = str(testH5[variable].attrs['internal_vic_name']) + '_' + str(depth)
                    band_headers.append(column_header)
                if agreement == False:
                    diffs_exist = True
                    diffs_band = abs(all_test_data[cell][variable][depth] - all_base_data[cell][variable][depth])
                    if csv_out == True:
                        #diffs_band_table = np.column_stack([diffs_band_table, diffs_band[depth]])        
                        diffs_band_table = np.column_stack([diffs_band_table, diffs_band])        
                        diffs_band_headers.append(column_header)
                    diffs_depths.append(depth)
                    max_temp = np.max(diffs_band) 
                    max_diff = max_temp if max_temp > max_diff else max_diff # running max of differences across all bands
                    num_diffs += len(diffs_band[diffs_band > tolerance]) # running number of differences across all bands
            if verbose:    
                if not diffs_depths: # no differences were found at any depth
                    print '    ' + str(testH5[variable].attrs['internal_vic_name']) + ': ' + str(agreement)
                else:
                    print '    ' + str(testH5[variable].attrs['internal_vic_name']) + ': False ' + '(Differences at depths ' + str(diffs_depths) + ' ',
                    print 'Number of different entries across all bands: {} '.format(num_diffs), 
                    print 'Maximum absolute difference across all bands: {})'.format(max_diff) 
    
if csv_out == True:
    if csv_diffs_only == False:
        test_csv_filename = 'tabular_cell_{}_{}.csv'.format(cell, os.path.basename(testfile))
        np.savetxt(test_csv_filename, test_table, delimiter=',', fmt='%3.22f', header=",".join(general_headers))
        base_csv_filename = 'tabular_cell_{}_{}_base.csv'.format(cell, os.path.basename(basefile))
        np.savetxt(base_csv_filename, base_table, delimiter=',', fmt='%3.22f', header=",".join(general_headers))

        test_4D_csv_filename = 'tabular_cell_{}_{}_band.csv'.format(cell, os.path.basename(testfile)) 
        np.savetxt(test_4D_csv_filename, test_band_table, delimiter=',', fmt='%3.22f', header=",".join(band_headers))
        base_4D_csv_filename = 'tabular_cell_{}_{}_band_base.csv'.format(cell, os.path.basename(basefile))
        np.savetxt(base_4D_csv_filename, base_band_table, delimiter=',', fmt='%3.22f', header=",".join(band_headers))

    if diffs_exist:
        diffs_3D_csv_filename = 'tabular_cell_{}_{}_differences_tol={}.csv'.format(cell, os.path.basename(testfile), tolerance)
        np.savetxt(diffs_3D_csv_filename, diffs_table, delimiter=',', fmt='%3.22f', header=",".join(diffs_headers))     
        diffs_4D_csv_filename = 'tabular_cell_{}_{}_band_differences_tol={}.csv'.format(cell, os.path.basename(testfile), tolerance)
        np.savetxt(diffs_4D_csv_filename, diffs_band_table, delimiter=',', fmt='%3.22f', header=",".join(diffs_band_headers))

if diffs_exist:
    print '\nDifferences exist between testfile and basefile at the given tolerance of {}\n'.format(tolerance)
else:
    print '\ntestfile and basefile are in agreement within the given tolerance of {}.\n'.format(tolerance)

print '\nvic_output_compare_netcdf_universal finished.'

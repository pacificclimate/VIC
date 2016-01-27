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

debug = False

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
parser.add_argument('--var-map-file', action="store", dest="var_map_file", type=str, help = 'file name and path of a file providing NetCDF output variable name mapping (e.g. Snow Melt could be named snm / SNOW_MELT / OUT_SNOW_MELT, depending on the version of VIC or user-defined mappings).  Each line of the file should start with the general OUT_ variable name, followed by all possible alternate names, tab-separated, e.g. OUT_SNOW_MELT   SNOW_MELT   snm')

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
var_map_file = options.var_map_file

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
    print 'basefile declared as cell-major with depth dimension in position {}\n'.format(pos_depth_dim_base)

print 'Absolute tolerance (-tolerance) for agreement between elements: {}\n'.format(tolerance)

if csv_out:
    print 'CSV output selected (--csv).'
    if csv_diffs_only:
        print 'Only differences between input files will be saved (-csv_diffs_only).'
    # print '\n'

# Open test output NetCDF file
#testH5 = h5py.File('../out/results.nc', 'r')
testH5 = h5py.File(testfile, 'r')

# Open out_base NetCDF file
#baseH5 = h5py.File('../out_base/results.nc')
baseH5 = h5py.File(basefile, 'r')

# grab NaN fill value of non-initialized NetCDF records
for var in baseH5:
    try:
        fill_value_base = baseH5[var].attrs['_FillValue']
    except:
        continue
    if fill_value_base is not None:
        break
print 'basefile NetCDF _FillValue: {}'.format(fill_value_base)
for var in testH5:
    try:
        fill_value_test = testH5[var].attrs['_FillValue']
    except:
        continue
    if fill_value_test is not None:
        break
print 'testfile NetCDF _FillValue: {}'.format(fill_value_test)
if csv_out:
    print 'Empty records in output CSV files will be assigned the basefile NetCDF _FillValue.'

# grab the number of time records
test_time_len = len(testH5['time'])
base_time_len = len(baseH5['time'])

if not test_end_rec:
    test_end_rec = test_time_len
    
if not base_end_rec:
    base_end_rec = base_time_len
    
num_test_recs = test_end_rec - test_start_rec
num_base_recs = base_end_rec - base_start_rec 

print '\ntestfile data time record range to be read: {}:{} ({} records total)'.format(test_start_rec, test_end_rec, num_test_recs)
print 'basefile data time record range to be read: {}:{} ({} records total)'.format(base_start_rec, base_end_rec, num_base_recs)
if num_test_recs != num_base_recs:
    print 'Number of time records selected for comparison between testfile and basefile must be equal.  Exiting.\n'
    sys.exit(0)

# grab lat and lon dimensions of grid cells
lats = testH5['lat'][:]
lons = testH5['lon'][:]

lat_to_idx = dict([(x[1],x[0]) for x in enumerate(lats.tolist())])
lon_to_idx = dict([(x[1],x[0]) for x in enumerate(lons.tolist())])

# initialize with keys from NetCDF file, and empty values
test_data_keys = dict.fromkeys(testH5.keys()) 
del test_data_keys['lat']
del test_data_keys['lon']
del test_data_keys['time']
del test_data_keys['bnds']
del test_data_keys['depth']
#del test_data_keys['OUT_WIND'] # uncomment this to test new vs. old disaggregated forcing files

# initialize with keys from NetCDF file, and empty values
base_data_keys = dict.fromkeys(baseH5.keys()) 
del base_data_keys['lat']
del base_data_keys['lon']
del base_data_keys['time']
del base_data_keys['bnds']
del base_data_keys['depth']

# define basefile to testfile variable name mapping (e.g. precipitation might be called OUT_PREC in one, and PREC or pr in another)
var_map = {}
if var_map_file:
    print '\nOutput variable name mapping defined in file {} will be used.'.format(var_map_file)
    with open(var_map_file, 'r') as f:
        for line in f:
            if (len(line) > 1):
                split_line = line.split()
                var = split_line[0]
                var_map[var] = []
                for idx, mapped_var_name in enumerate(split_line):
                    if idx > 0:
                        var_map[var].append(mapped_var_name)
                        #print 'appended {} to var_map[{}]'.format(mapped_var_name, var)
else: # use the variable names provided in the testfile, assuming the basefile variable naming agrees
    #for var in baseH5:
    for var in test_data_keys:
        var_map[var] = var
print 'var_map: {}'.format(var_map)
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
        for variable in test_data_keys:
            if debug:
                print 'loading variable from test file: {}'.format(variable)
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
            if fill_value_test != fill_value_base:
                # overwrite all testfile NetCDF _FillValue (NaNs) entries read in to all_test_data with the basefile _FillValue          
                if len(testH5[variable].shape) == 3: # 3D variable
                    all_test_data[cell_label][variable][all_test_data[cell_label][variable] == fill_value_test] = fill_value_base
                elif len(testH5[variable].shape) == 4: # 4D variable
                    for depth in range(0, testH5[variable].shape[pos_depth_dim_test]):
                        all_test_data[cell_label][variable][depth][all_test_data[cell_label][variable][depth] == fill_value_test] = fill_value_base
        for variable in base_data_keys:
            if debug:
                print 'loading variable from base file: {}'.format(variable)
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
                              

## compare test data against base data
print'\nChecking agreement between testfile variables data and those in the basefile...\n (Note: this will only check the set \
of variables existing in your testfile against those same ones in the basefile, which may have a larger set of variables than testfile)'

# global flag indicating whether there are *any* differences between the files
diffs_exist = False
# create labels for each cell to use to index into defaultdicts we created
cell_labels = []
for lat_label in lats:
    for lon_label in lons:
        cell_labels.append(repr(lat_label) + '_' + repr(lon_label))
for cell in cell_labels:
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
       # print 'var_map: {}'.format(var_map)
    # Have to handle the situation where the output variable names for the same thing may differ between base and test files
    for variable in var_map:
        test_var_found = False
        base_var_found = False
        agreement = False
        diffs = []
        diffs_depths = []

        if verbose:
            print 'checking {}...'.format(variable)
        # Determine the name of the variable in the testfile
        if variable in test_data_keys:
            test_variable = variable
            test_var_found = True
        else: # look for other matches in var_map loaded from file
            for var_name in var_map[variable]:
                if var_name in test_data_keys:
                    test_variable = var_name
                    test_var_found = True  
        if test_var_found:  
            if debug:
                    print 'Found match: {} = {} in testfile'.format(variable, test_variable)
        else: # TODO: add option to skip this variable if not present (and tell the user)
            print 'No entry for {} output variable found in testfile.  Exiting.\n'.format(variable)
            sys.exit(0)
        # Determine the name of the variable in the basefile
        if variable in base_data_keys:
            base_variable = variable
            base_var_found = True
        else: # look for other matches in var_map loaded from file
            for var_name in var_map[variable]:
                if var_name in base_data_keys:
                    base_variable = var_name
                    base_var_found = True
        if base_var_found:  
            if debug:
                    print 'Found match: {} = {} in basefile'.format(variable, base_variable)
        else: # TODO: add option to skip this variable if not present (and tell the user)
            print 'No entry for {} output variable found in basefile.  Exiting.\n'.format(variable)
            sys.exit(0)

        if len(testH5[test_variable].shape) == 3: # 3D variable
            if tolerance > 0:
                agreement = np.allclose(all_test_data[cell][test_variable], all_base_data[cell][base_variable], 0, tolerance)
            else:
                agreement = np.array_equal(all_test_data[cell][test_variable], all_base_data[cell][base_variable])
            if csv_out == True:
                column_header = variable
                general_headers.append(column_header)
                test_table = np.column_stack([test_table, all_test_data[cell][test_variable]])
                base_table = np.column_stack([base_table, all_base_data[cell][base_variable]])
            if agreement == False:
                diffs_exist = True
                diffs = abs(all_test_data[cell][test_variable] - all_base_data[cell][base_variable])
                if csv_out == True:
                    diffs_table = np.column_stack([diffs_table, diffs])        
                    diffs_headers.append(column_header) 
                num_diffs = len(diffs[diffs > tolerance])
                max_diff = np.max(diffs)
                sum_diffs = np.sum(diffs)
                if verbose:
                    print '    ' + variable + ': ' + str(agreement)
                    print '      Number of different entries: {} '.format(num_diffs),
                    print 'Maximum absolute difference: {} '.format(max_diff),
                    print 'Total sum of differences: {} '.format(sum_diffs) 
            else:
                if verbose:
                    print '    ' + variable + ': ' + str(agreement)
        elif len(testH5[test_variable].shape) == 4: # 4D variable
            max_diff = 0
            num_diffs = 0
            for depth in range(0, testH5[test_variable].shape[pos_depth_dim_test]):
                if tolerance > 0:
                    agreement = np.allclose(all_test_data[cell][test_variable][depth], all_base_data[cell][base_variable][depth], 0, tolerance)
                else:
                    agreement = np.array_equal(all_test_data[cell][test_variable][depth], all_base_data[cell][base_variable][depth])
                if csv_out == True:
                    test_band_table = np.column_stack([test_band_table, all_test_data[cell][test_variable][depth]])
                    base_band_table = np.column_stack([base_band_table, all_base_data[cell][base_variable][depth]])
                    column_header = variable + '_' + str(depth)
                    band_headers.append(column_header)
                if agreement == False:
                    #print '{} disagrees at depth {}'.format(test_variable, depth)
                    diffs_exist = True
                    diffs_band = abs(all_test_data[cell][test_variable][depth] - all_base_data[cell][base_variable][depth])
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
                    print '    ' + variable + ': ' + str(agreement)
                else:
                    print '    ' + variable + ': False '
                    print '      Differences at depths ' + str(diffs_depths) + ' ',
                    print 'Number of different entries across all bands: {} '.format(num_diffs), 
                    print 'Maximum absolute difference across all bands: {}'.format(max_diff) 
    
    if csv_out == True:
        if csv_diffs_only == False:
            test_csv_filename = 'tabular_cell_{}_{}_test.csv'.format(cell, os.path.basename(testfile))
            np.savetxt(test_csv_filename, test_table, delimiter=',', fmt='%3.22f', header=",".join(general_headers), comments='')
            base_csv_filename = 'tabular_cell_{}_{}_base.csv'.format(cell, os.path.basename(basefile))
            np.savetxt(base_csv_filename, base_table, delimiter=',', fmt='%3.22f', header=",".join(general_headers), comments='')

            test_4D_csv_filename = 'tabular_cell_{}_{}_band_test.csv'.format(cell, os.path.basename(testfile)) 
            np.savetxt(test_4D_csv_filename, test_band_table, delimiter=',', fmt='%3.22f', header=",".join(band_headers), comments='')
            base_4D_csv_filename = 'tabular_cell_{}_{}_band_base.csv'.format(cell, os.path.basename(basefile))
            np.savetxt(base_4D_csv_filename, base_band_table, delimiter=',', fmt='%3.22f', header=",".join(band_headers), comments='')

        if diffs_exist:
            diffs_3D_csv_filename = 'tabular_cell_{}_{}_differences_tol={}.csv'.format(cell, os.path.basename(testfile), tolerance)
            np.savetxt(diffs_3D_csv_filename, diffs_table, delimiter=',', fmt='%3.22f', header=",".join(diffs_headers), comments='')     
            diffs_4D_csv_filename = 'tabular_cell_{}_{}_band_differences_tol={}.csv'.format(cell, os.path.basename(testfile), tolerance)
            np.savetxt(diffs_4D_csv_filename, diffs_band_table, delimiter=',', fmt='%3.22f', header=",".join(diffs_band_headers), comments='')

if diffs_exist:
    print '\nEquivalence test FAILED. Differences exist between testfile and basefile at the given tolerance of {}\n'.format(tolerance)
else:
    print '\nEquivalence test PASSED. The testfile and basefile are in agreement within the given tolerance of {}.\n'.format(tolerance)

print '\nvic_output_compare_netcdf_universal finished.'

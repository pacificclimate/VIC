#!/usr/bin/env python

""" This script checks any two VIC NetCDF output files against each other, even if one is
 time-major and the other is cell-major (just set the command line parameters accordingly).
  
 For usage information: ./vic_output_compare_netcdf_universal.py --help """
 
import sys
import argparse
import os.path
# import h5py
import netCDF4
import numpy as np
from collections import defaultdict

# need this in order to nest defaultdict objects beyond 2 levels
def tree(): return defaultdict(tree)

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
parser.add_argument('--debug', action="store_true", dest="debug", default=False, help = 'for debug output')
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
debug = options.debug
var_map_file = options.var_map_file

print('\n------- vic_output_compare_netcdf_universal -------')
# print interpretation of input
print('testfile: {} \nbasefile: {}\n'.format(testfile, basefile))

if test_is_time_major:
    print('testfile declared as time-major with depth dimension in position {}'.format(pos_depth_dim_test))
else:
    print('testfile declared as cell-major with depth dimension in position {}'.format(pos_depth_dim_test))

if base_is_time_major:
    print('basefile declared as time-major with depth dimension in position {}\n'.format(pos_depth_dim_base))
else:
    print('basefile declared as cell-major with depth dimension in position {}\n'.format(pos_depth_dim_base))

print('Absolute tolerance (-tolerance) for agreement between elements: {}\n'.format(tolerance))

if csv_out:
    print('CSV output selected (--csv).')
    if csv_diffs_only:
        print('Only differences between input files will be saved (-csv_diffs_only).')

# Open test output NetCDF file
testNC = netCDF4.Dataset(testfile, 'r').variables

# Open out_base NetCDF file
baseNC = netCDF4.Dataset(basefile, 'r').variables

# grab NaN fill value of non-initialized NetCDF records
for var in baseNC:
    try:
        fill_value_base = baseNC[var]._FillValue
    except:
        continue
    if fill_value_base is not None:
        break
print('basefile NetCDF _FillValue: {}'.format(fill_value_base))
for var in testNC:
    try:
        fill_value_test = testNC[var]._FillValue
    except:
        continue
    if fill_value_test is not None:
        break
print('testfile NetCDF _FillValue: {}'.format(fill_value_test))
if csv_out:
    print('Empty records in output CSV files will be assigned the basefile NetCDF _FillValue.')

# grab the number of time records
test_time_len = len(testNC['time'])
base_time_len = len(baseNC['time'])

if not test_end_rec:
    test_end_rec = test_time_len
    
if not base_end_rec:
    base_end_rec = base_time_len
    
num_test_recs = test_end_rec - test_start_rec
num_base_recs = base_end_rec - base_start_rec 

print('\ntestfile data time record range to be read: {}:{} ({} records total)'.format(test_start_rec, test_end_rec, num_test_recs))
print('basefile data time record range to be read: {}:{} ({} records total)'.format(base_start_rec, base_end_rec, num_base_recs))
if num_test_recs != num_base_recs:
    print('Number of time records selected for comparison between testfile and basefile must be equal.  Exiting.\n')
    sys.exit(0)

# grab lat and lon dimensions of grid cells
lats = testNC['lat'][:]
lons = testNC['lon'][:]

lat_to_idx = dict([(x[1],x[0]) for x in enumerate(lats.tolist())])
lon_to_idx = dict([(x[1],x[0]) for x in enumerate(lons.tolist())])

# initialize with keys from NetCDF file, and empty values
test_data_keys = dict.fromkeys(testNC.keys()) 
del test_data_keys['lat']
del test_data_keys['lon']
del test_data_keys['time']
del test_data_keys['depth']
if 'bnds' in test_data_keys.keys():
    del test_data_keys['bnds']
if debug:        
    print('test_data_keys loaded: {}'.format(test_data_keys))

# initialize with keys from NetCDF file, and empty values
base_data_keys = dict.fromkeys(baseNC.keys()) 
del base_data_keys['lat']
del base_data_keys['lon']
del base_data_keys['time']
del base_data_keys['depth']
if 'bnds' in base_data_keys.keys():
    del base_data_keys['bnds']
if debug:        
    print('base_data_keys loaded: {}'.format(base_data_keys))

# define basefile to testfile variable name mapping (e.g. precipitation might be called
# OUT_PREC, PREC, or pr in different files) and build up a vars_to_test nested dict like:
# {OUT_PREC:{'test_var': 'PREC', 'base_var': 'pr'}, OUT_WIND: ...} using the "OUT_" form as the key
var_map = {}
vars_to_test = {}
if var_map_file:
    print('\nOutput variable name mapping options given in file {} will be used.'.format(var_map_file))
    with open(var_map_file, 'r') as f:
        for line in f:
            if (len(line) > 1):
                split_line = line.split()
                var = split_line[0]
                var_map[var] = []
                for idx, mapped_var_name in enumerate(split_line):
                    if idx > 0:
                        var_map[var].append(mapped_var_name)
    for variable in var_map:
        test_var_found = False
        base_var_found = False
        if variable in test_data_keys:
            vars_to_test[variable] = {}
            vars_to_test[variable]['test_var'] = variable
            test_var_found = True
        else: # look for other matches in var_map loaded from file
            for var_name in var_map[variable]:
                if var_name in test_data_keys:
                    vars_to_test[variable] = {}
                    vars_to_test[variable]['test_var'] = var_name
                    test_var_found = True
        if not test_var_found:
            continue
        if variable in base_data_keys:
            vars_to_test[variable]['base_var'] = variable
            base_var_found = True
        else: # look for other matches in var_map loaded from file
            for var_name in var_map[variable]:
                if var_name in base_data_keys:
                    vars_to_test[variable]['base_var'] = var_name
                    base_var_found = True
        # TODO: add an option to skip comparison of vars that do not appear in the basefile
        if not base_var_found:
            print('No entry for {} output variable found in basefile.  Exiting.\n'.format(variable))
            sys.exit(0)
else: # use the variable names provided in the testfile, which assumes the basefile variable naming agrees
    for variable in test_data_keys:
        vars_to_test[variable] = {}
        vars_to_test[variable]['test_var'] = variable
        vars_to_test[variable]['base_var'] = variable
if debug:        
    print('vars_to_test: {}'.format(vars_to_test))

# create a big nested dictionary with all data for all cells, one for test data and one for out_base data
all_test_data = tree()
all_base_data = tree()

# load up all_test_data
for lat in lats:
    for lon in lons:
        cell_label = '{}_{}'.format(lat, lon)
        print('\nLoading cell {} data...'.format(cell_label))
        for variable in test_data_keys:
            if debug:
                print('loading variable from test file: {}'.format(variable))
            if test_is_time_major == True: # test file is time-major format
                if len(testNC[variable].shape) == 4: # 4D variable
                    if pos_depth_dim_test == 1: # <time, depth, lat, lon>
                        for depth in range(0, testNC[variable].shape[pos_depth_dim_test]):
                            all_test_data[cell_label][variable][depth] = testNC[variable][test_start_rec:test_end_rec,depth,lat_to_idx[lat],lon_to_idx[lon]]
                    elif pos_depth_dim_test == 3: # <time, lat, lon, depth>
                        for depth in range(0, testNC[variable].shape[pos_depth_dim_test]):
                        #print 'cell_label: {} variable: {} depth: {}'.format(cell_label, variable, depth)
                            all_test_data[cell_label][variable][depth] = testNC[variable][test_start_rec:test_end_rec,lat_to_idx[lat],lon_to_idx[lon],depth]
                    else:
                        print('The declared depth dimension position {} for a time-major NetCDF test file is not supported.'.format(pos_depth_dim_test))
                        sys.exit(0)
                elif len(testNC[variable].shape) == 3: # 3D variable
                    all_test_data[cell_label][variable] = testNC[variable][test_start_rec:test_end_rec,lat_to_idx[lat],lon_to_idx[lon]]

            elif test_is_time_major == False: # test file is cell-major format
                if len(testNC[variable].shape) == 4: # 4D variable
                    if pos_depth_dim_test == 0: # <depth, lat, lon, time>
                        for depth in range(0, testNC[variable].shape[pos_depth_dim_test]):
                            all_test_data[cell_label][variable][depth] = testNC[variable][depth,lat_to_idx[lat],lon_to_idx[lon],test_start_rec:test_end_rec]
                    else:
                        print('The declared depth dimension position {} for a cell-major NetCDF test file is not supported.'.format(pos_depth_dim_test))
                        sys.exit(0)
                elif len(testNC[variable].shape) == 3: # 3D variable
                    all_test_data[cell_label][variable] = testNC[variable][lat_to_idx[lat],lon_to_idx[lon],test_start_rec:test_end_rec]
            if fill_value_test != fill_value_base:
                # overwrite all testfile NetCDF _FillValue (NaNs) entries read in to all_test_data with the basefile _FillValue          
                if len(testNC[variable].shape) == 3: # 3D variable
                    all_test_data[cell_label][variable][all_test_data[cell_label][variable] == fill_value_test] = fill_value_base
                elif len(testNC[variable].shape) == 4: # 4D variable
                    for depth in range(0, testNC[variable].shape[pos_depth_dim_test]):
                        all_test_data[cell_label][variable][depth][all_test_data[cell_label][variable][depth] == fill_value_test] = fill_value_base
        for variable in base_data_keys:
            if debug:
                print('loading variable from base file: {}'.format(variable))
            if base_is_time_major == True: # base file is time-major format
                if len(baseNC[variable].shape) == 4: # 4D variable
                    if pos_depth_dim_base == 1: # <time, depth, lat, lon>
                        for depth in range(0, baseNC[variable].shape[pos_depth_dim_base]):
                            all_base_data[cell_label][variable][depth] = baseNC[variable][base_start_rec:base_end_rec,depth,lat_to_idx[lat],lon_to_idx[lon]]
                    elif pos_depth_dim_base == 3: # <time, lat, lon, depth>
                        for depth in range(0, baseNC[variable].shape[pos_depth_dim_base]):
                            all_base_data[cell_label][variable][depth] = baseNC[variable][base_start_rec:base_end_rec,lat_to_idx[lat],lon_to_idx[lon],depth]
                    else:
                        print('The declared depth dimension position {} for a time-major NetCDF base file is not supported.'.format(pos_depth_dim_base))
                        sys.exit(0)
                elif len(baseNC[variable].shape) == 3: # 3D variable
                    all_base_data[cell_label][variable] = baseNC[variable][base_start_rec:base_end_rec,lat_to_idx[lat],lon_to_idx[lon]]

            elif base_is_time_major == False: # base file is cell-major format
                if len(baseNC[variable].shape) == 4: # 4D variable
                    if pos_depth_dim_base == 0: # <depth, lat, lon, time>
                        for depth in range(0, baseNC[variable].shape[pos_depth_dim_base]):
                            all_base_data[cell_label][variable][depth] = baseNC[variable][depth,lat_to_idx[lat],lon_to_idx[lon],base_start_rec:base_end_rec]
                    else:
                        print('The declared depth dimension position {} for a cell-major NetCDF base file is not supported.'.format(pos_depth_dim_base))
                        sys.exit(0)
                elif len(baseNC[variable].shape) == 3: # 3D variable
                    all_base_data[cell_label][variable] = baseNC[variable][lat_to_idx[lat],lon_to_idx[lon],base_start_rec:base_end_rec]
                              

## compare test data against base data
print('\nChecking agreement between testfile variables data and those in the basefile...\n (Note: this will only check the set \
of variables existing in your testfile against those same ones in the basefile, which may have a larger set of variables than testfile)')

# create labels for each cell to use to index into defaultdicts we created
cell_labels = []
for lat_label in lats:
    for lon_label in lons:
        cell_labels.append(repr(lat_label) + '_' + repr(lon_label))

# global flag indicating whether there are *any* differences between the files
diffs_exist = False

for cell in cell_labels:
    # if we want to output the test dataset to CSV format for inspection
    if csv_out == True: 
        # Create tables that will be written into CSV files
        # time record index read from file will form the first column of the output tables,
        # relative time record index (0 at test/base_start_rec position) in second column
        general_headers = ['time_rec','relative_time_rec']
        test_table = testNC['time'][test_start_rec:test_end_rec]
        test_table = np.column_stack([test_table, range(0,num_test_recs)])
        base_table = baseNC['time'][base_start_rec:base_end_rec]
        base_table = np.column_stack([base_table, range(0,num_base_recs)])
        band_headers = ['time_rec', 'relative_time_rec'] 
        test_band_table = testNC['time'][test_start_rec:test_end_rec]
        test_band_table = np.column_stack([test_band_table, range(0,num_test_recs)])
        base_band_table = baseNC['time'][base_start_rec:base_end_rec]
        base_band_table = np.column_stack([base_band_table, range(0,num_base_recs)])
        diffs_headers = ['time', 'relative_time_rec']
        diffs_table = baseNC['time'][base_start_rec:base_end_rec]
        diffs_table = np.column_stack([diffs_table, range(0,num_base_recs)])
        diffs_band_headers = ['time', 'relative_time_rec']
        diffs_band_table = baseNC['time'][base_start_rec:base_end_rec]
        diffs_band_table = np.column_stack([diffs_band_table, range(0,num_base_recs)])

    if verbose:    
        print('\n')
        print('Cell {}\n'.format(str(cell)))
       # print 'var_map: {}'.format(var_map)
        print('\tVariable\tAgrees\tNum diffs\t\tMax abs diff\t\tSum of diffs')
        print('\t'+'-'*145)

    # Check agreement between all variables present in the testfile and their equivalent in the basefile
    for variable in vars_to_test:
        agreement = False
        diffs = []
        diffs_depths = []

        if len(testNC[vars_to_test[variable]['test_var']].shape) == 3: # 3D variable
            if tolerance > 0:
                agreement = np.allclose(all_test_data[cell][vars_to_test[variable]['test_var']], all_base_data[cell][vars_to_test[variable]['base_var']], 0, tolerance)
            else:
                agreement = np.array_equal(all_test_data[cell][vars_to_test[variable]['test_var']], all_base_data[cell][vars_to_test[variable]['base_var']])
            if csv_out == True:
                column_header = variable
                general_headers.append(column_header)
                test_table = np.column_stack([test_table, all_test_data[cell][vars_to_test[variable]['test_var']]])
                base_table = np.column_stack([base_table, all_base_data[cell][vars_to_test[variable]['base_var']]])
            if agreement == False:
                diffs_exist = True
                diffs = abs(all_test_data[cell][vars_to_test[variable]['test_var']] - all_base_data[cell][vars_to_test[variable]['base_var']])
                if csv_out == True:
                    diffs_table = np.column_stack([diffs_table, diffs])        
                    diffs_headers.append(column_header) 
                num_diffs = len(diffs[diffs > tolerance])
                max_diff = np.max(diffs)
                sum_diffs = np.sum(diffs)
                if verbose:
                    print('\t{}\t{}\t{}\t\t{}\t\t{}'.format(variable, str(agreement), num_diffs, max_diff, sum_diffs))
            else:
                if verbose:
                    print('\t{}\t{}'.format(variable, str(agreement)))
        elif len(testNC[vars_to_test[variable]['test_var']].shape) == 4: # 4D variable
            max_diff = 0
            num_diffs = 0
            for depth in range(0, testNC[vars_to_test[variable]['test_var']].shape[pos_depth_dim_test]):
                if tolerance > 0:
                    agreement = np.allclose(all_test_data[cell][vars_to_test[variable]['test_var']][depth], all_base_data[cell][vars_to_test[variable]['base_var']][depth], 0, tolerance)
                else:
                    agreement = np.array_equal(all_test_data[cell][vars_to_test[variable]['test_var']][depth], all_base_data[cell][vars_to_test[variable]['base_var']][depth])
                if csv_out == True:
                    test_band_table = np.column_stack([test_band_table, all_test_data[cell][vars_to_test[variable]['test_var']][depth]])
                    base_band_table = np.column_stack([base_band_table, all_base_data[cell][vars_to_test[variable]['base_var']][depth]])
                    column_header = variable + '_' + str(depth)
                    band_headers.append(column_header)
                if agreement == False:
                    diffs_exist = True
                    diffs_band = abs(all_test_data[cell][vars_to_test[variable]['test_var']][depth] - all_base_data[cell][vars_to_test[variable]['base_var']][depth])
                    if csv_out == True:
                        diffs_band_table = np.column_stack([diffs_band_table, diffs_band])        
                        diffs_band_headers.append(column_header)
                    diffs_depths.append(depth)
                    max_temp = np.max(diffs_band) 
                    max_diff = max_temp if max_temp > max_diff else max_diff # running max of differences across all bands
                    num_diffs += len(diffs_band[diffs_band > tolerance]) # running number of differences across all bands
            if verbose:    
                if not diffs_depths: # no differences were found at any depth
                    print('\t{}\t{}'.format(variable, str(agreement)))
                else:
                    print('\t{}\t{}\t [4D] Diffs at bands {} Num diffs across bands: {} Max abs diff across bands: {}'.format(variable, str(agreement), str(diffs_depths), num_diffs, max_diff))
                    
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
    print('\nEquivalence test FAILED. Differences exist between testfile and basefile at the given tolerance of {}\n'.format(tolerance))
else:
    print('\nEquivalence test PASSED. The testfile and basefile are in agreement within the given tolerance of {}.\n'.format(tolerance))

print('\nvic_output_compare_netcdf_universal finished.')

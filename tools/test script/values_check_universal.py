#!/usr/bin/env python

""" This script checks any two VIC NetCDF output files against each other, even if one is
 time-major and the other is cell-major (just set the input parameters accordingly).
  
 Usage: ./values_check_universal.py <path to test file> <test time-major?> <path to base file> <base time-major?>
 where time-majorness is indicated by a 't' after each input file"""
 
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
#parser = argparse.ArgumentParser()
parser.add_argument('--check', action="store", dest="testfile", type=str, help = 'file name and path of the NetCDF file to be checked')
parser.add_argument('-ctm', action="store_true", dest="test_is_time_major", default=True, help = 'file to be checked is time-major (default=True)')
parser.add_argument('-ccm', action="store_false", dest="test_is_time_major", default=True, help = 'file to be checked is cell-major (default=False)')
parser.add_argument('-cdp', action="store", dest="pos_depth_dim_test", type=int, default=1, help = 'set the 0-based position of the depth dimension of 4D variables in the NetCDF file to be checked (default=1, assuming time-major test file option [-ctm] is set, i.e. <time, depth, lat, lon>)')
parser.add_argument('--base', action="store", dest="basefile", type=str, help = 'file name and path of the base NetCDF file to check against')
parser.add_argument('-btm', action="store_true", dest="base_is_time_major", default=False, help = 'base file is time-major (default=False)')
parser.add_argument('-bcm', action="store_false", dest="base_is_time_major", default=False, help = 'base file is cell-major (default=True)')
parser.add_argument('-bdp', action="store", dest="pos_depth_dim_base", type=int, default=0, help = 'set the 0-based position of the depth dimension of 4D variables in the base NetCDF file (default=0, assuming cell-major base file option [-bcm] is set, i.e. <depth, lat, lon, time>)')
parser.add_argument('--csv', action="store_true", dest="csv_out", default=False, help = 'set if you want output of data read from NetCDF to CSV')
parser.add_argument('-csv_diffs_only', action="store_true", dest="csv_diffs_only", default=False, help = 'set if you want to output only the differences to CSV and not the original data')
parser.add_argument('--v', action="store_true", dest="verbose", default=False, help = 'for verbose output')
parser.add_argument('-start', action="store", dest="start_range", type=int, default=0, help = 'when in verbose mode, use this to specify the start index of the time window you want to view data from')
parser.add_argument('-end', action="store", dest="end_range", type=int, default=5, help = 'when in verbose mode, use this to specify the end index of the time window you want to view data from')
parser.add_argument('-tol', action="store", dest="tolerance", type=float, default=0, help = 'set the absolute tolerance that all values must be within to be considered in agreement')
parser.add_argument('-depth_check', action="store", dest="depth_check", type=int, default=0, help = 'set to choose which depth of 4-dimensional variables to be shown when in verbose mode (default=0)')

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
verbose = options.verbose
start_range = options.start_range
end_range = options.end_range
tolerance = options.tolerance
depth_check = options.depth_check


### Get data from test output
#h5 = h5py.File('/home/mfischer/vic_dev/out/testing/results.nc', 'r')
inputH5 = h5py.File(testfile, 'r')

# print interpretation of input
print 'Checking test file {} against base file {}'.format(testfile, basefile)

if test_is_time_major:
  print 'Test file declared as time-major'
else:
  print 'Test file declared as cell-major'

# out_base file
#out_base_file = h5py.File('/home/mfischer/vic_dev/out_base/place/netcdf/false/results.nc')
baseH5 = h5py.File(basefile, 'r')

if base_is_time_major:
  print 'Base file declared as time-major'
else:
  print 'Base file declared as cell-major'

print 'Verbose output (--v): {}'.format(verbose)
if verbose == True:
  print 'Range of indices for verbose variable inspection (-start_range, -end_range): {}:{}'.format(start_range, end_range)
  print '4-dimensional data will be shown at this range for depth (-depth_check): {}'.format(depth_check)

print 'Absolute tolerance (-tol) for agreement between elements: {}'.format(tolerance)

if csv_out == True:
  print 'CSV output selected (--csv).'
  if csv_diffs_only == True:
    print 'Only differences between input files will be saved (-csv_diffs_only).'

# grab the number of time records
time_len = len(inputH5['time'])

# grab the depth of 4-dimension variables
#depth = len(inputH5['depth'])

# grab NaN fill value of non-initialized NetCDF records from one standard attribute
fill_value = baseH5['pr'].attrs['_FillValue']

# grab lat and lon dimensions of grid cells
lats = inputH5['lat'][:]
lons = inputH5['lon'][:]

lat_to_idx = dict([(x[1],x[0]) for x in enumerate(lats.tolist())])
lon_to_idx = dict([(x[1],x[0]) for x in enumerate(lons.tolist())])

# initialize with keys from NetCDF file, and empty values
cell_data_keys = dict.fromkeys(inputH5.keys()) 
del cell_data_keys['lat']
del cell_data_keys['lon']
del cell_data_keys['time']
del cell_data_keys['bnds']
del cell_data_keys['depth']

# need this in order to nest defaultdict objects beyond 2 levels
def tree(): return defaultdict(tree)

# create one big nested dictionary with all data for all cells
all_test_data = tree()
all_base_data = tree()

# load up all_test_data
for lat in lats:
    print 'depth_check: {}'.format(depth_check)
    for lon in lons:
        cell_label = '{}_{}'.format(lat, lon)
        if verbose:
             print 'Loading cell {} data...'.format(cell_label)
        for variable in cell_data_keys:
            if test_is_time_major == True: # test file is time-major format
                if len(inputH5[variable].shape) == 4: # 4D variable
		    if pos_depth_dim_test == 1: # <time, depth, lat, lon>
			for depth in range(0, inputH5[variable].shape[pos_depth_dim_test]):
                            all_test_data[cell_label][variable][depth] = inputH5[variable][:,depth,lat_to_idx[lat],lon_to_idx[lon]]
		    elif pos_depth_dim_test == 3: # <time, lat, lon, depth>
			for depth in range(0, inputH5[variable].shape[pos_depth_dim_test]):
			    #print 'cell_label: {} variable: {} depth: {}'.format(cell_label, variable, depth)
                            all_test_data[cell_label][variable][depth] = inputH5[variable][:,lat_to_idx[lat],lon_to_idx[lon],depth]
		    else:
			print 'The declared depth dimension position {} for a time-major NetCDF test file is not supported.'.format(pos_depth_dim_test)
			sys.exit(0)
                elif len(inputH5[variable].shape) == 3: # 3D variable
                    all_test_data[cell_label][variable] = inputH5[variable][:,lat_to_idx[lat],lon_to_idx[lon]]

            elif test_is_time_major == False: # test file is cell-major format
                if len(inputH5[variable].shape) == 4: # 4D variable
		    if pos_depth_dim_test == 0: # <depth, lat, lon, time>
			for depth in range(0, inputH5[variable].shape[pos_depth_dim_test]):
	                    all_test_data[cell_label][variable][depth] = inputH5[variable][depth,lat_to_idx[lat],lon_to_idx[lon],:]
		    else:
			print 'The declared depth dimension position {} for a cell-major NetCDF test file is not supported.'.format(pos_depth_dim_test)
			sys.exit(0)
                elif len(inputH5[variable].shape) == 3: # 3D variable
                    all_test_data[cell_label][variable] = inputH5[variable][lat_to_idx[lat],lon_to_idx[lon],:]
            if verbose:
		if len(inputH5[variable].shape) == 4: # 4D variable
                    print 'TEST cell: {} variable: {} data: {}'.format(cell_label, variable, all_test_data[cell_label][variable][depth_check][start_range:end_range])
		    pass
                elif len(inputH5[variable].shape) == 3: # 3D variable
		    print 'TEST cell: {} variable: {} data: {}'.format(cell_label, variable, all_test_data[cell_label][variable][start_range:end_range])
	    # overwrite all NaNs (fill values) with 0	    	
            if len(inputH5[variable].shape) == 3: # 3D variable
                all_test_data[cell_label][variable][all_test_data[cell_label][variable] == fill_value] = 0
            elif len(inputH5[variable].shape) == 4: # 4D variable
		for depth in range(0, inputH5[variable].shape[pos_depth_dim_test]):
                    all_test_data[cell_label][variable][depth][all_test_data[cell_label][variable][depth] == fill_value] = 0

            if base_is_time_major == True: # base file is time-major format
                if len(baseH5[variable].shape) == 4: # 4D variable
		    if pos_depth_dim_base == 1: # <time, depth, lat, lon>
			for depth in range(0, baseH5[variable].shape[pos_depth_dim_base]):
                            all_base_data[cell_label][variable][depth] = baseH5[variable][:,depth,lat_to_idx[lat],lon_to_idx[lon]]
		    elif pos_depth_dim_base == 3: # <time, lat, lon, depth>
			for depth in range(0, baseH5[variable].shape[pos_depth_dim_base]):
                            all_base_data[cell_label][variable][depth] = baseH5[variable][:,lat_to_idx[lat],lon_to_idx[lon],depth]
		    else:
			print 'The declared depth dimension position {} for a time-major NetCDF base file is not supported.'.format(pos_depth_dim_base)
			sys.exit(0)
                elif len(baseH5[variable].shape) == 3: # 3D variable
                    all_base_data[cell_label][variable] = baseH5[variable][:,lat_to_idx[lat],lon_to_idx[lon]]

            elif base_is_time_major == False: # base file is cell-major format
                if len(baseH5[variable].shape) == 4: # 4D variable
		    if pos_depth_dim_base == 0: # <depth, lat, lon, time>
			for depth in range(0, baseH5[variable].shape[pos_depth_dim_base]):
       	                #depth = baseH5[variable].shape[pos_depth_dim_base] # gets actual depth of this variable in the NetCDF
                            all_base_data[cell_label][variable][depth] = baseH5[variable][depth,lat_to_idx[lat],lon_to_idx[lon],:]
		    else:
			print 'The declared depth dimension position {} for a cell-major NetCDF base file is not supported.'.format(pos_depth_dim_base)
			sys.exit(0)
                elif len(baseH5[variable].shape) == 3: # 3D variable
                    all_base_data[cell_label][variable] = baseH5[variable][lat_to_idx[lat],lon_to_idx[lon],:]
            if verbose:
		if len(baseH5[variable].shape) == 4:
                    print 'BASE cell: {} variable: {} data: {}'.format(cell_label, variable, all_base_data[cell_label][variable][depth_check][start_range:end_range])
		elif len(baseH5[variable].shape) == 3:
		    print 'BASE cell: {} variable: {} data: {}'.format(cell_label, variable, all_base_data[cell_label][variable][start_range:end_range])
#                    raw_input("Press enter to continue.")
#                    continue
                  
	    # overwrite all NaNs (fill values) with 0	    	
            if len(baseH5[variable].shape) == 3: # 3D variable
                all_base_data[cell_label][variable][all_base_data[cell_label][variable] == fill_value] = 0
            elif len(baseH5[variable].shape) == 4: # 4D variable
		for depth in range(0, baseH5[variable].shape[pos_depth_dim_base]):
                    all_base_data[cell_label][variable][depth][all_base_data[cell_label][variable][depth] == fill_value] = 0



# create labels for each cell to use to index into defaultdicts we created
cell_labels = []
for lat_label in lats:
    for lon_label in lons:
        cell_labels.append(repr(lat_label) + '_' + repr(lon_label))
    
# sanity check
#var = 'pr'
#for cell_label in cell_labels:
#    print 'sanity check: variable: {} test values: {} base values: {}'.format(var, all_test_data[cell_label][var][0:5], all_base_data[cell_label][var][0:5])

## compare test data against base data
print'Checking agreement between test output data with base values:'
for cell in cell_labels:
    diffs_exist = False
    # if we want to output the test dataset to CSV format for inspection
    if csv_out == True:
        # time will form the first column of the output tables
        general_headers = ['time']
        test_table = inputH5['time'] 
        base_table = baseH5['time']
	band_headers = ['time'] 
        test_band_table = inputH5['time'] 
        base_band_table = baseH5['time'] 
        diffs_table = baseH5['time']
        diffs_headers = ['time']
        diffs_band_table = baseH5['time']
        diffs_band_headers = ['time']
        
    print'\n'
    print'Cell ' + str(cell) + ': \n'
    print ''
    for variable in cell_data_keys:
#	print 'checking agreement on cell: {} variable: {}'.format(cell, variable)
	if len(inputH5[variable].shape) == 3: # 3D variable
            if tolerance > 0:
                agreement = np.allclose(all_test_data[cell][variable], all_base_data[cell][variable], 0, tolerance)
	    else:
                agreement = np.array_equal(all_test_data[cell][variable], all_base_data[cell][variable])
            if agreement == False:
		diffs_exist = True
                diffs = abs(all_test_data[cell][variable] - all_base_data[cell][variable])
		num_diffs = len(diffs)
		max_diff = np.max(diffs)
                print '    ' + str(inputH5[variable].attrs['internal_vic_name']) + ': ' + str(agreement),
	        print '(Number of different entries: {} '.format(num_diffs), 
	        print 'Maximum absolute difference: {})'.format(max_diff) 
	    else:
                print '    ' + str(inputH5[variable].attrs['internal_vic_name']) + ': ' + str(agreement)

	elif len(inputH5[variable].shape) == 4: # 4D variable
	    diffs_depths = []
 	    max_diff = 0
            num_diffs = 0
	    for depth in range(0, inputH5[variable].shape[pos_depth_dim_test]):
                if tolerance > 0:
                    agreement = np.allclose(all_test_data[cell][variable][depth], all_base_data[cell][variable][depth], 0, tolerance)
	        else:
                    agreement = np.array_equal(all_test_data[cell][variable][depth], all_base_data[cell][variable][depth], 0, tolerance)
                if agreement == False:
	            diffs_exist = True
	            diffs_depths.append(depth)
                    diffs_band = abs(all_test_data[cell][variable][depth] - all_base_data[cell][variable][depth])
		    diff = np.max(diffs_band) 
                    max_diff = diff if diff > max_diff else max_diff # running max of differences across all bands
		    num_diffs += len(diffs_band) # running number of differences across all bands
	    if not diffs_depths: # no differences were found at any depth
                print '    ' + str(inputH5[variable].attrs['internal_vic_name']) + ': ' + str(agreement)
	    else:
                print '    ' + str(inputH5[variable].attrs['internal_vic_name']) + ': False ' + '(Differences at depths ' + str(diffs_depths) + ' ',
	        print 'Number of different entries across all bands: {} '.format(num_diffs), 
	        print 'Maximum absolute difference across all bands: {})'.format(max_diff) 

        if csv_out == True:
            if len(baseH5[variable].shape) == 4: # write out all layers to another CSV file
		# Note: this assumes that 4D variables in the test and base files have the same depths
		for depth in range(0, inputH5[variable].shape[pos_depth_dim_test]):
               	    test_band_table = np.column_stack([test_band_table, all_test_data[cell][variable][depth]])
	            base_band_table = np.column_stack([base_band_table, all_base_data[cell][variable][depth]])
		    column_header = str(inputH5[variable].attrs['internal_vic_name']) + '_' + str(depth)
		    band_headers.append(column_header)
		    if agreement == False:	
        	        diffs_band_table = np.column_stack([diffs_band_table, diffs_band[depth]])        
	                diffs_band_headers.append(column_header)
            elif len(baseH5[variable].shape) == 3:
		column_header = str(inputH5[variable].attrs['internal_vic_name'])
                general_headers.append(column_header)
                test_table = np.column_stack([test_table, all_test_data[cell][variable]])
                base_table = np.column_stack([base_table, all_base_data[cell][variable]])
                if agreement == False:
                    diffs_table = np.column_stack([diffs_table, diffs])        
                    diffs_headers.append(column_header) 
                
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

        if diffs_exist == True:
            diffs_3D_csv_filename = 'tabular_cell_{}_{}_differences_tol={}.csv'.format(cell, os.path.basename(testfile), tolerance)
            np.savetxt(diffs_3D_csv_filename, diffs_table, delimiter=',', fmt='%3.22f', header=",".join(diffs_headers))	    
            diffs_4D_csv_filename = 'tabular_cell_{}_{}_band_differences_tol={}.csv'.format(cell, os.path.basename(testfile), tolerance)
            np.savetxt(diffs_4D_csv_filename, diffs_band_table, delimiter=',', fmt='%3.22f', header=",".join(diffs_band_headers))

print 'Finished.'

#!/usr/bin/env python
from datetime import datetime
from datetime import timedelta
import struct
import numpy as np
from numpy import dtype
import netCDF4
from netCDF4 import stringtochar
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs='*')
    parser.add_argument("-o", "--output", help="output file name", required=True)
    args = parser.parse_args()
    converter = Merger()
    converter.convert(args)

class Merger() :
    def convert(self, args):
        self.files = args.files
        self.output = args.output
        self.create_output()

    def create_output(self) :
        self.create_base()
        self.output_file()

    def create_base(self) :
        #
        # create base ouput file from first input file
        #
        infile = self.files[0]
        self.nsweep = len(self.files)
        self.nc = netCDF4.Dataset(self.output, "w", format="NETCDF4")
        # set global attribute
        self.write_global_attributes()
        self.write_global_variables()

    def write_global_attributes(self):
        nc = self.nc
        src = netCDF4.Dataset(self.files[0], "r", format="NETCDF4")
        for name in src.ncattrs():
            nc.setncattr(name, src.getncattr(name))
        src.close()

    def write_global_variables(self):
        nc = self.nc
        nc.createDimension("string_length", None)
        nc.createDimension("frequency", 1)
        global_variable_name = [ "volume_number", " time_reference", "latitude", "longitude", "altitude", "time_reference", "time_coverage_start", "time_coverage_end", "range" ]
        nc = self.nc
        src = netCDF4.Dataset(self.files[0], "r", format="NETCDF4")
        self.nrange= len(src.dimensions["range"])
        nc.createDimension("range", self.nrange )
        for name, variable in src.variables.items():
            if name in global_variable_name :
                x = nc.createVariable(name, variable.datatype, variable.dimensions)
                nc.variables[name][:] = src.variables[name][:]
                nc[name].setncatts(src[name].__dict__)
        for name, variable in src.variables.items():
            if name in ["frequency"]:
                x = nc.createVariable(name, variable.datatype, variable.dimensions)
                nc.variables[name][:] = src.variables[name][:]
                nc[name].setncatts(src[name].__dict__)
        src.close()
        # set end volume scan
        src = netCDF4.Dataset(self.files[-1], "r", format="NETCDF4")
        nc["time_coverage_end"][:] = src["time_coverage_end"][:]
        src.close()
        #

    def output_file(self):
        nc = self.nc
        nsweep = self.nsweep
        nc.createDimension("sweep", nsweep)
        ntime = 0
        for n in range(nsweep):
            src = netCDF4.Dataset(self.files[n], "r", format="NETCDF4")
            ntime = ntime +  len(src.dimensions["time"])
            if self.nrange != len(src.dimensions["range"])  :
                    print('Error: range is not same in files', file=sys.stderr)
                    sys.exit(1)
            src.close()
        nc.createDimension("time", ntime)
        ns = 0
        ne = 0
        time_variable_name = ["azimuth", "elevation", "time" ]
        sweep_variable_name = ["fixed_angle", "sweep_mode", "sweep_number","sweep_start_ray_index", "sweep_end_ray_index"  ]
        varnames = sweep_variable_name + time_variable_name
        for n in range(nsweep) :
            src = netCDF4.Dataset(self.files[n], "r", format="NETCDF4")
            # copy attribute and create variables
            if n == 0 :
                for name, variable in src.variables.items():
                    if name in varnames :
                        x = nc.createVariable(name, variable.datatype, variable.dimensions)
                        nc[name].setncatts(src[name].__dict__)
                    if variable.dimensions == ('time', 'range') :
                        x = nc.createVariable(name, variable.datatype, variable.dimensions)
                        nc[name].setncatts(src[name].__dict__)
            ne +=  len(src.dimensions["time"])
            for name, variable in src.variables.items():
                if name in time_variable_name :
                    nc.variables[name][ns:ne] = src.variables[name][:]
                if name in sweep_variable_name :
                    if name == "sweep_number":
                        nc.variables[name][n] = n
                    elif name == "sweep_start_ray_index" :
                        nc.variables[name][n] = ns
                    elif name == "sweep_end_ray_index" :
                        nc.variables[name][n] = ne - 1
                    else:
                        nc.variables[name][n] = src.variables[name][:]
                if variable.dimensions == ('time', 'range') :
                    nc.variables[name][ns:ne, : ] = src.variables[name][:]
            ns += len(src.dimensions["time"])
            #print(ns)
            src.close()
        nc.close()

if __name__ == "__main__":
    main()


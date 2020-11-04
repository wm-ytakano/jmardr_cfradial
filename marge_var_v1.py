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
        self.append_file()

    def create_base(self) :
        #
        # create base ouput file from first input file
        #
        infile = self.files[0]
        self.nfile = len(self.files)
        self.nc = netCDF4.Dataset(self.output, "w", format="NETCDF4" )
        # set global attribute
        self.copy_file()


    def copy_file(self):
        nc = self.nc
        src = netCDF4.Dataset(self.files[0], "r", format="NETCDF4")
        # copy global attribute
        for name in src.ncattrs():
            nc.setncattr(name, src.getncattr(name))
        # copy dimension
        for name, dimension in src.dimensions.items():
            nc.createDimension( name, (len(dimension) if not dimension.isunlimited() else None))
        # copy variable 
        for name, variable in src.variables.items():
            x = nc.createVariable(name, variable.datatype, variable.dimensions, zlib=True)
            nc[name].setncatts(src[name].__dict__)
            nc.variables[name][:] = src.variables[name][:]
        for name in src.ncattrs():
            nc.setncattr(name, src.getncattr(name))
        src.close()
        self.nrange = len(nc.dimensions["range"])

    def append_file(self):
        nc = self.nc
        nvar = self.nfile
        for n in range(1,nvar) :
            src = netCDF4.Dataset(self.files[n], "r", format="NETCDF4")
            # copy attribute and create variables
            for name, variable in src.variables.items():
                # error check
                nrange = len(src.dimensions["range"])
                if nrange != self.nrange :
                    print("range is not same in files", file=sys.stderr)
                    sys.exit(1)
                if variable.dimensions == ('time', 'range') :
                    x = nc.createVariable(name, variable.datatype, variable.dimensions, zlib=True)
                    nc[name].setncatts(src[name].__dict__)
                    nc.variables[name][:] = src.variables[name][:]
            src.close()
        nc.close()

if __name__ == "__main__":
    main()


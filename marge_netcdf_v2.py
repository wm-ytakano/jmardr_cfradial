#!/usr/bin/env python3
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
    parser.add_argument("-w", "--work", help="path to working directory") 
    parser.add_argument("-o", "--output", help="output file name", required=True)    
    args = parser.parse_args()
    converter = Merger()
    converter.convert(args)
    
class Merger() :
    def convert(self, args):
        self.files = args.files
        if args.work:
            self.workdir = args.work
        else:
            self.workdir = "."
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
        # change version
        nc.setncattr("version", "2.0")
        
    def write_global_variables(self):
        nc = self.nc        
        nc.createDimension("sweep", self.nsweep)
        nc.createDimension("frequency", 1)        
        nc.createDimension("string_length", None)
        global_variable_name = [ "volume_number",  " time_reference", "latitude", "longitude", "altitude", "time_coverage_start", "time_coverage_end"]
        nc = self.nc        
        src = netCDF4.Dataset(self.files[0], "r", format="NETCDF4")
        for name, variable in src.variables.items():
            if name in global_variable_name :
                x = nc.createVariable(name, variable.datatype, variable.dimensions)
                nc.variables[name][:] = src.variables[name][:]          
                nc[name].setncatts(src[name].__dict__)
        # treate  radar_parameters
        gr = nc.createGroup("radar_parameters")
        for name, variable in src.variables.items():
            if name in ["frequency"]:
                x = gr.createVariable(name, variable.datatype, variable.dimensions)
                gr.variables[name][:] = src.variables[name][:]          
                gr[name].setncatts(src[name].__dict__)
        # set end volume scan
        src = netCDF4.Dataset(self.files[-1], "r", format="NETCDF4")        
        nc["time_coverage_end"][:] = src["time_coverage_end"][:]
        src.close()
        # 
        x = nc.createVariable("sweep_group_name", "str", ("sweep"))
        x.long_name = "group_name_for_sweep"
        for n in range(self.nsweep) :
            x[n] = "sweep_"+str(n).zfill(4) 
    def output_file(self) :
        nc = self.nc
        dimension_name = [ "time", "range", "string_length" ]
        var_name = ["time", "range", "azimuth", "elevation" ]
        sweep_fixed_angle = np.empty([self.nsweep])
        
        for n in range(self.nsweep) :
            src = netCDF4.Dataset(self.files[n], "r", format="NETCDF4")     
            # create group of each sweep
            gr = nc.createGroup("sweep_"+str(n).zfill(4) )
            # treat sweep related variables
            name = "sweep_number"
            gr.createVariable( name,  dtype("int32").char, () )
            gr[name].setncatts(src[name].__dict__)
            gr.variables[name][:] = src.variables[name][0]         
            name = "sweep_mode"
            gr.createVariable(name, "S1", ("string_length"))
            gr[name].setncatts(src[name].__dict__)
            gr.variables[name][:] = src.variables[name][0]         
            name = "fixed_angle" 
            gr.createVariable(name, dtype("float32").char, ())            
            gr[name].setncatts(src[name].__dict__)
            gr.variables[name][:] = src.variables[name][0]
            sweep_fixed_angle[n]  = src.variables[name][0]  
            # 
            for name, dimension in src.dimensions.items():
                if name in dimension_name :
                    gr.createDimension(
                        name, (len(dimension) if not dimension.isunlimited() else None))           
            for name, variable in src.variables.items():
                if name in var_name :
                    x = gr.createVariable(name, variable.datatype, variable.dimensions)
                    gr.variables[name][:] = src.variables[name][:]          
                    gr[name].setncatts(src[name].__dict__)

                if variable.dimensions == ('time', 'range') : 
                    x = gr.createVariable(name, variable.datatype, variable.dimensions)
                    gr[name].setncatts(src[name].__dict__)                                        
                    gr.variables[name][:] = src.variables[name][:]          
            src.close()
        # write sweep fixed angle
        x = nc.createVariable("sweep_fixed_angle",dtype("float32").char,("sweep"))            
        x.long_name = "fixed_angle_for_sweep"
        x.units = "degrees"
        x[:] = sweep_fixed_angle[:]
        nc.close() 

if __name__ == "__main__":
    main()


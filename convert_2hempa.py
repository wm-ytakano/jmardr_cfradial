import os
import subprocess
from datetime import datetime
from datetime import timedelta
import numpy as np
from numpy import dtype
import netCDF4
from netCDF4 import stringtochar
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("grib")
    parser.add_argument("netcdf")
    parser.add_argument("-w", "--work", help="path to working directory")
    args = parser.parse_args()
    converter = Converter()
    converter.convert(args)


class Converter:
    _FillValueF32 = np.array(9.999e20, "float32")
    _FillValueF64 = np.array(9.999e20, "float64")
    _FillValueI32 = np.array(-9999, "int32")

    def convert(self, args):
        self.gribpath = args.grib
        self.ncpath = args.netcdf
        if args.work:
            self.workdir = args.work
        else:
            self.workdir = "."
        self.read_grib()
        self.write_netcdf()

    def read_grib(self):
        self.read_header()
        self.read_values()

    def read_header(self):
        with open(self.gribpath, "rb") as f:
            data = bytearray(f.read())

        offset = 0

        # Section 0
        offset += 16

        # Section 1
        year = read_int(data, 13, 14, offset)
        month = read_int(data, 15, 15, offset)
        day = read_int(data, 16, 16, offset)
        hour = read_int(data, 17, 17, offset)
        minute = read_int(data, 18, 18, offset)
        second = read_int(data, 19, 19, offset)
        tstr = f"{year:0>4}-{month:0>2}-{day:0>2}T{hour:0>2}:{minute:0>2}:{second:0>2}"
        self.time_reference = datetime.strptime(tstr, "%Y-%m-%dT%H:%M:%S")
        offset += read_int(data, 1, 4, offset)

        # The time(time) coordinate variable stores the time of each ray, in seconds, from a reference time,
        # which is normally the start of the volume (time_coverage_start)
        # but may be a specified reference time (time_reference)
        self.time = []

        # The elevation(time) coordinate variable stores the elevation angle for each ray.
        self.elevation = []

        # The azimuth(time) coordinate variable stores the azimuth angle for each ray
        self.azimuth = []

        # The number of the sweep, in the volume scan. 0-based.
        self.sweep_number = []

        # Options are: "sector", "coplane", rhi", "vertical_pointing", "idle", "azimuth_surveillance", "elevation_surveillance", "sunscan", "pointing", "manual_ppi", "manual_rhi"
        self.sweep_mode = []

        # Target angle for the sweep. elevation in most modes. azimuth in RHI mode.
        self.fixed_angle = []

        # Index of first ray in sweep, relative to start of volume. 0-based
        self.sweep_start_ray_index = []

        # Index of last ray in sweep, relative to start of volume. 0-based
        self.sweep_end_ray_index = []

        # nbの最大値
        self.max_nb = 0
        self.nb_list = []
        self.nr_list = []
        self.nb_list_all = []
        self.nr_list_all = []
        sweep_index = 0

        while True:
            # Section 8
            if data[offset:offset + 4] == b"7777":
                break

            # Section 3
            if data[offset + 4] == 3:
                template_no = read_int(data, 13, 14, offset)
                if template_no != 50120:
                    raise GRIBDecodeError(
                        f"template 3.{template_no}には対応していません")
                nb = read_int(data, 15, 18, offset)  # 径線に沿った資料ビン(data bins)の数
                nr = read_int(data, 19, 22, offset)  # 径線の数
                dx = read_int(data, 31, 34, offset) * 1e-3  # 径線に沿ったビンの間隔
                dstart = read_int(data, 35, 38, offset)
                radar_range = dstart + np.arange(nb) * dx + dx / 2
                azi = read_int(data, 40, 41, offset) * 1e-2  # 開始方位
                if nb > self.max_nb:
                    self.max_nb = nb
                    self.radar_range = radar_range

            # Section 4
            if data[offset + 4] == 4:
                template_no = read_int(data, 8, 9, offset)
                if template_no != 51022:
                    raise GRIBDecodeError(
                        f"template 4.{template_no}には対応していません")
                self.nb_list_all.append(nb)
                self.nr_list_all.append(nr)
                self.parameter_number = read_int(data, 11, 11, offset)
                self.latitude = read_int(data, 15, 18, offset) * 1e-6
                self.longitude = read_int(data, 19, 22, offset) * 1e-6
                self.altitude = read_int(data, 23, 24, offset) * 1e-1
                self.site_id = read_int(data, 29, 30, offset)
                self.frequency = read_int(data, 33, 36, offset) * 1e3
                time_start = read_int_sgn(data, 51, 52, offset)
                time_end = read_int_sgn(data, 53, 54, offset)

                self.time.extend(np.linspace(time_start, time_end, nr))
                self.nb_list.append(nb)
                self.nr_list.append(nr)

                self.sweep_start_ray_index.append(len(self.azimuth))
                azimuth = [(azi + 360 * i / nr) % 360 for i in range(nr)]
                self.azimuth.extend(azimuth)
                self.sweep_end_ray_index.append(len(self.azimuth) - 1)

                for x in range(nr):
                    elv = read_int_sgn(
                        data, 61 + 4 * x, 62 + 4 * x, offset)
                    self.elevation.append(elv * 1e-2)

                self.sweep_number.append(len(self.sweep_number))
                fixed_angle = read_int_sgn(data, 42, 43, offset) * 1e-2
                self.fixed_angle.append(fixed_angle)
                if fixed_angle < 90:
                    self.sweep_mode.append("azimuth_surveillance")
                else:
                    self.sweep_mode.append("vertical_pointing")
                sweep_index += 1

            offset += read_int(data, 1, 4, offset)

    def read_values(self):
        gribname = os.path.basename(self.gribpath)
        grdname = f"{self.workdir}/{gribname}.grd"
        subprocess.run(["wgrib2", self.gribpath, "-no_header", "-bin", grdname],
                       stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        buffer = np.fromfile(grdname, dtype="<f")
        subprocess.run(["rm", "-f", grdname],
                       stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)

        self.data = np.full((np.sum(self.nr_list), self.max_nb),
                            self._FillValueF32,
                            dtype="float32")
        br_offset = 0
        r_offset = 0
        for sweep_index, (nb, nr) in enumerate(zip(self.nb_list_all, self.nr_list_all)):
            bstr = br_offset
            bend = br_offset + nb * nr
            dstr = r_offset
            dend = r_offset + nr
            self.data[dstr:dend, 0:nb] = buffer[bstr:bend].reshape(
                (nr, nb))
            r_offset += nr
            br_offset += nb * nr

    def write_netcdf(self):
        self.nc = netCDF4.Dataset(self.ncpath, "w", format="NETCDF4")
        self.write_global_attributes()  # Section 4.1
        self.write_dimensions()  # Section 4.2
        self.write_global_variables()  # Section 4.3
        self.write_coordinate_variables()  # Section 4.4
        # Section 4.5 Ray dimension variables: ommited
        self.write_location_variables()  # Section 4.6
        self.write_sweep_variables()  # Section 4.7
        self.write_sensor_pointing_variables()  # Section 4.8
        # Section 4.9 Moving platform geo-reference variables: omitted
        self.write_moments_field_data_variables()  # Section 4.10
        self.write_instrument_parameters()  # Section 5.1
        self.nc.close()

    def write_global_attributes(self):
        nc = self.nc

        # Conventions string will specify CF/Radial, plus selected sub-conventions as applicable
        nc.setncattr("Conventions", "CF/Radial instrument_parameters")

        # [optional] CF/Radial version number
        nc.setncattr("version", "1.3")

        # Short description of file contents
        nc.setncattr("title", "")

        # Where the original data were produced
        nc.setncattr("institution", "Japan Meteorological Agency")

        # Method of production of the original data
        nc.setncattr("source", "")

        # List of modifications to the original data
        nc.setncattr("history", "")

        # Miscellaneous information
        nc.setncattr("comment", "")

        # Name of radar or lidar
        nc.setncattr("instrument_name", "")

        # [optional] Name of site where data were gathered
        nc.setncattr("site_name", str(self.site_id))

        # [optional] Name of scan strategy used, if applicable
        nc.setncattr("scan_name", "")

        # [optional] Scan strategy id, if applicable. Assumed 0 if missing
        nc.setncattr("scan_id", "0")

        # [optional] "true" or "false" Assumed "false" if missing.
        nc.setncattr("platform_is_mobile", "false")

        # [optional] "true" or "false" Assumed "false" if missing.
        nc.setncattr("n_gates_vary", "false")

        # [optional] "true" or "false" Assumed "true”" if missing. This is set to false if the rays are not stored in time order.
        nc.setncattr("ray_times_increase", "true")

        # [optional] Comma-delimited list of field names included in this file.
        nc.setncattr("field_names", "")

    def write_dimensions(self):
        nc = self.nc

        # The number of rays. This dimension is optionally UNLIMITED
        nc.createDimension("time", len(self.time))

        # The number of range bin
        nc.createDimension("range", self.max_nb)

        # The number of sweeps
        nc.createDimension("sweep", len(self.sweep_number))

        # [optional] Number of frequencies used
        nc.createDimension("frequency", 1)

        nc.createDimension("string_length", None)

    def write_global_variables(self):
        nc = self.nc

        # Volume numbers are sequential, relative to some arbitrary start time, and may wrap.
        volume_number = nc.createVariable("volume_number", dtype("int32").char)
        volume_number[:] = 0  # 暫定的に0を格納
        volume_number.long_name = "data_volume_index_number"
        volume_number.units = "unitless"

        # TC time of first ray in file. Resolution is integer seconds.
        # The time(time) variable is computed relative to this time.
        # Format is: yyyy-mm-ddThh:mm:ssZ
        time_coverage_start = nc.createVariable(
            "time_coverage_start", "S1", ('string_length'))
        t = self.time_reference + timedelta(seconds=self.time[0])
        tstr = t.strftime("%Y-%m-%dT%H:%M:%SZ")
        datain = np.array(tstr, dtype="S20")
        time_coverage_start[:] = stringtochar(datain)
        time_coverage_start.long_name = "data_volume_start_time_utc"
        time_coverage_start.units = "unitless"

        # UTC time of last ray in file. Resolution is integer seconds.
        # Format is: yyyy-mm-ddThh:mm:ssZ
        time_coverage_end = nc.createVariable(
            "time_coverage_end", "S1", ('string_length'))
        t = self.time_reference + timedelta(seconds=self.time[-1])
        tstr = t.strftime("%Y-%m-%dT%H:%M:%SZ")
        datain = np.array(tstr, dtype="S20")
        time_coverage_end[:] = stringtochar(datain)
        time_coverage_end.long_name = "data_volume_end_time_utc"
        time_coverage_end.units = "unitless"

        # UTC time reference. Resolution is integer seconds.
        # If defined, the time(time) variable is computed relative to this time instead of relative to time_coverage_start.
        # Format is: yyyy-mm-ddThh:mm:ssZ
        time_reference = nc.createVariable(
            "time_reference", "S1", ('string_length'))
        tstr = self.time_reference.strftime("%Y-%m-%dT%H:%M:%SZ")
        datain = np.array(tstr, dtype="S20")
        time_reference[:] = stringtochar(datain)
        time_reference.long_name = "time_reference"
        time_reference.units = "unitless"

    def write_coordinate_variables(self):
        nc = self.nc

        # Coordinate variable for time.
        # Time at center of each ray, in fractional seconds since time_coverage_start.
        time = nc.createVariable("time",
                                 dtype("double").char,
                                 ("time"))
        time[:] = np.array(self.time)
        time.long_name = "time_in_seconds_since_volume_start"
        time.units = f"seconds since {self.time_reference.strftime('%Y-%m-%dT%H:%M:%SZ')}"
        time.calendar = "gregorian"

        # Coordinate variable for range. Range to center of each bin.
        radar_range = nc.createVariable("range",
                                        dtype("float32").char,
                                        ("range"))
        radar_range[:] = self.radar_range.astype("float32")
        radar_range.standard_name = "projection_range_coordinate"
        radar_range.long_name = "range_to_measurement_volume"
        radar_range.units = "meters"
        radar_range.spacing_is_constant = "true"
        radar_range.meters_to_center_of_first_gate = radar_range[0]
        radar_range.meters_between_gates = radar_range[1] - radar_range[0]
        radar_range.axis = "radial_range_coordinate"

    def write_location_variables(self):
        nc = self.nc

        # Latitude of instrument.
        # For a stationary platform, this is a scalar.
        # For a moving platform, this is a vector.
        latitude = nc.createVariable("latitude", dtype('double').char)
        latitude[:] = self.latitude
        latitude.long_name = "latitude"
        latitude.units = "degrees_north"

        # Longitude of instrument.
        # For a stationary platform, this is a scalar.
        # For a moving platform, this is a vector.
        longitude = nc.createVariable("longitude", dtype('double').char)
        longitude[:] = self.longitude
        longitude.long_name = "longitude"
        longitude.units = "degrees_east"

        # Altitude of instrument above mean sea level.
        # For a stationary platform, this is a scalar.
        # For a moving platform, this is a vector.
        altitude = nc.createVariable("altitude", dtype('double').char)
        altitude[:] = self.altitude
        altitude.long_name = "altitude"
        altitude.units = "meters"

    def write_sweep_variables(self):
        nc = self.nc

        # The number of the sweep, in the volume scan. 0-based.
        sweep_number = nc.createVariable(
            "sweep_number", dtype("int32").char, ("sweep"))
        sweep_number[:] = np.array(self.sweep_number, dtype="int32")
        sweep_number.long_name = "sweep_index_number_0_based"
        sweep_number.units = "unitless"

        # Options are: "sector", "coplane", rhi", "vertical_pointing", "idle", "azimuth_surveillance", "elevation_surveillance", "sunscan", "pointing", "manual_ppi", "manual_rhi"
        sweep_mode = nc.createVariable(
            "sweep_mode", "S1", ("sweep", "string_length"))
        datain = np.array(self.sweep_mode, dtype='S22')
        sweep_mode[:] = stringtochar(datain)
        sweep_mode.long_name = "scan_mode_for_sweep"
        sweep_mode.unit = "unitless"

        # Target angle for the sweep. elevation in most modes. azimuth in RHI mode.
        fixed_angle = nc.createVariable(
            "fixed_angle", dtype("float32").char, ("sweep"))
        fixed_angle[:] = np.array(self.fixed_angle, dtype="float32")
        fixed_angle.long_name = "target_fixed_angle"
        fixed_angle.units = "degrees"

        # Index of first ray in sweep, relative to start of volume. 0-based
        sweep_start_ray_index = nc.createVariable(
            "sweep_start_ray_index", dtype("int32").char, ("sweep"))
        sweep_start_ray_index[:] = np.array(
            self.sweep_start_ray_index, dtype="int32")
        sweep_start_ray_index.long_name = "index_of_first_ray_in_sweep"
        sweep_start_ray_index.units = "unitless"

        # Index of last ray in sweep, relative to start of volume. 0-based
        sweep_end_ray_index = nc.createVariable(
            "sweep_end_ray_index", dtype("int32").char, ("sweep"))
        sweep_end_ray_index[:] = np.array(
            self.sweep_end_ray_index, dtype="int32")
        sweep_end_ray_index.long_name = "index_of_last_ray_in_sweep"
        sweep_end_ray_index.units = "unitless"

    def write_sensor_pointing_variables(self):
        nc = self.nc

        # Azimuth of antenna, relative to true north.
        azimuth = nc.createVariable("azimuth", dtype("float32").char, ("time"))
        azimuth[:] = np.array(self.azimuth, dtype="float32")
        azimuth.standard_name = "ray_azimuth_angle"
        azimuth.long_name = "azimuth_angle_from_true_north"
        azimuth.units = "degrees"
        azimuth.axis = "radial_azimuth_coordinate"

        # Elevation of antenna, relative to the horizontal plane.
        elevation = nc.createVariable(
            "elevation", dtype("float32").char, ("time"))
        elevation[:] = np.array(self.elevation, dtype="float32")
        elevation.standard_name = "ray_elevation_angle"
        elevation.long_name = "elevation_angle_from_horizontal_plane"
        elevation.units = "degrees"
        elevation.axis = "radial_elevation_coordinate"

    def write_moments_field_data_variables(self):
        nc = self.nc
        if self.parameter_number == 1:
            short_name = "DBZ"
            standard_name = "equivalent_reflectivity_factor"
            units = "dBZ"
        elif self.parameter_number == 2:
            short_name = "VEL"
            standard_name = "radial_velocity_of_scatterers_away_from_instrument"
            units = "m/s"
        variable = nc.createVariable(short_name,
                                     dtype("float32").char,
                                     ("time", "range"),
                                     fill_value=self._FillValueF32)
        variable[:] = self.data
        variable.standard_name = standard_name
        variable.units = units

    def write_instrument_parameters(self):
        nc = self.nc

        # List of operating frequencies, in Hertz. In most cases, only a single frequency is used.
        frequency = nc.createVariable(
            "frequency", dtype("float32").char, ("frequency"))
        frequency[:] = np.array([self.frequency], dtype="float32")
        frequency.long_name = "radiation_frequency"
        frequency.units = "s-1"


def read_int(ary, i0, i1, offset):
    return int.from_bytes(ary[offset + i0 - 1:offset + i1], "big")


def read_int_sgn(ary, i0, i1, offset):
    l = 8 * (i1 - i0 + 1)
    result = read_int(ary, i0, i1, offset)
    if result & 1 << (l - 1) > 0:
        result = -(result & ~(1 << (l - 1)))
    return result


class GRIBDecodeError(Exception):
    pass


if __name__ == "__main__":
    main()

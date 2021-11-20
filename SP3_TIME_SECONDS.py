import wget
from ftplib import FTP, error_perm
import gzip
import shutil
# from unlzw import unlzw
import zlib
import unlzw3
from pathlib import Path
import datetime as dt
import os
import numpy as np
import pandas as pd
import pymap3d as pm

password = 'thakrarc@my.erau.edu'
username = 'anonymous'

ftp = FTP('garner.ucsd.edu')
ftp.login(user=username, passwd=password)


def projalt(site, az, el, proj_alt=300.):
    print(site)
    lat0, lon0, alt0 = site
    az = az * np.pi / 180.
    el = el * np.pi / 180.

    x, y, z = pm.geodetic2ecef(lat0, lon0, alt0)
    vx, vy, vz = pm.enu2uvw(np.cos(el) * np.sin(az), np.cos(el) * np.cos(az), np.sin(el), lat0, lon0)

    earth = pm.Ellipsoid()
    a2 = (earth.semimajor_axis + proj_alt * 1000.) ** 2
    b2 = (earth.semimajor_axis + proj_alt * 1000.) ** 2
    c2 = (earth.semiminor_axis + proj_alt * 1000.) ** 2

    A = vx ** 2 / a2 + vy ** 2 / b2 + vz ** 2 / c2
    B = x * vx / a2 + y * vy / b2 + z * vz / c2
    C = x ** 2 / a2 + y ** 2 / b2 + z ** 2 / c2 - 1

    alpha = (np.sqrt(B ** 2 - A * C) - B) / A

    lat, lon, alt = pm.ecef2geodetic(x + alpha * vx, y + alpha * vy, z + alpha * vz)

    return lat, lon, alt / 1000.


def parse_sp3(filename, output_root=r"E:\SP3"):
    sat_ephem = {}
    time = []

    with open(filename, 'r') as sp3file:
        lines = sp3file.readlines()
        num_sat = int(lines[2].split()[1])
        sat_lin = lines[2]
        sat_lin2 = lines[3]

        print("num of sat ", num_sat)
        for s in range(32):
            sat = '{:02}'.format(s + 1)
            if (sat in sat_lin) or (sat in sat_lin2):
                sat_ephem[sat] = []
        # print(sat_ephem.keys())
        for j in range(22, len(lines) - 1, num_sat + 1):
            # print(j)
            t = lines[j][2:-6]
            # print(t)
            time.append(dt.datetime.strptime(t, ' %Y %m %d %H %M %S.%f'))

            for i in range(num_sat):
                s = lines[j + i + 1].split()
                # print(lines[j+i+1])
                sat_ephem[s[0][2:]].append([float(s[1]), float(s[2]), float(s[3])])
    d = {}
    for ind, sat in enumerate(sat_ephem):

        temp = np.array(sat_ephem[sat]).T * 1000.
        if ind == 0:

            # array_com = np.column_stack((np.array(time),np.array(sat_ephem[sat])*1000.))
            time_df = np.array(time).T

            d["Time"] = pd.DataFrame(columns=['Time'],
                                     data=np.array(time))
            d["G{0}".format(sat)] = pd.DataFrame(columns=['x', 'y', 'z'],
                                                 data=np.array(sat_ephem[sat]) * 1000.)
        else:
            d["G{0}".format(sat)] = pd.DataFrame(columns=['x', 'y', 'z'],
                                                 data=np.array(sat_ephem[sat]) * 1000.)

        final_df = pd.concat(d, axis=1)
        sat_ephem[sat] = temp
    final_df.to_csv(os.path.join(output_root, filename.replace("sp3", "csv")), index=False)
    return sat_ephem, np.array(time), final_df


def download_sp3(date, output_root=r"E:\SP3"):
    gps_start = dt.datetime(1980, 1, 6)
    day_of_week = (date.weekday() + 1) % 7
    week_number = (date - gps_start).days // 7
    filename = "igs{0}{1}.sp3.Z".format(week_number, day_of_week)
    output_dir = os.path.join(output_root, date.strftime("%Y%m%d"))
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    file_output = os.path.join(output_dir, filename)
    remotefile = "pub/products/{0}/{1}".format(week_number, filename)
    try:
        with open(file_output, 'wb') as fp:
            ftp.retrbinary('RETR {}'.format(remotefile), fp.write)
    except:
        os.remove(output_dir)
        return None
    uncompressed_data = unlzw3.unlzw(Path(file_output))

    with open(str(file_output)[:-2], 'wb') as f_out:
        # uncompressed_data = zlib.decompress(compressed_data)
        f_out.write(uncompressed_data)
    return file_output[:-2]


######################################################################User Inputs#############################
sit = [29.20357, -81.0394, 100]  # Location of receiver. Currently set for DAB
satlite = 2  # PRN number
starttime = dt.datetime(2021, 6, 18, 16, 30)  # date, and start time. Format year,month, day, hour, minute
endtime = dt.datetime(2021, 6, 18, 17, 30)  # date, and end time. Format year,month, day, hour, minute
output_root = r"E:\SP3"  # root directroy where data will be stored.
height = 120  # Height of IPP
######################################################################User Inputs#############################

sp3file = download_sp3(starttime, output_root)
ecef, time, data_df = parse_sp3(sp3file, output_root)

satlite = (2 - len(str(satlite))) * '0' + str(satlite)

ind = np.where((time >= starttime) & (time <= endtime))
x = ecef[satlite][0][ind]
y = ecef[satlite][1][ind]
z = ecef[satlite][2][ind]
time = time[ind]

delta = int((time[1] - time[0]).seconds / 60)
# print(delta)

az, el, srange = pm.ecef2aer(x, y, z, sit[0], sit[1], sit[2])

data_axelv = pd.DataFrame({"Time": time, "Azimuth": az, "Elevation": el, "PRN": satlite})
data_axelv.to_csv(sp3file.replace(".sp3", "azelv.csv"), index=False)

print(len(az))
# print((len(az)-1)*delta)
# new_len = (len(az)-1)*delta
# az = np.linspace(az[0],az[-1],new_len)
# el = np.linspace(el[0],el[-1],new_len)
# print(time)
print(el)
prev = time[0]
prev_ind = 0
time_new = []
el_new = np.array([])
az_new = np.array([])
for index, time_sing in enumerate(time[1:]):

    minute = int((time_sing - prev).seconds / 60)

    if (time_sing == time[-1]):
        time_new = time_new + [prev + dt.timedelta(minutes=i) for i in range(minute + 1)]
        el_new = np.append(el_new, np.linspace(el[prev_ind], el[index + 1], minute + 1))
        az_new = np.append(az_new, np.linspace(az[prev_ind], az[index + 1], minute + 1))
    else:
        time_new = time_new + [prev + dt.timedelta(minutes=i) for i in range(minute)]
        el_new = np.append(el_new, np.linspace(el[prev_ind], el[index + 1], minute + 1)[:-1])
        az_new = np.append(az_new, np.linspace(az[prev_ind], az[index + 1], minute + 1)[:-1])

    prev_ind = index + 1

    prev = time_sing
    print("______________")
data_axelv_lin = pd.DataFrame({"Time": time_new, "Azimuth": az_new, "Elevation": el_new, "PRN": satlite})
data_axelv_lin.to_csv(sp3file.replace(".sp3", "azelv_linspace.csv"), index=False)

lat_ipp, long_ipp, alt_ipp = projalt(sit, az_new, el_new, height)

data_ipp = pd.DataFrame({"Time": time_new, "lat": lat_ipp, "long": long_ipp, "PRN": satlite})
data_ipp.to_csv(sp3file.replace(".sp3", "_IPP_{0}.csv").format(height), index=False)

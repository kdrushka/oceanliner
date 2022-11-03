# Native packages
from math import radians, degrees, sin, cos, asin, acos, sqrt
import datetime
import sys
import os
import requests
from pathlib import Path

# Third-party packages for data manipulation
import numpy as np
import pandas as pd
import xarray as xr

# Third-party packages for data interpolation
from scipy import interpolate
from xgcm import Grid


from netrc import netrc
from urllib import request
from platform import system
from getpass import getpass
from http.cookiejar import CookieJar
from os.path import expanduser, join
from datetime import datetime, date, time, timedelta
import gsw as sw
import numpy as np
import xgcm.grid
import netCDF4 as nc4




# ***This library includes*** 
# - setup_earthdata_login_auth
# - download_llc4320_data
# - get_survey_track
# - survey_interp
# - great_circle
# - compute_derived_fields


# NOTE: Code to download data from PO.DAAC adapted from this demo : https://github.com/podaac/tutorials/blob/master/notebooks/Pre-SWOT_Numerical_Simulation_Demo.ipynb by Jinbo Wang & Jack McNelis (JPL) 



def setup_earthdata_login_auth(endpoint: str='urs.earthdata.nasa.gov'):
    """Set up the the Earthdata login authorization for downloading data.

    Extended description of function.

    Returns:
        bool: True if succesful, False otherwise.
        
    Raises: 
        FileNotFoundError: If the Earthdata account details entered are incorrect.
    

    """
    return True
    netrc_name = "_netrc" if system()=="Windows" else ".netrc"
    try:
        username, _, password = netrc(file=join(expanduser('~'), netrc_name)).authenticators(endpoint)
    except (FileNotFoundError, TypeError):
        print('Please provide your Earthdata Login credentials for access.')
        print('Your info will only be passed to %s and will not be exposed in Jupyter.' % (endpoint))
        username = input('Username: ')
        password = getpass('Password: ')
    manager = request.HTTPPasswordMgrWithDefaultRealm()
    manager.add_password(None, endpoint, username, password)
    auth = request.HTTPBasicAuthHandler(manager)
    jar = CookieJar()
    processor = request.HTTPCookieProcessor(jar)
    opener = request.build_opener(auth, processor)
    request.install_opener(opener)
    

def rotate_vector_to_EN(U, V, AngleCS, AngleSN):
    """Rotate vector to east north direction.
    
    Assumes that AngleCS and AngleSN are already of same dimension as V and U (i.e. already interpolated to cell center)
                
    Args:
        U (xarray Dataarray): Zonal vector component
        V (array Dataarray): Meridonal vector component

    Returns:
        uE (xarray Dataarray): Rotated zonal component
        vN (xarray Dataarray): Rotated meridonial component
        
    Raises: 
        FileNotFoundError: If the Earthdata account details entered are incorrect.
    
    Note: adapted from https://github.com/AaronDavidSchneider/cubedsphere/blob/main/cubedsphere/regrid.py
    

    """
               
    # rotate the vectors:
    uE = AngleCS * U - AngleSN * V
    vN = AngleSN * U + AngleCS * V

    return uE, vN
            

def download_llc4320_data(RegionName, datadir, start_date, ndays):
    """Download the MITgcm LLC4320 data from PODAAC Earthdata website.

    Creates a http access for each target file using the setup_earthdata_login_auth function. Checks for existing llc4320 files in 'datadir' and downloads them in the datadir if not found.

    Args:
        RegionName (str): Can be selected from the list of regions with pre-SWOT llc4320 data
        datadir (str): Directory where input model files are stored
        start_date (datetime): Starting date for downloading data
        ndays (int): Number of days to be downloaded from the start date

    Returns:
        None
        
    Raises: 
        FileNotFoundError: If the Earthdata account details entered are incorrect
        error-skipping this file: If the file already exists
    

    """
   
    ShortName = "MITgcm_LLC4320_Pre-SWOT_JPL_L4_" + RegionName + "_v1.0"
    date_list = [start_date + timedelta(days=x) for x in range(ndays)]
    target_files = [f'LLC4320_pre-SWOT_{RegionName}_{date_list[n].strftime("%Y%m%d")}.nc' for n in range(ndays)] # list of files to check for/download
    setup_earthdata_login_auth()
    
    # https access for each target_file
    url = "https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-protected"
    https_accesses = [f"{url}/{ShortName}/{target_file}" for target_file in target_files]
   
    # create datadir if it doesn't exist
    Path(datadir).mkdir(parents=True, exist_ok=True) 

    # list of dataset objects
    dds = []
    for https_access,target_file in zip(https_accesses,target_files):        
        # if no file locally, download it
        if not(os.path.isfile(datadir + target_file)):
            print('downloading ' + target_file) # print file name
            try:
                filename_dir = os.path.join(datadir, target_file)
                request.urlretrieve(https_access, filename_dir)
            except:
                print(' ---- error - skipping this file')



def compute_derived_fields(RegionName, datadir, start_date, ndays, DERIVED_VARIABLES):
    """ Check for derived files in {datadir}/derived and compute if the files don't exist


    Args:
        RegionName (str): It can be selected from the list of regions with pre-SWOT llc4320 data
        datadir (str): Directory where input model files are stored
        start_date (datetime): Starting date for computing derived fields
        ndays (int): Number of days from the start date to compute derived fields
        DERIVED_VARIABLES (str list): specifies which variables to derive (steric_height and/or vorticity)

    Returns:
        None
        
    Raises: 
        TBD: TBD

    """
    
    # directory to save derived data to - create if doesn't exist
    derivedir = datadir + 'derived/'
    if not(os.path.isdir(derivedir)):
        os.mkdir(derivedir)
        
    # files to load:
    date_list = [start_date + timedelta(days=x) for x in range(ndays)]
    target_files = [f'{datadir}LLC4320_pre-SWOT_{RegionName}_{date_list[n].strftime("%Y%m%d")}.nc' for n in range(ndays)] # list target files
    
    # list of derived files:
    derived_files = [f'{derivedir}LLC4320_pre-SWOT_{RegionName}_derived-fields_{date_list[n].strftime("%Y%m%d")}.nc' for n in range(ndays)] # list target files

        
    # loop through input files, then do the following:
    # - compute steric height
    # - interpolate vector quantities (velocity and wind) to the tracer grid
    # - compute vorticity
    fis = range(len(target_files))
    
    cnt = 0 # count
    for fi in fis:
        # input filename:
        thisf=target_files[fi]
        # output filename:
        fnout = thisf.replace(RegionName + '_' , RegionName + '_derived-fields_')
        fnout = fnout.replace(RegionName + '/' , RegionName + '/derived/')
        # check if output file already exists
        if (not(os.path.isfile(fnout))):   
            print(f'computing {DERIVED_VARIABLES} for {thisf}') 
            # load file:
            ds = xr.open_dataset(thisf)
            
            if 'steric_height' in DERIVED_VARIABLES:
                # -------
                # first time through the loop, load reference profile:
                # load a single file to get coordinates
                if cnt==0:
                    # mean lat/lon of domain
                    xav = ds.XC.isel(j=0).mean(dim='i')
                    yav = ds.YC.isel(i=0).mean(dim='j')

                    # for transforming U and V, and for the vorticity calculation, build the xgcm grid:
                    # see https://xgcm.readthedocs.io/en/latest/xgcm-examples/02_mitgcm.html
                    grid = xgcm.Grid(ds, coords={'X':{'center': 'i', 'left': 'i_g'}, 
                                 'Y':{'center': 'j', 'left': 'j_g'},
                                 'T':{'center': 'time'},
                                 'Z':{'center': 'k'}})


                    # --- load reference file of argo data
                    # here we use the 3x3 annual mean Argo product on standard produced by IRPC & distributed by ERDDAP
                    # https://apdrc.soest.hawaii.edu/erddap/griddap/hawaii_soest_defb_b79c_cb17.html
                    # - download the profile closest to xav,yav once (quick), use it, then delete it.

                    # URL gets temp & salt at all levels
                    argofile = f'https://apdrc.soest.hawaii.edu/erddap/griddap/hawaii_soest_625d_3b64_cc4d.nc?temp[(0000-12-15T00:00:00Z):1:(0000-12-15T00:00:00Z)][(0.0):1:(2000.0)][({yav.data}):1:({yav.data})][({xav.data}):1:({xav.data})],salt[(0000-12-15T00:00:00Z):1:(0000-12-15T00:00:00Z)][(0.0):1:(2000.0)][({yav.data}):1:({yav.data})][({xav.data}):1:({xav.data})]'

                    # delete the argo file if it exists 
                    if os.path.isfile('argo_local.nc'):
                        os.remove('argo_local.nc')
                    # use requests to get the file, and write locally:
                    r = requests.get(argofile)
                    file = open('argo_local.nc','wb')
                    file.write(r.content)
                    file.close()
                    # open the argo file:
                    argods = xr.open_dataset('argo_local.nc',decode_times=False)
                    # get rid of time coord/dim/variable, which screws up the time in ds if it's loaded
                    argods = argods.squeeze().reset_coords(names = {'time'}, drop=True) 
                    # reference profiles: annual average Argo T/S using nearest neighbor
                    Tref = argods["temp"]
                    Sref = argods["salt"]
                    # SA and CT from gsw:
                    # see example from https://discourse.pangeo.io/t/wrapped-for-dask-teos-10-gibbs-seawater-gsw-oceanographic-toolbox/466
                    Pref = xr.apply_ufunc(sw.p_from_z, -argods.LEV, yav)
                    Pref.compute()
                    SAref = xr.apply_ufunc(sw.SA_from_SP, Sref, Pref, xav, yav,
                                           dask='parallelized', output_dtypes=[Sref.dtype])
                    SAref.compute()
                    CTref = xr.apply_ufunc(sw.CT_from_pt, Sref, Tref, # Theta is potential temperature
                                           dask='parallelized', output_dtypes=[Sref.dtype])
                    CTref.compute()
                    Dref = xr.apply_ufunc(sw.density.rho, SAref, CTref, Pref,
                                        dask='parallelized', output_dtypes=[Sref.dtype])
                    Dref.compute()


                    cnt = cnt+1
                    print()

                # -------
                # --- COMPUTE STERIC HEIGHT IN STEPS ---
                # 0. create datasets for variables of interest:
                ss = ds.Salt
                tt = ds.Theta
                pp = xr.DataArray(sw.p_from_z(ds.Z,ds.YC))

                # 1. compute absolute salinity and conservative temperature
                sa = xr.apply_ufunc(sw.SA_from_SP, ss, pp, xav, yav, dask='parallelized', output_dtypes=[ss.dtype])
                sa.compute()
                ct = xr.apply_ufunc(sw.CT_from_pt, sa, tt, dask='parallelized', output_dtypes=[ss.dtype])
                ct.compute()
                dd = xr.apply_ufunc(sw.density.rho, sa, ct, pp, dask='parallelized', output_dtypes=[ss.dtype])
                dd.compute()
                # 2. compute specific volume anomaly: gsw.density.specvol_anom_standard(SA, CT, p)
                sva = xr.apply_ufunc(sw.density.specvol_anom_standard, sa, ct, pp, dask='parallelized', output_dtypes=[ss.dtype])
                sva.compute()
                # 3. compute steric height = integral(0:z1) of Dref(z)*sva(z)*dz(z)
                # - first, interpolate Dref to the model pressure levels
                Drefi = Dref.interp(LEV=-ds.Z)
                dz = -ds.Z_bnds.diff(dim='nb').drop_vars('nb').squeeze() # distance between interfaces

                # steric height computation (summation/integral)
                # - increase the size of Drefi and dz to match the size of sva
                Db = Drefi.broadcast_like(sva)
                dzb = dz.broadcast_like(sva)
                dum = Db * sva * dzb
                sh = dum.cumsum(dim='k') 
                # this gives sh as a 3-d variable, (where the depth dimension 
                # represents the deepest level from which the specific volume anomaly was interpolated)
                # - but in reality we just want the SH that was determined by integrating over
                # the full survey depth, which gives a 2-d output:
                sh_true = dum.sum(dim='k') 
            
                # make into dataset:
                sh_ds = sh.to_dataset(name='steric_height')
                sh_true_ds = sh_true.to_dataset(name='steric_height_true')            
                # add/rename the Argo reference profile variables to dout:
                tref = Tref.to_dataset(name='Tref')
                tref = tref.merge(Sref).rename({'salt': 'Sref'}).\
                    rename({'LEV':'zref','latitude':'yav','longitude':'xav'})
            
            if 'vorticity' in DERIVED_VARIABLES:                
                # --- COMPUTE VORTICITY using xgcm and interpolate to X, Y
                # see https://xgcm.readthedocs.io/en/latest/xgcm-examples/02_mitgcm.html
                vorticity = (grid.diff(ds.V*ds.DXG, 'X') - grid.diff(ds.U*ds.DYG, 'Y'))/ds.RAZ
                vorticity = grid.interp(grid.interp(vorticity, 'X', boundary='extend'), 'Y', boundary='extend')
                # make into dataset
                v_ds =vorticity.to_dataset(name='vorticity')
     
    
    
            
            # --- save derived fields in a new file
            if 'steric_height' in DERIVED_VARIABLES:
                dout = sh_ds
                dout = dout.merge(sh_true_ds)
                # - add ref profiles to dout and drop uneeded vars/coords
                dout = dout.merge(tref).drop_vars({'longitude','latitude','LEV'})
                # - add attributes for all variables
                dout.steric_height.attrs = {'long_name' : 'Steric height',
                                        'units' : 'm',
                                        'comments_1' : 'Computed by integrating the specific volume anomaly (SVA) multiplied by a reference density, where the reference density profile is calculated from temperature & salinity profiles from the APDRC 3x3deg gridded Argo climatology product (accessed through ERDDAP). The profile nearest to the center of the domain is selected, and T & S profiles are averaged over one year before computing ref density. SVA is computed from the model T & S profiles. the Gibbs Seawater Toolbox is used compute reference density and SVA. steric_height is given at all depth levels (dep): steric_height at a given depth represents steric height signal generated by the water column above that depth - so the deepest steric_height value represents total steric height (and is saved in steric_height_true'
                                           }
                dout.steric_height_true.attrs = dout.steric_height.attrs
                dout.Tref.attrs = {'long_name' : f'Reference temperature profile at {yav.data}N,{xav.data}E',
                                    'units' : 'degree_C',
                                    'comments_1' : 'From Argo 3x3 climatology produced by APDRC'}
                dout.Sref.attrs = {'long_name' : f'Reference salinity profile at {yav.data}N,{xav.data}E',
                                        'units' : 'psu',
                                        'comments_1' : 'From Argo 3x3 climatology produced by APDRC'}

                dout.zref.attrs = {'long_name' : f'Reference depth for Tref and Sref',
                                        'units' : 'm',
                                        'comments_1' : 'From Argo 3x3 climatology produced by APDRC'}
                
                # merge vorticity 
                if 'vorticity' in DERIVED_VARIABLES:  
                    dout = dout.merge(v_ds)
                    
            # if we only computed vorticity, dout = v_ds
            elif 'vorticity' in DERIVED_VARIABLES:  
                dout = v_ds
                
                
            # if vorticity, add the attrs:
            if 'vorticity' in DERIVED_VARIABLES:  
                dout.vorticity.attrs = {'long_name' : 'Vertical component of the vorticity',
                                    'units' : 's-1',
                                    'comments_1' : 'computed on DXG,DYG then interpolated to X,Y'}         
               
            # - save netcdf file with derived fields
            netcdf_fill_value = nc4.default_fillvals['f4']
            dv_encoding = {}
            for dv in dout.data_vars:
                dv_encoding[dv]={'zlib':True,  # turns compression on\
                            'complevel':1,     # 1 = fastest, lowest compression; 9=slowest, highest compression \
                            'shuffle':True,    # shuffle filter can significantly improve compression ratios, and is on by default \
                            'dtype':'float32',\
                            '_FillValue':netcdf_fill_value}
            # save to a new file
            print(' ... saving to ', fnout)
            # TROUBLESHOOTING::::: DELETE THE RETURN LINE
            #return dout, dv_encoding
            dout.to_netcdf(fnout,format='netcdf4',encoding=dv_encoding)

            
            
    # release Argo file at the end of all files
    if 'argods' in locals():
        argods.close()


    
def set_defaults(sampling_details):
    
    """Calculates the survey indices and track based on the sampling details for the dataset for all days.


    Args:
        ds (xarray.core.dataset.Dataset): MITgcm LLC4320 data for all days
        sampling_details (dict): It includes number of days, waypoints, and depth range, horizontal and vertical platform speed. These can typical (default) or user-specified, in the                                      case where user specfies only some of the details the default values will be used for rest.

    Returns:
        survey_track (xarray.core.dataset.Dataset): Returns the track (lat, lon, depth, time) of the sampling trajectory based on the type of sampling                               
        survey_indices (xarray.core.dataset.Dataset): Returns the indices (i, j, k, time) of the sampling trajectory based on the type of sampling
        sampling_details (dict): Returns the modified sampling_details by filling in the missing parameters with defaults.
        
    Raises: 
        Sampling strategy is invalid: If a sampling strategy is not specified or different from the available strategies - sim_utcd, sim_glider, sim_mooring, wave_glider, sail_drone
    

    """       
    
    # ------ default sampling parameters: in the dict named "defaults" -----
    # these are used when these parameters are not specified by the user
    defaults = {}
    # default values depend on the sampling type
    # typical speeds and depth ranges based on platform 
    SAMPLING_PLATFORM = sampling_details['SAMPLING_PLATFORM']
    if SAMPLING_PLATFORM == 'uctd':
        # typical values for uctd sampling:
        defaults['zrange'] = [-5, -500] # depth range of profiles (down is negative)
        defaults['z_res'] = 1 # the vertical sampling rate in meters
        defaults['hspeed'] = 5 # platform horizontal speed in m/s
        defaults['vspeed'] = 1 # platform vertical (profile) speed in m/s 
        defaults['PATTERN'] = 'lawnmower'
        defaults['AT_END'] = 'terminate'  # behaviour at and of trajectory: 'repeat', 'reverse', or 'terminate'
    elif SAMPLING_PLATFORM == 'glider':
        defaults['zrange'] = [-1, -1000] # depth range of profiles (down is negative)
        defaults['z_res'] = 1 # the vertical sampling rate in meters
        defaults['hspeed'] = 0.25 # platform horizontal speed in m/s
        defaults['vspeed'] = 0.1 # glider vertical (profile) speed in m/s     
        defaults['AT_END'] = 'terminate'  # behaviour at and of trajectory: 'repeat', 'reverse', or 'terminate'
        defaults['PATTERN'] = 'lawnmower'
    elif SAMPLING_PLATFORM == 'wave_glider':
        defaults['zrange'] = [-1, -1.5] # depth range of profiles (down is negative)
        defaults['z_res'] = 0 # the vertical sampling rate in meters - zero for wave glider
        defaults['hspeed'] = 1 # platform horizontal speed in m/s
        defaults['vspeed'] = 0 # surface vehicle, no vertical motion  
        defaults['AT_END'] = 'terminate'  # behaviour at and of trajectory: 'repeat', 'reverse', or 'terminate'
        defaults['PATTERN'] = 'back-forth'
    elif SAMPLING_PLATFORM == 'saildrone':
        defaults['zrange'] = [-1, -3] # depth range of profiles (down is negative)
        defaults['z_res'] = 0 # the vertical sampling rate in meters - zero for saildrone
        defaults['hspeed'] = 2.57 # platform horizontal speed in m/s
        defaults['vspeed'] = 0 # surface vehicle, no vertical motion       
        defaults['AT_END'] = 'terminate'  # behaviour at and of trajectory: 'repeat', 'reverse', or 'terminate'
        defaults['PATTERN'] = 'back-forth'
    elif SAMPLING_PLATFORM == 'mooring':
        defaults['WAYPOINTS'] = {'x':model_xav, 'y':model_yav}  # default lat/lon is the center of the domain
        defaults['zmooring'] = [-1, -10, -50, -100] # depth of T/S/U/V instruments
    else:
        # if SAMPLING_PLATFORM not specified, return an error
        print('error: SAMPLING_PLATFORM ' + SAMPLING_PLATFORM + ' invalid')
        return -1
    
    #
    defaults['SAVE_PRELIMINARY'] = False    
    
    # merge defaults & sampling_details
    # - by putting sampling_details second, items that appear in both dicts are taken from sampling_details: 
    sampling_details = {**defaults, **sampling_details}
    return sampling_details
    
    
def get_survey_track(ds, sampling_details_in):
    
    """Calculates the survey indices and track based on the sampling details for the dataset for all days.


    Args:
        ds (xarray.core.dataset.Dataset): MITgcm LLC4320 data for all days
        sampling_details_in (dict): It includes number of days, waypoints, and depth range, horizontal and vertical platform speed. These can typical (default) or user-specified, in the                                      case where user specfies only some of the details the default values will be used for rest.

    Returns:
        survey_track (xarray.core.dataset.Dataset): Returns the track (lat, lon, depth, time) of the sampling trajectory based on the type of sampling                               
        survey_indices (xarray.core.dataset.Dataset): Returns the indices (i, j, k, time) of the sampling trajectory based on the type of sampling
        sampling_details (dict): Returns the modified sampling_details by filling in the missing parameters with defaults.
        
    Raises: 
        Sampling strategy is invalid: If a sampling strategy is not specified or different from the available strategies - sim_utcd, glider, sim_mooring, wave_glider, saildrone
    

    """
    
    survey_time_total = (ds.time.values.max() - ds.time.values.min()) # (timedelta) - limits the survey to a total time
    survey_end_time = ds.time.isel(time=0).data + survey_time_total # end time of survey
    # Convert lon, lat and z to index i, j and k with f_x, f_y and f_z
    # XC, YC and Z are the same at all times, so select a single time
    X = ds.XC.isel(time=0) 
    Y = ds.YC.isel(time=0)
    i = ds.i
    j = ds.j
    z = ds.Z.isel(time=0)
    k = ds.k
    f_x = interpolate.interp1d(X[0,:].values, i)
    f_y = interpolate.interp1d(Y[:,0].values, j)
    f_z = interpolate.interp1d(z, k, bounds_error=False)

    # Get boundaries and center of model region
    model_boundary_n = Y.max().values
    model_boundary_s = Y.min().values
    model_boundary_w = X.min().values
    model_boundary_e = X.max().values
    model_xav = ds.XC.isel(time=0, j=0).mean(dim='i').values
    model_yav = ds.YC.isel(time=0, i=0).mean(dim='j').values
    
    # --------- define sampling, i.e., the x/y/z/t points to interpolate to
    #  1) call the set_defaults function, which sets details that have not been specified by the user
    #  and stores the correct and complete sampling_details
    sampling_details = set_defaults(sampling_details_in)
    
    # 2) get the waypoints:
    WAYPOINTS = sampling_details['WAYPOINTS']
    if not(WAYPOINTS):
        print('no waypoints specified - using defaults')
        # define some reasonable sampling pattern within the domain based on SAMPLING_PLATFORM and PATTERN and AT_END
        if sampling_details['PATTERN'] == 'lawnmower':
            # "mow the lawn" pattern - define all waypoints
            xwaypoints = model_boundary_w + 1 + [0, 0, 0.5, 0.5, 1, 1, 1.5, 1.5, 2, 2]
            ywaypoints = model_boundary_s + [1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2]
        elif sampling_details['PATTERN'] == 'back-forth':
            # repeated back & forth transects - define the end-points
            xwaypoints = model_xav + [-1, 1]
            ywaypoints = model_yav + [-1, 1]
            # repeat waypoints based on total # of transects: 
            dkm_per_transect = great_circle(xwaypoints[0], ywaypoints[0], xwaypoints[1], ywaypoints[1]) # distance of one transect in km
#           # time per transect, seconds, as a np.timedelta64 value
            t_per_transect = np.timedelta64(int(dkm_per_transect * 1000 / sampling_details['hspeed']), 's')    
            num_transects = np.round(survey_time_total / t_per_transect)
            for n in np.arange(num_transects):
                xwaypoints = np.append(xwaypoints, xwaypoints[-2])
                ywaypoints = np.append(ywaypoints, ywaypoints[-2])
        # if the survey pattern repeats, add the first waypoint to the end of the list of waypoints:
        if sampling_details['AT_END'] == 'repeat': 
            xwaypoints = np.append(xwaypoints, xwaypoints[0])
            ywaypoints = np.append(ywaypoints, ywaypoints[0])         
    elif ['x' in WAYPOINTS]:
        # x and y specified by the user
        xwaypoints = WAYPOINTS['x']
        ywaypoints = WAYPOINTS['y']
    elif ['waypoint_file' in WAYPOINTS]:
        # load file
        waypoints = xr.open_dataset(WAYPOINTS['waypoint_file'])
        xwaypoints = waypoints.xwaypoints.values
        ywaypoints = waypoints.ywaypoints.values
    
               
    
    
    # 3) define the x/y/z/t points to interpolate to
    SAMPLING_PLATFORM = sampling_details['SAMPLING_PLATFORM']
    if SAMPLING_PLATFORM == 'mooring':
        # for moorings, location is fixed and a set of waypoints isn't needed
        ts = ds.time.values # in hours
        # same sampling for T/S/U/V 
        zs = sampling_details['zmooring'] 
        xs = xwaypoints
        ys = ywaypoints
    else:
        # vertical resolution
        # TODO: CHECK z_res=1000........print error!
        zresolution = sampling_details['z_res'] # in meters
        # max depth can't be deeper than the max model depth in this region
        # (note, this code is not aware enough to look at the local depth for each profile or dive)
        sampling_details['zrange'][1] = -np.min([-sampling_details['zrange'][1], ds.Depth.isel(time=1).max(...).values])        
        zprofile = np.arange(sampling_details['zrange'][0],sampling_details['zrange'][1],-zresolution) # depths for one profile
        ztwoway = np.append(zprofile,zprofile[-1::-1])
        # time resolution of sampling (dt):
        dt = zresolution / sampling_details['vspeed'] # sampling resolution in seconds
        dt_td64 = np.timedelta64(int(dt), 's') # np.timedelta64 format
        # for each timestep dt, how far horizontally and vertically does the platform travel:
        deltah = sampling_details['hspeed']*dt # horizontal distance traveled per sample
        deltav = sampling_details['vspeed']*dt # vertical distance traveled per sample

        # determine the sampling locations in 4-d space: xs, ys, zs, ts
        # - initialize sample locations xs, ys, zs, ts
        xs = []
        ys = []
        zs = []
        ts = []
        dkm_total = 0  # total horizontal distance travelled in km    

        # ---- first, horizontal sampling (xs, ys)
        # loop through each waypoint:
        for w in np.arange(len(xwaypoints)-1):
            # ----- interpolate between this and the following waypoint
            # distance between waypoints in km:
            dkm = great_circle(xwaypoints[w], ywaypoints[w], xwaypoints[w+1], ywaypoints[w+1])
            # number of time steps (measurements) between this and the next waypoint
            nstep = int(dkm*1000 / deltah) 
            # interpolated x/y positions between the two waypoints
            yi = np.linspace(ywaypoints[w], ywaypoints[w+1], nstep)
            xi = np.linspace(xwaypoints[w], xwaypoints[w+1], nstep)
            # remove last point, which is the next waypoint
            xi = xi[0:-1] 
            yi = yi[0:-1] 
            # append to xs, ys - the list of sample points
            xs = np.append(xs, xi)
            ys = np.append(ys, yi)
            # update the total distance traveled
            dkm_total = dkm_total + dkm           
            # cumulative survey time to this point, in seconds, as a np.timedelta64 value
            t_total = np.timedelta64(int(dkm_total * 1000 / sampling_details['hspeed']), 's')
            # cut off the survey after survey_time_total
            if t_total > survey_time_total:
                break 
               
        # km for one lap of the survey
        dkm_once = dkm_total 

        # if time is less than survey_time_total, trigger AT_END behavior:
        if t_total < survey_time_total:
            if sampling_details['AT_END'] == 'repeat': 
                # start at the beginning again
                # - determine how many times the survey repeats:
                num_transects = np.round(survey_time_total / t_total)
                x_once = xs
                y_once = ys
                # append x_once and y_once as many times as needed (ie, num_transects)
                for n in np.arange(num_transects):
                    xs = np.append(xs, x_once)
                    ys = np.append(ys, y_once)
                    dkm_total += dkm_once
            elif sampling_details['AT_END'] == 'reverse': 
                # turn around & go in the opposite direction
                # - determine how many times the survey repeats:
                num_transects = np.round(survey_time_total / t_total)
                x_once = xs
                y_once = ys
                # append both a backward & another forward transect, as many times as needed (num_transects/2)
                for n in np.arange(np.ceil(num_transects/2)):
                    xs = np.append(np.append(xs, x_once[-2:1:-1]), x_once)
                    ys = np.append(np.append(ys, y_once[-2:1:-1]), y_once)
                    dkm_total += dkm_once*2

        
        # ---- next, vertical sampling (zs) and time steps (ts)
        # compute the number of vertical profiles we make during the survey:
        n_profiles = np.ceil(xs.size / ztwoway.size)
        # - repeat (tile) the two-way sampling depths 
        zs = np.tile(ztwoway, int(n_profiles))
        zs = zs[0:xs.size] # limit to # of samples
        # timesteps = dt * number of samples (ie, size of xs)
        ts = ds.time.isel(time=0).data + dt_td64 * np.arange(xs.size)
        # get rid of points with sample time > survey_time_total
        if survey_time_total > 0:
            idx = np.argmin(np.abs(ts - survey_end_time))# index of ts closest to survey_end_time
            print('originally, ', idx, ' points')
            # make sure this is multiple of the # of profiles:
            idx = int(np.floor((idx+1)/len(ztwoway)) * (len(ztwoway)))
            xs = xs[:idx]
            ys = ys[:idx]
            ts = ts[:idx]
            zs = zs[:idx]
            n_profiles = np.ceil(xs.size / ztwoway.size)
            # update t_total
            t_total = np.diff(ts[[0,-1]])
            t_total_seconds = int(t_total)/1e9 # convert from nanoseconds to seconds
            # use the speed to determine dkm_total (time * hspeed)
            dkm_total = t_total_seconds * sampling_details['hspeed'] / 1000
            print('limited to ', idx, 'points: n_profiles=', n_profiles, ', ', len(zprofile), 'depths per profile, ', len(ztwoway), 'depths per two-way')
        
        # write total time and distance to sampling_details
        sampling_details['distance_total_km'] = dkm_total
        sampling_details['time_total_s'] = t_total_seconds  
        # -- end if not a mooring
        
    # ----- Assemble the xs/ys/zs/ts sampling locations into "survey_track" 
    # (same regardless of sampling strategy except "mooring")
    if not SAMPLING_PLATFORM == 'mooring':
        # - true (lat/lon) coordinates:
        survey_track = xr.Dataset(
            dict(
                lon = xr.DataArray(xs,dims='points'),
                lat = xr.DataArray(ys,dims='points'),
                dep = xr.DataArray(zs,dims='points'),
                time = xr.DataArray(ts,dims='points'),
                n_profiles = n_profiles
            )
        )
        
    elif SAMPLING_PLATFORM == 'mooring':
        # x/y/z don't change
        survey_track = xr.Dataset(
            dict(
                lon = xr.DataArray(xs*[1], dims='position'),
                lat = xr.DataArray(ys*[1], dims='position'),
                dep = xr.DataArray(zs, dims='depth'),
                time = xr.DataArray(ts, dims='time')
            )
        )
    # Transform survey_track to survey_indices (i,j,k coordinates) with the functions (f_*) defined earlier
    survey_indices= xr.Dataset(
        dict(
            i = xr.DataArray(f_x(survey_track.lon), dims='points'),
            j = xr.DataArray(f_y(survey_track.lat), dims='points'),
            k = xr.DataArray(f_z(survey_track.dep), dims='points'),
            time = xr.DataArray(survey_track.time, dims='points'),
        )
    )
    
    
    
    # store SAMPLING_PLATFORM in survey_track so they can be used later
    survey_track['SAMPLING_PLATFORM'] = SAMPLING_PLATFORM
    return survey_track, survey_indices, sampling_details
 
    
def survey_interp(ds, survey_track, survey_indices, sampling_details):
    """Interpolates dataset 'ds' along the survey track given by the sruvey coordinates.


    Args:
        ds (xarray.core.dataset.Dataset): MITgcm LLC4320 data for all days
        survey_track (xarray.core.dataset.Dataset): lat,lon,dep,time of the survey used for the interpolation
        survey_indices (xarray.core.dataset.Dataset): i,j,k coordinates used for the interpolation
        sampling_details (dict):Includes number of days, waypoints, and depth range, horizontal and vertical platform speed. These can typical (default) or user-specified, in the                                      case where user specfies only some of the details the default values will be used for rest.
        

    Returns:
        subsampled_data: all field interpolated onto the track
        sh_true: 'true' steric height along the track
        
    Raises: 
        Sampling strategy is invalid: If a sampling strategy is not specified or different from the available strategies - sim_utcd, glider, sim_mooring, wave_glider, saildrone
    

    """
      
        
    ## Create a new dataset to contain the interpolated data, and interpolate
    # for 'mooring', skip this step entirely - return an empty array for 'subsampled_data'
    SAMPLING_PLATFORM = survey_track['SAMPLING_PLATFORM']
    if SAMPLING_PLATFORM == 'mooring':
        subsampled_data = []
        
        # zgridded and times are simply zs, ta (i.e., don't interpolate to a finer grid than the mooring sampling gives)
        zgridded = survey_track['dep']
        times = survey_track['time']
        
        # -- initialize the dataset:
        sgridded = xr.Dataset(
            coords = dict(depth=(["depth"],zgridded),
                      time=(["time"],times))
        )
        # variable names (if DERIVED_VARIABLES is not set, don't load the vector quantities)
        if sampling_details['DERIVED_VARIABLES']:
            vbls3d = ['Theta','Salt','vorticity','steric_height', 'U', 'V']
            vbls2d = ['steric_height_true', 'Eta', 'KPPhbl', 'PhiBot', 'oceTAUX', 'oceTAUY', 'oceFWflx', 'oceQnet', 'oceQsw', 'oceSflux']
        else:
            vbls3d = ['Theta','Salt']
            vbls2d = ['Eta', 'KPPhbl', 'PhiBot', 'oceFWflx', 'oceQnet', 'oceQsw', 'oceSflux']
        
        
        # loop through 3d variables & interpolate:
        for vbl in vbls3d:
            print(vbl)
            sgridded[vbl]=ds[vbl].interp(survey_indices).compute().transpose()

        # loop through 2d variables & interpolate:
        # create 2-d survey track by removing the depth dimension
        survey_indices_2d =  survey_indices.drop_vars('k')    
        for vbl in vbls2d:
            print(vbl)
            sgridded[vbl]=ds[vbl].interp(survey_indices_2d).compute()
            
        
        # clean up sgridded: get rid of the dims we don't need and rename coords        
        if sampling_details['DERIVED_VARIABLES']:
            sgridded = sgridded.reset_coords(names = {'i', 'j', 'k'}).squeeze().rename_vars({'xav' : 'lon','yav' : 'lat'}).drop_vars(names={'i', 'j', 'k'})
        
            # for sampled steric height, we want the value integrated from the deepest sampling depth:
            sgridded['steric_height'] = (("time"), sgridded['steric_height'].isel(depth=int(len(zgridded))-1))
            # rename to "sampled" for clarity
            sgridded.rename_vars({'steric_height':'steric_height_sampled'})
        else:
            sgridded = sgridded.reset_coords(names = {'i', 'j', 'k'}).squeeze().drop_vars(names={'i', 'j', 'k'})
    
    else:
        subsampled_data = xr.Dataset(
            dict(
                t = xr.DataArray(survey_track.time, dims='points'), # call this time, for now, so that the interpolation works
                lon = xr.DataArray(survey_track.lon, dims='points'),
                lat = xr.DataArray(survey_track.lat, dims='points'),
                dep = xr.DataArray(survey_track.dep, dims='points'),
                points = xr.DataArray(survey_track.points, dims='points')
            )
        )

        # variable names (if DERIVED_VARIABLES is not set, don't load the vector quantities)
        if sampling_details['DERIVED_VARIABLES']:
            vbls3d = ['Theta','Salt','vorticity','steric_height', 'U', 'V']
            vbls2d = ['steric_height_true', 'Eta', 'KPPhbl', 'PhiBot', 'oceTAUX', 'oceTAUY', 'oceFWflx', 'oceQnet', 'oceQsw', 'oceSflux']
        else:
            vbls3d = ['Theta','Salt']
            vbls2d = ['Eta', 'KPPhbl', 'PhiBot', 'oceFWflx', 'oceQnet', 'oceQsw', 'oceSflux']        
        
        print('Interpolating model fields to the sampling track...')
        # loop & interpolate through 3d variables:
        for vbl in vbls3d:
            subsampled_data[vbl]=ds[vbl].interp(survey_indices)       

        # loop & interpolate through 2d variables:
        # create 2-d survey track by removing the depth dimension
        survey_indices_2d =  survey_indices.drop_vars('k')
        for vbl in vbls2d:
            subsampled_data[vbl]=ds[vbl].interp(survey_indices_2d)   

        # fix time, which is currently a coordinate (time) & a variable (t)
        subsampled_data = subsampled_data.reset_coords('time', drop=True).rename_vars({'t':'time'})

        # make xav and yav variables instead of coords, and rename
        if sampling_details['DERIVED_VARIABLES']:
            subsampled_data = subsampled_data.reset_coords(names = {'xav','yav'}).rename_vars({'xav' : 'lon_average','yav' : 'lat_average'})
            
        
        if sampling_details['SAVE_PRELIMINARY']:
            # ----- save preliminary data
            # (not applicable to mooring data)
            # add metadata to attributes
            attrs = sampling_details
            attrs['start_date'] = sampling_details['start_date'].strftime('%Y-%m-%d')
            end_date = subsampled_data['time'].data[-1]
            attrs['end_date'] = np.datetime_as_string(end_date,unit='D')
            attrs.pop('DERIVED_VARIABLES')    
            attrs.pop('SAVE_PRELIMINARY')
            
            # filename:
            filename_out = sampling_details['filename_out_base'] + '_subsampled.nc'
            print(f'saving to {filename_out}')
            subsampled_data.attrs = attrs
            netcdf_fill_value = nc4.default_fillvals['f4']
            dv_encoding={'zlib':True,  # turns compression on\
                        'complevel':9,     # 1 = fastest, lowest compression; 9=slowest, highest compression \
                        'shuffle':True,    # shuffle filter can significantly improve compression ratios, and is on by default \
                        'dtype':'float32',\
                        '_FillValue':netcdf_fill_value}
            # save to a new file
            subsampled_data.to_netcdf(filename_out,format='netcdf4')            
            
            
        # -----------------------------------------------------------------------------------
        # ------Regrid the data to depth/time (3-d fields) or subsample to time (2-d fields)
        print('Gridding the interpolated data...')
        
        # get times associated with profiles:
        if SAMPLING_PLATFORM == 'sim_mooring':
            # - for mooring, use the subsampled time grid:
            times = np.unique(subsampled_data.time.values)
        else:
            # -- for glider/uctd, take the shallowest & deepest profiles (every second value, since top/bottom get sampled twice for each profile)
            time_deepest = subsampled_data.time.where(subsampled_data.dep == subsampled_data.dep.min(), drop=True).values[0:-1:2]
            time_shallowest = subsampled_data.time.where(subsampled_data.dep == subsampled_data.dep.max(), drop=True).values[0:-1:2]
            times = np.sort(np.concatenate((time_shallowest, time_deepest)))
            # this results in a time grid that may not be uniformly spaced, but is correct
            # - for a uniform grid, use the mean time spacing - may not be perfectly accurate, but is evenly spaced
            dt = np.mean(np.diff(time_shallowest))/2 # average spacing of profiles (half of one up/down, so divide by two)
            times_uniform = np.arange(survey_track.n_profiles.values*2) * dt

        # nt is the number of profiles (times):
        nt = len(times)  
        # xgr is the vertical grid; nz is the number of depths for each profile
        # depths are negative, so sort in reverse order using flip
        zgridded = np.flip(np.unique(subsampled_data.dep.data))
        nz = int(len(zgridded))

        # -- initialize the dataset:
        sgridded = xr.Dataset(
            coords = dict(depth=(["depth"],zgridded),
                      time=(["time"],times))
        )
        # -- 3-d fields: loop & reshape 3-d data from profiles to a 2-d (depth-time) grid:
        # first, extract each variable, then reshape to a grid        
        for vbl in vbls3d:
            print(f'  {vbl}')
            if sampling_details['SAVE_PRELIMINARY']:
                # not a dask array, so no "compute" command needed
                this_var = subsampled_data[vbl].data.copy() 
            else:
                this_var = subsampled_data[vbl].data.compute().copy() 
            # reshape to nz,nt
            this_var_reshape = np.reshape(this_var,(nz,nt), order='F') # fortran order is important!
            # for platforms with up & down profiles (uCTD and glider),
            # every second column is upside-down (upcast data)
            # starting with the first column, flip the data upside down so that upcasts go from top to bottom
            if SAMPLING_PLATFORM != 'sim_mooring':
                this_var_fix = this_var_reshape.copy()
                #this_var_fix[:,0::2] = this_var_fix[-1::-1,0::2] 
                this_var_fix[:,1::2] = this_var_fix[-1::-1,1::2]  # Starting with SECOND column
                sgridded[vbl] = (("depth","time"), this_var_fix)
            elif SAMPLING_PLATFORM == 'sim_mooring':
                sgridded[vbl] = (("depth","time"), this_var_reshape)                
                
        if sampling_details['DERIVED_VARIABLES']:
            # for sampled steric height, we want the value integrated from the deepest sampling depth:
            sgridded['steric_height'] = (("time"), sgridded['steric_height'].isel(depth=nz-1).data)
            # rename to "steric_height_sampled" for clarity
            sgridded.rename_vars({'steric_height':'steric_height_sampled'})
  

        #  -- 2-d fields: loop & reshape 2-d data to the same time grid 
        for vbl in vbls2d:
            
            if sampling_details['SAVE_PRELIMINARY']:
                # not a dask array, so no "compute" command needed
                this_var = subsampled_data[vbl].data.copy() 
            else:
                this_var = subsampled_data[vbl].data.compute().copy() 
            # subsample to nt
            this_var_sub = this_var[0:-1:nz]
            sgridded[vbl] = (("time"), this_var_sub)

    # ------------ RETURN INTERPOLATED & GRIDDED DATA ------------

    # -- add variable attributes from ds
    if SAMPLING_PLATFORM == 'mooring':
        sgridded.attrs = ds.attrs
    else:
        # - find which variables in ds are also in our interpolated dataset:
        vars_ds = list(ds.keys())
        vars_sdata = list(subsampled_data.keys())
        vars_both = list(set(vars_ds) & set(vars_sdata))
        for var in vars_both:
            # copy over the attribute from ds:
            subsampled_data[var].attrs = ds[var].attrs
            sgridded[var].attrs = ds[var].attrs
    
    
    
    return subsampled_data, sgridded


# great circle distance (from Jake Steinberg) 
def great_circle(lon1, lat1, lon2, lat2):
    """Computes great-circle distance


    Args:
        

    Returns:
        
    Raises: 
    

    """
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    return 6371 * (acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2)))




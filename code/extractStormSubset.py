#!/usr/bin/env python
# Import python modules
import psycopg2, subprocess, datetime, re, sys, pdb
from psycopg2.extensions import AsIs
from ast import literal_eval
import pandas as pd
import geopandas as gpd
import netCDF4 as nc
import numpy as np

# This function extracts a single timestamp of data for the region defined in the roi.json file
def getData(stormvartable, timestamp, roi):
    # Get tile coordinates base on region defined in the roi.json file
    std_result = subprocess.run(['geodex', roi, 
                                 str(13), '--output-format', '({z}, {x}, {y})'], stdout=subprocess.PIPE)
    tile_inds_string = std_result.stdout.decode('utf-8').split('\n')
    tile_inds_tuple = [literal_eval(ti) for ti in tile_inds_string if len(ti) > 0]

    # Create empty Geopandas Dataframe to append data 
    gdfall = gpd.GeoDataFrame()

    # Loop through tile coordinates and use to extract data from database.
    for tile in tile_inds_tuple:
        # Defind zoom, x, and y coordinates from tile
        zoom = str(tile[0])
        x = str(tile[1])
        y = str(tile[2])

        # Try connecting to database
        try:
            # Connect to database using the ajay account
            conn = psycopg2.connect("dbname='reg3sim' user='ajay' host='localhost' port='5432' password='adcirc'")
            cur = conn.cursor()

            # If storm table being queried is a fort63 table run query
            if stormvartable[18:] == 'fort63':
                cur.execute("""WITH
                                   bounds AS (
                                       SELECT ST_TileEnvelope(%(vzoom)s, %(vx)s, %(vy)s) AS geomclip,
                                       ST_Expand(ST_TileEnvelope(%(vzoom)s, %(vx)s, %(vy)s), 0) AS geombuf
                                   ),
                                   node_ids AS
                                       (SELECT G.node, G.geom, G.bathymetry
                                        FROM r3sim_fort_geom AS G, bounds
                                        WHERE ST_Intersects(G.geom, ST_Transform(bounds.geombuf, 4326)))
                                   SELECT ST_X(T.geom) AS longitude, 
                                          ST_Y(T.geom) AS latitude, 
                                          CAST(T.bathymetry AS DOUBLE PRECISION) AS bathymetry,
                                          S.node AS node, CAST(S.zeta AS DOUBLE PRECISION) AS zeta 
                                   FROM
                                       (SELECT node, zeta, timestamp
                                        FROM %(vstormvartable)s
                                        WHERE timestamp = %(vtimestamp)s
                                        AND node IN (SELECT node FROM node_ids)) S
                                       LEFT JOIN
                                       (SELECT G.geom AS geom, 
                                        G.node AS node, G.bathymetry AS bathymetry
                                        FROM node_ids G, bounds
                                        WHERE ST_Intersects(G.geom, ST_Transform(bounds.geombuf, 4326))) T
                                        ON S.node = T.node
                                        ORDER BY S.node;
                            """,
                            {'vzoom':AsIs(zoom),'vx':AsIs(x),'vy':AsIs(y),'vstormvartable':AsIs(stormvartable),'vtimestamp':timestamp})

                # Convert data queried to DataFrame
                df = pd.DataFrame(cur.fetchall(), columns=['longitude','latitude','bathymetry','node','zeta'])

            # If storm table being queried is a swan63 table run query
            elif stormvartable[18:] == 'swan63':
                cur.execute("""WITH
                                   bounds AS (
                                       SELECT ST_TileEnvelope(%(vzoom)s, %(vx)s, %(vy)s) AS geomclip,
                                       ST_Expand(ST_TileEnvelope(%(vzoom)s, %(vx)s, %(vy)s), 0) AS geombuf
                                   ),
                                   node_ids AS
                                       (SELECT G.node, G.geom, G.bathymetry
                                        FROM r3sim_fort_geom AS G, bounds
                                        WHERE ST_Intersects(G.geom, ST_Transform(bounds.geombuf, 4326)))
                                   SELECT ST_X(T.geom) AS longitude, 
                                          ST_Y(T.geom) AS latitude, 
                                          CAST(T.bathymetry AS DOUBLE PRECISION) AS bathymetry,
                                          S.node AS node, CAST(S.hs AS DOUBLE PRECISION) AS hs, 
                                          CAST(S.tps AS DOUBLE PRECISION) AS tps,CAST(S.dir AS DOUBLE PRECISION) AS dir
                                   FROM
                                       (SELECT node, hs, tps, dir, timestamp
                                        FROM %(vstormvartable)s
                                        WHERE timestamp = %(vtimestamp)s
                                        AND node IN (SELECT node FROM node_ids)) S
                                       LEFT JOIN
                                       (SELECT G.geom AS geom, 
                                        G.node AS node, G.bathymetry AS bathymetry
                                        FROM node_ids G, bounds
                                        WHERE ST_Intersects(G.geom, ST_Transform(bounds.geombuf, 4326))) T
                                        ON S.node = T.node
                                        ORDER BY S.node;
                            """,
                            {'vzoom':AsIs(zoom),'vx':AsIs(x),'vy':AsIs(y),'vstormvartable':AsIs(stormvartable),'vtimestamp':timestamp})

                # Convert data queried to DataFrame
                df = pd.DataFrame(cur.fetchall(), columns=['longitude','latitude','bathymetry','node','hs','tps','dir'])

            # If storm table being queried is a fort64 table run query
            elif stormvartable[18:] == 'fort64':
                cur.execute("""WITH
                                   bounds AS (
                                       SELECT ST_TileEnvelope(%(vzoom)s, %(vx)s, %(vy)s) AS geomclip,
                                       ST_Expand(ST_TileEnvelope(%(vzoom)s, %(vx)s, %(vy)s), 0) AS geombuf
                                   ),
                                   node_ids AS
                                       (SELECT G.node, G.geom, G.bathymetry
                                        FROM r3sim_fort_geom AS G, bounds
                                        WHERE ST_Intersects(G.geom, ST_Transform(bounds.geombuf, 4326)))
                                   SELECT ST_X(T.geom) AS longitude,
                                          ST_Y(T.geom) AS latitude,
                                          CAST(T.bathymetry AS DOUBLE PRECISION) AS bathymetry,
                                          S.node AS node, CAST(S.u_vel AS DOUBLE PRECISION) AS u_vel,
                                          CAST(S.v_vel AS DOUBLE PRECISION) AS v_vel
                                   FROM
                                       (SELECT node, u_vel, v_vel, timestamp
                                        FROM %(vstormvartable)s
                                        WHERE timestamp = %(vtimestamp)s
                                        AND node IN (SELECT node FROM node_ids)) S
                                       LEFT JOIN
                                       (SELECT G.geom AS geom,
                                        G.node AS node, G.bathymetry AS bathymetry
                                        FROM node_ids G, bounds
                                        WHERE ST_Intersects(G.geom, ST_Transform(bounds.geombuf, 4326))) T
                                        ON S.node = T.node
                                        ORDER BY S.node;
                            """,
                            {'vzoom':AsIs(zoom),'vx':AsIs(x),'vy':AsIs(y),'vstormvartable':AsIs(stormvartable),'vtimestamp':timestamp})

                # Convert data queried to DataFrame
                df = pd.DataFrame(cur.fetchall(), columns=['longitude','latitude','bathymetry','node','u_vel', 'v_vel'])

            # Convert DataFrame to a Geopandas DataFrame
            gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude))

            # Append DataFrame
            gdfall = pd.concat([gdfall,gdf])

        except (Exception, psycopg2.DatabaseError) as error:
            print(error)
        finally:
            if cur is not None:
                cur.close()
            if conn is not None:
                conn.close()

    # Sort data by node and reset index
    gdfall = gdfall.sort_values(by=['node'], ascending=True)
    gdfall = gdfall.reset_index(drop=True)

    # Check for duplicated node values
    duplicateRowGDF = gdfall[gdfall.duplicated(['node'])]
    
    if len(duplicateRowGDF) > 0:
        print(len(duplicateRowGDF))

    # Return DataFrame
    return(gdfall)

# Write queried data to netCDF file
def writeNetCDF(stormvartable, timesteps, roi):
    # Define output path and filename
    path = '/home/ajay/work/data/'
    filename = path+stormvartable+'_subset.nc'

    # Create netCDF data set
    ds = nc.Dataset(filename, 'w', format='NETCDF4')

    # Get data for first timestamp 
    gdf = getData(stormvartable, timesteps[0], roi)

    # Define and load timestamp data for netCDF file
    ts = list(map(int, re.split('T|-|:', timesteps[0])))
    timestamp = (datetime.datetime(ts[0],ts[1],ts[2],ts[3],ts[4],ts[5]) - 
                 datetime.datetime(2000,9,1,0,0,0)).total_seconds()
   
    # Define time and node dimensions for netCDF file
    time = ds.createDimension('time', None)
    nodes = ds.createDimension('node', gdf.shape[0])

    # Defind time field and load timestamp data in netCDF file
    time = ds.createVariable('time', np.float64, ('time',))
    time.long_name = 'model time'
    time.standard_name = 'time'
    time.units = 'seconds since 2000-09-01 00:00:00'
    time.base_date = '2000-09-01 00:00:00'
    time[:] = timestamp

    # Define longitude field and load longitude data in netCDF file
    lons = ds.createVariable('x', np.float32, ('node',))
    lons.long_name = 'longitude'
    lons.standard_name = 'longitude'
    lons.units = 'degrees_east'
    lons[:] = gdf.geometry.x[:]

    # Define latitude field and load latitude data in netCDF file
    lats = ds.createVariable('y', np.float32, ('node',))
    lats.long_name = 'latitude'
    lats.standard_name = 'latitude'
    lats.units = 'degrees_north'
    lats[:] = gdf.geometry.y[:]

    # Define depth field and load depth data in netCDF file 
    bathymetry = ds.createVariable('depth', np.float64, ('node',))
    bathymetry.long_name = 'Relative Peak Period'
    bathymetry.standard_name = 'depth_below_geoid'
    bathymetry.units = 'm'
    bathymetry[:] = gdf.bathymetry[:]

    # If data type is fort63 load data
    if stormvartable[18:] == 'fort63':
        # Define zeta field and load zeta data in netCDF file
        zeta = ds.createVariable('zeta', np.float64, ('time', 'node',))
        zeta.long_name = 'water surface elevation above geoid'
        zeta.standard_name = 'sea_surface_height_above_geoid'
        zeta.units = 'm'
        zeta[0, :] = gdf.zeta[:]
   
    # If data type is swan63 load data
    elif stormvartable[18:] == 'swan63':
        # Define hs field and load hs data in netCDF file
        hs = ds.createVariable('hs', np.float64, ('time', 'node',))
        hs.long_name = 'significant wave height'
        hs.standard_name = 'significant_wave_height'
        hs.units = 'm'
        hs[0, :] = gdf.hs[:]

        # Define tps field and load tps data in netCDF file
        tps = ds.createVariable('tps', np.float64, ('time', 'node',))
        tps.long_name = 'Relative Peak Period'
        tps.standard_name = 'sea_surface_wave_period_at_variance_spectral_density_maximum'
        tps.units = 's'
        tps[0, :] = gdf.tps[:]

        # Define dir field and load dir data in netCDF file
        wavedir = ds.createVariable('dir', np.float64, ('time', 'node',))
        wavedir.long_name = 'wave direction'
        wavedir.standard_name = 'wave_direction'
        wavedir.units = 'degrees_CW_from_East'
        wavedir[0, :] = gdf.dir[:]

    # If data type is fort64 load data
    elif stormvartable[18:] == 'fort64':
        # Define u_vel field and load u_vel data in netCDF file
        u_vel = ds.createVariable('u_vel', np.float64, ('time', 'node',))
        u_vel.long_name = 'water column vertically averaged east/west velocity'
        u_vel.standard_name = 'eastward_water_velocity'
        u_vel.units = 'm s-1'
        u_vel[0, :] = gdf.u_vel[:]

        # Define v_vel field and load v_vel data in netCDF file
        v_vel = ds.createVariable('v_vel', np.float64, ('time', 'node',))
        v_vel.long_name = 'water column vertically averaged north/south velocity'
        v_vel.standard_name = 'northward_water_velocity'
        v_vel.units = 'm s-1'
        v_vel[0, :] = gdf.v_vel[:]

    # Defind i for indexing
    i = 1

    # Loop through rest of timestamps and load data
    for timestep in timesteps[1:]:
        # Get data for other timestamps
        gdf = getData(stormvartable, timestep, roi)

        # Define and load timestamp data
        ts = list(map(int, re.split('T|-|:', timestep)))
        timestamp = (datetime.datetime(ts[0],ts[1],ts[2],ts[3],ts[4],ts[5]) - 
                     datetime.datetime(2000,9,1,0,0,0)).total_seconds()
        time[:] = timestamp
   
        # Load latitude and longitude data
        lons[:] = gdf.geometry.x[:]
        lats[:] = gdf.geometry.y[:]

        # If data type is fort63 load data
        if stormvartable[18:] == 'fort63':
            # Load zeta data
            zeta[i, :] = gdf.zeta[:]

        # If data type is swan63 load data
        elif stormvartable[18:] == 'swan63':
            # Load hs data
            hs[i, :] = gdf.hs[:]

            #Load tsp data
            tps[i, :] = gdf.tps[:]

            # Load dir data
            wavedir[i, :] = gdf.dir[:]

        # If data type is fort64 load data
        elif stormvartable[18:] == 'fort64':
            # Load u_vel data
            u_vel[i, :] = gdf.u_vel[:]

            # Load v_vel data
            v_vel[i, :] = gdf.v_vel[:]

        # Itterate i
        i = i + 1

# Defind roi file path
roi = '/home/ajay/work/code/roi.geojson'

# Read timestamp csv file and output to list
dft = pd.read_csv('/home/ajay/work/code/timestamps.csv', header=None)
timesteps = dft.values.tolist()[0]

# Defind stormvartable from sys.argv input
stormvartable = sys.argv[1]

# Run writeNetCDF
writeNetCDF(stormvartable, timesteps, roi)

import numpy as np
import json
import os
from pysheds.grid import Grid
from pysheds import io
from geojson import FeatureCollection, LineString
import shapely.geometry
import rasterio
import rioxarray
from scipy.spatial import Delaunay
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pyvista as pv


def geojson2shapely(line):
    assert (isinstance(line, LineString))
    return shapely.geometry.LineString(line["coordinates"])


def generate_watersheds(elevation):

    grid = Grid.from_raster(elevation)
    dem = grid.read_raster(elevation)
    xmin, xmax, ymin, ymax = dem.extent

    padding = 2000
    xmin = xmin + padding
    xmax = xmax - padding
    ymax = ymax - padding
    ymin = ymin + padding

    line_top = LineString(coordinates = [(xmin, ymax), (xmax, ymax)])
    line_top = geojson2shapely(line_top)

    line_bottom = LineString(coordinates = [(xmin, ymin), (xmax, ymin)])
    line_bottom = geojson2shapely(line_bottom)

    line_left = LineString(coordinates = [(xmin, ymin), (xmin, ymax)])
    line_left = geojson2shapely(line_left)

    line_right = LineString(coordinates = [(xmax, ymin), (xmax, ymax)])
    line_right = geojson2shapely(line_right)

    shapely_extent = shapely.geometry.Polygon([

            (xmin, ymin), (xmin, ymax)
        ,

            (xmin, ymax), (xmax, ymax)
        ,

            (xmax, ymax), (xmax, ymin)
        ,

            (xmax, ymin), (xmin, ymin)

    ])
    with open('boundary.geojson', 'w') as fptr:
        fptr.write(json.dumps(shapely_extent.__geo_interface__))


    # Condition DEM
    # ----------------------
    # Fill pits in DEM
    pit_filled_dem = grid.fill_pits(dem)

    # Fill depressions in DEM
    flooded_dem = grid.fill_depressions(pit_filled_dem)

    # Resolve flats in DEM
    inflated_dem = grid.resolve_flats(flooded_dem)

    # Determine D8 flow directions from DEM
    # ----------------------
    # Specify directional mapping
    dirmap = (64, 128, 1, 2, 4, 8, 16, 32)

    # Compute flow directions
    # -------------------------------------
    fdir = grid.flowdir(inflated_dem, dirmap=dirmap)

    # Calculate flow accumulation
    # --------------------------
    acc = grid.accumulation(fdir, dirmap=dirmap)

    # Extract river network
    # ---------------------
    branches = grid.extract_river_network(fdir, acc > 50, dirmap=dirmap)

    with open('branches.geojson', 'w') as fptr:
        fptr.write(json.dumps(branches.__geo_interface__))

    # extract outlets
    outlets = []
    shapely_branches = []
    for branch in branches['features']:
        shapely_branch = geojson2shapely(branch['geometry'])
        for line in [line_bottom, line_right, line_left, line_top]:
            intersection = line.intersection(shapely_branch)
            if intersection:
                outlets.append(intersection.centroid)

    i = 0
    for outlet in outlets:
        try:
            catch = grid.catchment(x=outlet.x, y=outlet.y, fdir=fdir, xytype='coordinate')
            # save
            inflated_dem.mask = catch
            ofile = f'w_{i}.tif'
            if os.path.exists(ofile):
                os.remove(ofile)
            io.to_raster(inflated_dem, ofile, dtype=np.float)

        except ValueError as ex:
            # pour point out of bounds
            print(str(ex))
            pass
        i += 1


def generate_tin(wshed):
    rds = rioxarray.open_rasterio(wshed)
    rds = rds.squeeze().drop("spatial_ref").drop("band")
    rds.name = "data"
    df = rds.to_dataframe().reset_index()
    sub_df = df[df.data >= 0.0]
    sub_df.to_csv("out.csv", index=False)

    points = list(zip(sub_df['x'], sub_df['y'], sub_df['data']))
    mesh = pv.PolyData(points)
    surf = mesh.delaunay_2d()
    surf = surf.decimate(0.75)
    surf.plot(show_edges=True)

    # mesh.plot(point_size=10, style='points')

    # pl = pv.Plotter()
    # pl.add_mesh(surf, show_edges=True, color='white')
    # # pl.add_points(mesh.points, color='red',
    # #               point_size=20)
    # # pl.camera_position = cpos
    # pl.show()

    # points = list(zip(sub_df['x'], sub_df['y']))
    # tri = Delaunay(points)
    # import matplotlib.pyplot as plt
    # plt.triplot(sub_df['x'], sub_df['y'], tri.simplices)
    # plt.plot(sub_df['x'], sub_df['y'], 'o')
    # plt.show()


if __name__ == '__main__':
    # elevation = '../../sw_sub_epsg2163_eq_area.tif'
    # generate_watersheds(elevation)
    generate_tin('w_18.tif')
    # generate_tin('w_37.tif')

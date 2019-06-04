# pylint: disable=missing-docstring,invalid-name,undefined-loop-variable

import datetime
import io
import json
from math import radians, cos, sin, asin, sqrt

import click
import fiona
import geopandas as gpd
import matplotlib.pyplot as plt  # NOQA
import numpy
import pandas as pd
import shapely


def to_miles(m):
    return m * 0.000621371


def json_default(v):
    if isinstance(v, numpy.int64):
        return int(v)
    if isinstance(v, pd.Timestamp):
        return '{}'.format(v)
    if isinstance(v, pd._libs.tslibs.nattype.NaTType):
        return None

    return json.JSONEncoder().default(v)


def read_file(path):
    src = fiona.open(path, layer='tracks')
    meta = src.meta
    meta['driver'] = 'GeoJSON'

    with io.BytesIO() as buf:
        with fiona.open(buf, 'w', **meta) as dst:
            for feature in src:
                if len(feature['geometry']['coordinates'][0]) > 1:
                    dst.write(feature)

        buf.seek(0)
        tracks = gpd.read_file(buf, driver='GeoJSON')

    points = gpd.read_file(path, layer='track_points')
    points['time'] = pd.to_datetime(points['time'])

    return tracks, points


def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))
    r = 6371000  # Radius of earth in m.
    return c * r


def stats_for(points, trackid):
    track = points[points.track_fid == trackid]

    i1 = track.itertuples(index=True)
    i2 = track.itertuples(index=True)
    next(i2)

    d_total = 0
    t_total = datetime.timedelta()
    t_moving = datetime.timedelta()

    for point in zip(i1, i2):
        delta_d = haversine(point[0].geometry.y, point[0].geometry.x,
                            point[1].geometry.y, point[1].geometry.x)
        delta_t = point[1].time - point[0].time
        delta_e = point[1].ele - point[0].ele
        speed = delta_d/delta_t.total_seconds()

        points.loc[point[1].Index, 'delta_t'] = delta_t
        points.loc[point[1].Index, 'delta_d'] = delta_d
        points.loc[point[1].Index, 'delta_e'] = delta_e
        points.loc[point[1].Index, 'speed'] = speed
        points.loc[point[1].Index, 'segment'] = (
            shapely.geometry.LineString((
                (point[0].geometry.x, point[0].geometry.y),
                (point[1].geometry.x, point[1].geometry.y)
            ))
        )

        if speed >= 1:
            t_moving += delta_t
        t_total += delta_t
        d_total += delta_d


@click.command()
@click.option('-m', '--min-length', type=float)
@click.argument('gpxfiles', nargs=-1)
def main(min_length=None, gpxfiles=None):
    for gpxfile in gpxfiles:
        tracklines, points = read_file(gpxfile)

        points['segment'] = shapely.geometry.LineString()
        points['speed'] = 0
        points['delta_d'] = 0
        points['delta_t'] = 0
        points['delta_e'] = 0

        for trackid in points.track_fid.unique():
            stats = stats_for(points, trackid)
            if not stats:
                continue
            if min_length and stats['d_total'] < min_length:
                continue

        print(gpxfile)
        import code
        code.interact(local=locals())


if __name__ == '__main__':
    main()

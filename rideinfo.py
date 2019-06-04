# pylint: disable=missing-docstring,invalid-name

import datetime
import io
import json
from math import radians, cos, sin, asin, sqrt

import click
import fiona
import geopandas as gpd
import numpy
import pandas as pd


def to_miles(m):
    return m * 0.000621371


def json_default(v):
    if isinstance(v, numpy.int64):
        return int(v)
    elif isinstance(v, pd.Timestamp):
        return '{}'.format(v)

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


def stats_for(trackid, track):
    d_total = 0
    i1 = track.itertuples()
    i2 = track.itertuples()
    next(i2)

    d_total = 0
    t_total = datetime.timedelta()
    t_moving = datetime.timedelta()

    stats = {'trackid': trackid}
    for i, point in enumerate(zip(i1, i2)):
        if i == 0:
            stats['t_start'] = point[0].time

        d_delta = haversine(point[0].geometry.y, point[0].geometry.x,
                            point[1].geometry.y, point[1].geometry.x)
        t_delta = point[1].time - point[0].time
        speed = d_delta/t_delta.total_seconds()
        if speed >= 1:
            t_moving += t_delta
        t_total += t_delta
        d_total += d_delta

    stats['t_end'] = point[1].time
    stats['t_total'] = t_total.total_seconds()
    stats['t_moving'] = t_moving.total_seconds()
    stats['d_total'] = d_total

    stats['avg_speed'] = d_total/stats['t_total']
    stats['mov_speed'] = d_total/stats['t_moving'] if stats['t_moving'] else 0

    stats['avg_speed_mph'] = stats['avg_speed'] * 2.23694
    stats['mov_speed_mph'] = stats['mov_speed'] * 2.23694

    return stats

@click.command()
@click.argument('gpxfiles', nargs=-1)
def main(gpxfiles):
    for gpxfile in gpxfiles:
        tracklines, points = read_file(gpxfile)

        trackpoints = {}
        for trackid in points.track_fid.unique():
            trackpoints[trackid] = points.query(
                'track_fid == {}'.format(trackid))

        for trackid, track in trackpoints.items():
            if len(track) <= 1:
                continue

            stats = stats_for(trackid, track)
            print(json.dumps(stats, indent=2, default=json_default))

if __name__ == '__main__':
    main()

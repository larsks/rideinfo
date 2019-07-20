import datetime
import io
import json
from math import radians, cos, sin, asin, sqrt

import click
import fiona
import geopandas as gpd
import matplotlib.pyplot as plt  # NOQA
import numpy as np
import pandas as pd
import shapely


def to_miles(m):
    return m * 0.000621371


def to_mph(m):
    return m * 2.237


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


def haversine_np(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(np.radians,
                                 [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2
    c = 2 * np.arcsin(np.sqrt(a))
    r = 6371000  # Radius of earth in m.
    return c * r


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
    track = points.track_fid == trackid

    t1 = points.loc[track, :].shift().iloc[1:]
    t2 = points.loc[track, :].iloc[1:]
    res = haversine_np(t1.geometry.x, t1.geometry.y,
                       t2.geometry.x, t2.geometry.y)

    points.loc[track, 'delta_d'] = res
    points.loc[track, 'delta_t'] = t2.time - t1.time
    points.loc[track, 'delta_e'] = t2.ele - t1.ele

    points.loc[track, 'speed'] = (
        points.loc[track, 'delta_d'] /
        points.loc[track, 'delta_t'].dt.total_seconds()
    )
#        points.loc[point[1].Index, 'segment'] = (
#            shapely.geometry.LineString((
#                (point[0].geometry.x, point[0].geometry.y),
#                (point[1].geometry.x, point[1].geometry.y)
#            ))
#        )


@click.command()
@click.option('-l', '--min-length', type=float, default=0)
@click.option('-m', '--moving-speed', type=float, default=0)
@click.argument('gpxfiles', nargs=-1)
def main(min_length, moving_speed, gpxfiles):
    agg_tracks = []
    for gpxfile in gpxfiles:
        tracklines, points = read_file(gpxfile)

        tracks = points.groupby(['track_fid'])
        for track_fid in tracks.groups.keys():
            stats_for(points, track_fid)

        for track_fid in tracks.groups.keys():
            if tracks.sum()['delta_d'][track_fid] < min_length:
                continue

            agg_tracks.append((
                track_fid,
                tracks.first().time[track_fid],
                tracks.sum()['delta_d'][track_fid],
                tracks.mean()['speed'][track_fid],
                tracks.get_group(track_fid)[
                    tracks.get_group(track_fid).speed > moving_speed
                ].mean()['speed'],
            ))

    agg_tracks_df = pd.DataFrame(
        agg_tracks,
        columns=['track', 'time', 'distance', 'speed', 'moving_speed'])

    breakpoint()

    agg_tracks_df.to_csv('stats.csv')
    
    weeks = agg_tracks_df.groupby(pd.Grouper(key='time', freq='W-Sun'))
    for week in weeks.groups.keys():
        print(week.strftime('%Y-%m-%d'),
              to_miles(weeks.sum()['distance'][week]),
              to_mph(weeks.mean()['speed'][week]),
              to_mph(weeks.mean()['moving_speed'][week]),
              )

if __name__ == '__main__':
    main()

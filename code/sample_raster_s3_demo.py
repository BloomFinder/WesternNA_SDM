#!/usr/bin/env python

import re
import boto3
import rasterio

##Grabs summaries of rasters in S3.
prefix = f'cog/'
summaries = list(objsum for objsum in boto3.resource('s3').Bucket('bloomfindersdm').objects.filter(Prefix=prefix))

keys = [None] * len(summaries)

for i in range(len(summaries)):
        keys[i] = summaries[i].key
del keys[0]
#print(keys)

##Samples all raster values at coordinates.

coords = [(-13460634, 6095878)]

for k in keys:
    with rasterio.open(f's3://bloomfindersdm/{k}') as src:
        
        vals = src.sample(coords)
        print(k)
        for val in vals:
            print(list(val))
        src.close()

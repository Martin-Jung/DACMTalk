#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Import necessary libraries
try:
    from osgeo import ogr
except ImportError:
    import ogr
try:
    from osgeo import osr
except ImportError:
    import osr
try:
    from osgeo import gdal
except ImportError:
    import gdal
import numpy as np
from datetime import date
import tempfile, os
try:
    import fiona
except ImportError:
    print "This function needs fiona to be installed"
try:
    import shapely
    from shapely.geometry import shape
except ImportError:
    print "This function needs shapely to be installed"



# --------------------------------------------------------------- #

def reprojectSHP(shpPath,outputFile=True,inputEPSG = 4326, outputproj4 = "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext"):
    """
    Reproject the input shapefile
    """
    
    driver = ogr.GetDriverByName('ESRI Shapefile')
    # create coordinate transformation
    inSpatialRef = osr.SpatialReference()
    inSpatialRef.ImportFromEPSG(inputEPSG)
    # The output
    outSpatialRef = osr.SpatialReference()
    outSpatialRef.ImportFromProj4(outputproj4)

    # Do the transformation
    coordTransform = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)
    
    # get the input layer
    inDataSet = driver.Open(shpPath)
    inLayer = inDataSet.GetLayer()
    
    # create the output layer
    if outputFile:
        p = tempfile.NamedTemporaryFile(suffix=".shp",prefix="reproj_",delete=False).name
    else:
        p = outputFile
    pn = os.path.basename(p)
    if os.path.exists(p):
        driver.DeleteDataSource(p)
    outDataSet = driver.CreateDataSource(p)
    outLayer = outDataSet.CreateLayer(pn, geom_type=ogr.wkbPoint)
    
    # add fields
    inLayerDefn = inLayer.GetLayerDefn()
    for i in range(0, inLayerDefn.GetFieldCount()):
        fieldDefn = inLayerDefn.GetFieldDefn(i)
        outLayer.CreateField(fieldDefn)
    
    # get the output layer's feature definition
    outLayerDefn = outLayer.GetLayerDefn()
    
    # loop through the input features
    inFeature = inLayer.GetNextFeature()
    while inFeature:
        # get the input geometry
        geom = inFeature.GetGeometryRef()
        # reproject the geometry
        geom.Transform(coordTransform)
        # create a new feature
        outFeature = ogr.Feature(outLayerDefn)
        # set the geometry and attribute
        outFeature.SetGeometry(geom)
        for i in range(0, outLayerDefn.GetFieldCount()):
            outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))
        # add the feature to the shapefile
        outLayer.CreateFeature(outFeature)
        # destroy the features and get the next input feature
        outFeature.Destroy()
        inFeature.Destroy()
        inFeature = inLayer.GetNextFeature()
    
    # Finally write a prj file
    outSpatialRef.MorphToESRI()
    file = open(str( os.path.splitext(p)[0] ) + ".prj", 'w')
    file.write(outSpatialRef.ExportToWkt())
    file.close()
    
    # close the shapefiles
    inDataSet.Destroy()
    outDataSet.Destroy()
    
    return p


def pointLatLong2Sinu(point, inputEPSG = 4326, outputproj4 = "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext"):
    # create coordinate transformation
    inSpatialRef = osr.SpatialReference()
    inSpatialRef.ImportFromEPSG(inputEPSG)
    # The output
    outSpatialRef = osr.SpatialReference()
    outSpatialRef.ImportFromProj4(outputproj4)

    coordTransform = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)

    # transform point
    try:
        point.Transform(coordTransform)
    except Exception:
        print("Error. Latitude or Longitude exceeded limits")
        return None
    
    # Return transformed point in resulting projection
    return point



def mapToPixel(self,geoMatrix, x, y):
    """
    Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate
    the pixel location of a geospatial coordinate 
    """
    ulX = geoMatrix[0]
    ulY = geoMatrix[3]
    xDist = geoMatrix[1]
    yDist = geoMatrix[5]
    rtnX = geoMatrix[2]
    rtnY = geoMatrix[4]
    pixel = int((x - ulX) / xDist)
    line = int((ulY - y) / xDist)
    return (pixel, line) 

# Save results to CSV
def saveToCSV( results, titles, filePath = "/research/ecocon_d/mj291/" ):
    """
    Saves a result table row by row into a csv
    """
    f = open(filePath, "wb" )
    writer = csv.writer(f,delimiter=';',quotechar="",quoting=csv.QUOTE_NONE)
    writer.writerow(titles)
    for item in results:
        writer.writerow(item)
    f.close()

def decomposeMODISFilename(file):
    """
    Decomposes a given original MODIS filename into its parameters
    Returns:
        A dictionary object containing Product, Tilename, h,v and the Date
    """
    theYear = file[9:13]
    theDay = file[13:16]
    prefix = file[0:9]
    tile = file[17:27]
    hi = tile.split('.')[0].find('h')
    vi = tile.split('.')[0].find('v')
    h = tile.split('.')[0][hi+1:hi+3]
    v = tile.split('.')[0][vi+1:vi+3]

    # convert day to month and year
    theYear = int(theYear)
    theDay = int(theDay)
    myDate = date.fromordinal(date(theYear, 1, 1).toordinal() + theDay - 1)
    
    myYear = myDate.year
    myMonth = myDate.month
    myDay = myDate.day
        
    outFile = prefix + "."  + str(myYear) + "." + str(myMonth) + "." + str(myDay) + "." + tile
    result = dict()
    result["Product"] = prefix
    result["FileName"] = os.path.splitext(file)[0]
    result["TileName"] = tile
    result["h"] = h
    result["v"] = v
    result["Year"] = myYear
    result["Month"] = myMonth
    result["Day"] = myDay
    
    return result


def length(x):
    """
    Short wrapper for length
    """
    return(len(x))
    

def ModisTileQuery(point,path2shp = 'modis_sinusoidal_grid_world.shp'):
    """
    Do a query of a given point shapefile on the MODISTiles
    Returns:
        A list with the respective h / v tiles of interest
    """

    tiles = [pol for pol in fiona.open(path2shp)]
    points = [pt for pt in fiona.open(point)]
    
    res = []
    for i, pt in enumerate(points):
        p = shape(pt['geometry'])
        #iterate through tiles
        for j, poly in enumerate(tiles):
            if p.within(shape(poly['geometry'])):
                # sum of attributes values
                h = int( tiles[j]['properties']['h'] )
                v = int( tiles[j]['properties']['v'] )
                res.append((h,v))

    return res


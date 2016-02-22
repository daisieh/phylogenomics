#!/usr/bin/env python

import re
import xmltodict
import sys
import json
import os
from svg.path import Path, Line, Arc, CubicBezier, QuadraticBezier, parse_path

total_width = 0
total_height = 0
def main():
    filename = sys.argv[1]
    outputfile = 'test.svg'
    
    outdict = {}
    outdict['svg'] = []
    file_name, extension = os.path.splitext(sys.argv[1])
    if extension != '.svg':
        print "can't open this file"
        return
    f = open(filename, 'r')
    xmldict = ''
    xml = f.read()
    xmldict = xmltodict.parse(xml)['svg']
    global total_width
    global total_height
    total_width = float(xmldict['@width'].replace('pt',''))
    total_height = float(xmldict['@height'].replace('pt',''))
    outdict['svg'] = {}
    outdict['svg']['width'] = xmldict['@width']
    outdict['svg']['height'] = xmldict['@height']
    outdict['svg']['path'] = []
    outdict['svg']['circle'] = []
    xmlpaths = []
    style = {}
    if 'path' in xmldict:
        xmlpaths = xmldict['path']
    elif 'g' in xmldict:
        xmlpaths = xmldict['g']['path']
        del xmldict['g']['path']
        style = "fill:#FF0000; stroke:none;"
        print style
    for path in xmlpaths:
#         boxpath = {}
        oldpath = {}
        if '@style' in path:
            matcher = re.match("^(.*fill:#)(.+?)(;.*$)", path['@style'])
            if matcher.group(0) is not None:
                fill = matcher.group(2)
                avg_col = average_color(fill)
                if avg_col < 0x50:
                    oldpath['@style'] = matcher.group(1) + 'FF0000' + matcher.group(3)
        else:
            oldpath['@style'] = style
        polygon = path_to_polygon(path['@d'])
        next = bounding_box(polygon)
        if threshold_area(bounding_box(polygon), 0.4):
            oldpath['@d'] = polygon_to_path(polygon)
        if ('@style' in oldpath) and ('@d' in oldpath):
            outdict['svg']['path'].append(oldpath) 
            outdict['svg']['circle'].extend(polygon_to_circles(polygon))    
#     print json.dumps(outdict)    
    outf = open(outputfile,'w')
    outf.write(xmltodict.unparse(outdict, pretty=True))
    outf.close()
    
def path_to_polygon_old(path):
    polygon = []
    nodes = re.findall('([MLHVCSQTAZ][-*\d\.\s]+)', path, re.I)
    start_x = 0
    start_y = 0
    print "hi"
    for n in nodes:
        print n
        movematch = re.match('(M)(.+?) (.+?)$', n, re.I)
        linematch = re.match('(L)(.+?) (.+?)$', n, re.I)
        curvematch = re.match('(C)(.+?) (.+?) (.+?) (.+?) (.+?) (.+?)$', n, re.I)
        if curvematch is not None:
            if curvematch.group(1) == 'C':
                polygon.append('%s %s' % (curvematch.group(6), curvematch.group(7)))
            else:
                print 'relative'
                polygon.append('%s %s' % (float(curvematch.group(6)) + start_x, float(curvematch.group(7)) + start_y))
        elif movematch is not None:
            if movematch.group(1) == 'M':
                start_x = float(movematch.group(2))
                start_y = float(movematch.group(3))
                polygon.append('%s %s' % (movematch.group(2), movematch.group(3)))
            else:
                polygon.append('%s %s' % (float(movematch.group(2)) + start_x, float(movematch.group(3)) + start_y))
        elif linematch is not None:
            print linematch.group(0)
            if linematch.group(1) == 'L':
                polygon.append('%s %s' % (linematch.group(2), linematch.group(3)))
            else:
                polygon.append('%s %s' % (float(linematch.group(2)) + start_x, float(linematch.group(3)) + start_y))
    return polygon

def path_to_polygon(path):
    polygon = []
    new_path = Path()    
    for segment in parse_path(path):
        new_path.append(Line(segment.start, segment.end))
    new_path.closed = True
    nodes = re.findall('[ML]\s*(\d+\.*\d*,\d+\.*\d*)\s*', new_path.d())
    for n in nodes:
        coords = n.split(',')
        polygon.append('%s %s' % (coords[0], coords[1]))
    return polygon

def polygon_to_path(polygon):
    return 'M' + 'L'.join(polygon) + 'Z'
    
def rotate_polygon(polygon):
    curr_min = re.split(' ', polygon[0])
    for i in range(len(polygon)):
        polygon.append(polygon.pop(0))
        curr_point = re.split(' ', polygon[0])
      
def polygon_to_circles(polygon):
    circlelist = []
    for i in range(len(polygon)):
        circledict = {}
        coords = re.split(' ', polygon[i])
        circledict['@r'] = '3'
        circledict['@stroke'] = 'none'
        circledict['@fill'] = '#c6c6%(number)02x' % {"number":(i*16)}
        circledict['@cx'] = coords[0]
        circledict['@cy'] = coords[1]
        circlelist.append(circledict)
#     print circledict
    return circlelist

def average_color(col):
    matcher = re.match("(..)(..)(..)", col)
    r = matcher.group(1)
    g = matcher.group(2)
    b = matcher.group(3)
    average = (int(r,16) + int(g,16) + int(b,16))/3
    return average

def threshold_area(rect, threshold):
    mins = re.split(' ', rect[0])
    maxs = re.split(' ', rect[2])
    min_x = float(mins[0])
    min_y = float(mins[1])
    max_x = float(maxs[0])
    max_y = float(maxs[1])
    if abs(float(max_y - min_y) / float(total_height)) > float(threshold):
        return True
    if abs(float(max_x - min_x) / float(total_width)) > float(threshold):
        return True
    return False
    
def bounding_box(polygon):
#     print polygon
    x_points = []
    y_points = []
    for point in polygon:
        coord = re.split(' ', point)
        x_points.append(coord[0])
        y_points.append(coord[1])
    max_x = max(x_points)
    max_y = max(y_points)
    min_x = min(x_points)
    min_y = min(y_points)
#     print total_width
#     print abs((float(max_x) - float(min_x)) / float(total_width))
#     if (abs((float(max_x) - float(min_x)) / float(total_width)) > 0.04) or (abs((float(max_y) - float(min_y)) / float(total_height)) > 0.04):
    return ['%s %s' % (min_x, min_y), '%s %s' % (min_x, max_y), '%s %s' % (max_x, max_y), '%s %s' % (max_x, min_y)] 
    
if __name__ == '__main__':
    main()

#!/usr/bin/env python

import re
import xmltodict
import sys
import json
import os

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
    xmldict = xmltodict.parse(xml)
    global total_width
    global total_height
    total_width = xmldict['svg']['@width']
    total_height = xmldict['svg']['@height']
    outdict['svg'] = {}
    outdict['svg']['width'] = xmldict['svg']['@width']
    outdict['svg']['height'] = xmldict['svg']['@height']
    outdict['svg']['path'] = []
    for path in xmldict['svg']['path']:
        boxpath = {}
        oldpath = {}
        matcher = re.match("^(.*fill:#)(.+?)(;.*$)", path['@style'])
        if matcher.group(0) is not None:
            fill = matcher.group(2)
            avg_col = average_color(fill)
            if avg_col < 0x50:
                boxpath['@style'] = matcher.group(1) + '000000' + matcher.group(3)
                oldpath['@style'] = matcher.group(1) + 'FF0000' + matcher.group(3)
        polygon = path_to_polygon(path['@d'])
        oldpath['@d'] = polygon_to_path(polygon)
        next = bounding_box(polygon)
        if threshold_area(next, 0.04):
            boxpath['@d'] = polygon_to_path(next)
        if ('@style' in boxpath) and ('@d' in boxpath):
            outdict['svg']['path'].append(boxpath)
            outdict['svg']['path'].append(oldpath)            
    outf = open(outputfile,'w')
    outf.write(xmltodict.unparse(outdict, pretty=True))
    outf.close()
    
def path_to_polygon(path):
    polygon = []
    nodes = re.findall('([MLHVCSQTAZ][\d\.\s]+)', path, re.I)
    for n in nodes:
        movematch = re.match('M(.+?) (.+?)$', n, re.I)
        linematch = re.match('L(.+?) (.+?)$', n, re.I)
        curvematch = re.match('C(.+?) (.+?) (.+?) (.+?) (.+?) (.+?)$', n, re.I)
        if curvematch is not None:
            polygon.append('%s %s' % (curvematch.group(5), curvematch.group(6)))
        elif movematch is not None:
            polygon.append('%s %s' % (movematch.group(1), movematch.group(2)))
        elif linematch is not None:
            polygon.append('%s %s' % (linematch.group(1), linematch.group(2)))
    return polygon

def polygon_to_path(polygon):
    return 'M' + 'L'.join(polygon) + 'Z'

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
#     print str(int(max_x)) + '-' + str(int(min_x)) +'='+ str(abs((int(max_x) - int(min_x))))
    print str(max_y) + " " + str(min_y)
#     print str(abs(float(max_y - min_y) / float(total_height))) + ">" + str(threshold)
    if abs(float(max_y - min_y) / float(total_height)) > float(threshold):
        print "y"
        return True
    if abs(float(max_x - min_x) / float(total_width)) > float(threshold):
        print "x"
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

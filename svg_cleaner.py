#!/usr/bin/env python

import re
import xmltodict
import sys
import json
import os
from svg.path import Path, Line, Arc, CubicBezier, QuadraticBezier, parse_path

total_width = 0
total_height = 0
scale_width = 1.0
scale_height = 1.0
max_x = 0
max_y = 0
min_x = 0
min_y = 0

def main():
    filename = sys.argv[1]
    numtaxa = sys.argv[2]
    # potrace -o outputfile -s -k 0.8 -W 10 -H 10 raw_pbm_file
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
    global total_width, total_height, scale_width, scale_height
    total_width = float(xmldict['@width'].replace('pt',''))
    total_height = float(xmldict['@height'].replace('pt',''))
    xmlpaths = []
    style = {}
    transform = ''
    
    # parse original paths
    if 'path' in xmldict:
        xmlpaths = xmldict['path']
    elif 'g' in xmldict:
        xmlpaths = xmldict['g']['path']
        del xmldict['g']['path']
        if '@transform' in xmldict['g']:
            transform = xmldict['g']['@transform']
            if transform != '':
                parse_transform(transform)
    polygons = []
    for path in xmlpaths:
        polygon = path_to_polygon(path['@d'])
        if threshold_area(bounding_box(polygon), 0.6):
            polygons.append(polygon)

    global max_x, max_y, min_x, min_y
    global scale_width, scale_height
    max_x = abs(max_x * scale_width)
    max_y = abs(max_y * scale_height)
    for polygon in polygons:
        polygon = scale(polygon, scale_width, scale_height)
        polygon = translate(polygon, 0, max_y)

    lines = []
    segments = []    
    radius = max_y / (int(numtaxa)*5)
    polygon_total = []
    for polygon in polygons:
        polygon = cleanup_polygon(polygon, radius)
        polygon = simplify_polygon(polygon)
        polygon = cleanup_polygon(polygon, radius*1.8)
        polygon = simplify_polygon(polygon)
        segments.extend(lineify_polygon(polygon))
        segments = remove_dups(segments)
        lines.extend(segments_to_lines(segments))
#         polygon_total.extend(polygon)
#         circles.extend(polygon_to_circles(polygon))    
#         path = {}
#         path['@d'] = polygon_to_path(polygon)
#         path['@style'] = "fill:#FF0000; stroke:#FF0000;"
#         paths.append(path) 
    
    make_tree(segments)
    outdict['svg'] = {}
    outdict['svg']['width'] = xmldict['@width']
    outdict['svg']['height'] = xmldict['@height']
    outdict['svg']['g'] = [{'line':lines}]

    outf = open(outputfile,'w')
    outf.write(xmltodict.unparse(outdict, pretty=True))
    outf.close()

def scale(polygon, scale_width, scale_height):
    for node in polygon:
        node[0] = int(float(node[0]) * scale_width)
        node[1] = int(float(node[1]) * scale_height)
    return polygon
    
def translate(polygon, x, y):
    for node in polygon:
        node[0] = int(float(node[0]) + x)
        node[1] = int(float(node[1]) + y)
    return polygon

def parse_transform(transform):
    global scale_width, scale_height
    scalematcher = re.search('scale\(([0-9\-\.]+),*(.*?)\)', transform)
    if scalematcher is not None:
        scale_width = float(scalematcher.group(1))
        if scalematcher.group(2) is not None:
            scale_height = float(scalematcher.group(2))
            
def make_tree(segments):
    vert_lines = set()
    horiz_lines = set()
    nodes = set()
    levels = set()   
    for seg in segments:
        nodes.add(int(seg[1]))
        nodes.add(int(seg[3]))
        levels.add(int(seg[0]))
        levels.add(int(seg[2]))
        seg_str = '%s %s %s %s' % (seg[0], seg[1], seg[2], seg[3])
        if seg[0] == seg[2]:
            vert_lines.add(seg_str)
        elif seg[1] == seg[3]:
            horiz_lines.add(seg_str)
    
    # how many different levels are there?
    levels = list(levels)
    # sort them backwards because we want to work from the leaves back
    levels.sort(cmp=lambda x,y: cmp(int(y), int(x)))
    print "there are %s levels" % levels
    # how many different nodes are there?
    nodes = list(nodes)
    nodes.sort(cmp=lambda x,y: cmp(int(x), int(y)))
    print "there are %s nodes" % nodes        
#     
#     result_segments = []
#     # starting from the leaves
#     for i in range(len(levels)-1):
#         out_x = levels[i]
#         in_x = levels[i+1]
#         # what are the node-lines that represent nodes at this level?
#         node_lines = set()
#         for line in vert_lines:
#             if 
                
            
        
def remove_dups(segments):
    seg_set = set()
    vert_set = set()
    for seg in segments:
        # clean up horizontal lines
        if seg[1] == seg[3]:
            if int(seg[0]) < int(seg[2]):
                seg_set.add('%d %d %d %d' % (seg[0], seg[1], seg[2], seg[3]))
            else:
                seg_set.add('%d %d %d %d' % (seg[2], seg[1], seg[0], seg[3]))
        # clean up vertical lines
        if seg[0] == seg[2]:
            if int(seg[1]) < int(seg[3]):
                vert_set.add('%d %d %d %d' % (seg[0], seg[1], seg[2], seg[3]))
            else:
                vert_set.add('%d %d %d %d' % (seg[0], seg[3], seg[2], seg[1]))
    
    # remove the smaller vertical lines, if they overlap with a longer line.
    sorted_verts = list(vert_set)
    sorted_verts.sort(cmp=lambda x,y: cmp(y, x))
    pruned_verts = []
    
    seg_list = []
    for seg in seg_set:
        seg_list.append(re.split(' ',seg))
    for seg in vert_set:
        seg_list.append(re.split(' ',seg))
    return seg_list
    
def segments_to_lines(segments):
    lines = []
    for seg in segments:
        lines.append({'@x1':str(seg[0]), '@y1':str(seg[1]), '@x2':str(seg[2]), '@y2':str(seg[3]), '@stroke-width':'5', '@stroke':'green'})
    return lines

def lineify_polygon(polygon):
    lines = []
    segment_hash = {}
    last_node = polygon[0]
    for node in polygon:
        lines.append([last_node[0], last_node[1], node[0], node[1]])
        last_node = node
    return lines
       
def cleanup_polygon(polygon, radius):
    # find all the horizontal points
    x_sort_dict = {}
    x_sort_points = [x for x in polygon]
    x_sort_points.sort(cmp=lambda x,y: cmp(float(x[0]), float(y[0])))
    
    curr_x = 0
    for point in x_sort_points:
        if float(point[0]) > (float(curr_x) + float(radius)):
            curr_x = point[0]
        x_sort_dict['%d %d' % (point[0],point[1])] = [curr_x,point[1]]

    new_polygon = []
    for point in polygon:
        new_polygon.append(x_sort_dict['%d %d' % (point[0],point[1])])
        
    # find all the vertical points
    sort_dict = {}
    y_sort_points = [x for x in new_polygon]
    y_sort_points.sort(cmp=lambda x,y: cmp(float(x[1]), float(y[1])))
    curr_y = 0
    for point in y_sort_points:
        if float(point[1]) > (float(curr_y) + float(radius)):
            curr_y = point[1]
        sort_dict['%d %d' % (point[0],point[1])] = [point[0],curr_y]
    polygon = []
    for point in new_polygon:
        polygon.append(sort_dict['%d %d' % (point[0],point[1])])
    return polygon   

# remove all in-between singletons from a cleaned-up polygon
def simplify_polygon(polygon):
    new_polygon = [polygon[0]]
    for i in range(len(polygon)-2):
        add_me = True
        # if the three y-vals are equal
        if (polygon[i][1] == polygon[i+1][1]) and (polygon[i+1][1] == polygon[i+2][1]):
            # if polygon[i+1][0] is between polygon[i][0] and polygon[i+2][0], do not add
            if (polygon[i][0] < polygon[i+1][0] and polygon[i+1][0] < polygon[i+2][0]) or (polygon[i][0] > polygon[i+1][0] and polygon[i+1][0] > polygon[i+2][0]):
                add_me = False
        # if the three x-vals are equal
        if (polygon[i][0] == polygon[i+1][0]) and (polygon[i+1][0] == polygon[i+2][0]):
            # if polygon[i+1][1] is between polygon[i][1] and polygon[i+2][1], do not add
            if (polygon[i][1] < polygon[i+1][1] and polygon[i+1][1] < polygon[i+2][1]) or (polygon[i][1] > polygon[i+1][1] and polygon[i+1][1] > polygon[i+2][1]):
                add_me = False
        if add_me == True:
            new_polygon.append(polygon[i+1])
    return new_polygon

def path_to_polygon(path):
    polygon = []
    global max_x, max_y, min_x, min_y
    new_path = Path()    
    for segment in parse_path(path):
        new_path.append(Line(segment.start, segment.end))
    new_path.closed = True
    nodes = re.findall('[ML]\s*(\d+\.*\d*,\d+\.*\d*)\s*', new_path.d())
    for n in nodes:
        coords = n.split(',')
        if max_x < int(coords[0]):
            max_x = int(coords[0])
        if max_y < int(coords[1]):
            max_y = int(coords[1])
        if min_x > int(coords[0]):
            min_x = int(coords[0])
        if min_y > int(coords[1]):
            min_y = int(coords[1])
        polygon.append([int(coords[0]), int(coords[1])])
    return polygon
    
def polygon_to_path(polygon):
    path_points = []
    for point in polygon:
        path_points.append('%d %d' % (point[0],point[1]))
    return 'M' + 'L'.join(path_points) + 'Z'
          
def polygon_to_circles(polygon):
    circlelist = []
    for i in range(len(polygon)):
        circledict = {}
        coords = polygon[i]
        circledict['@r'] = '3'
        circledict['@stroke'] = 'none'
        circledict['@fill'] = '#c6c6%(number)02x' % {"number":(i)}
        circledict['@cx'] = str(coords[0])
        circledict['@cy'] = str(coords[1])
        circlelist.append(circledict)
    return circlelist

def average_color(col):
    matcher = re.match("(..)(..)(..)", col)
    r = matcher.group(1)
    g = matcher.group(2)
    b = matcher.group(3)
    average = (int(r,16) + int(g,16) + int(b,16))/3
    return average

def threshold_area(rect, threshold):
    mins = rect[0]
    maxs = rect[2]
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
    x_points = []
    y_points = []
    for point in polygon:
        coord = point
        x_points.append(coord[0])
        y_points.append(coord[1])
    max_x = max(x_points)
    max_y = max(y_points)
    min_x = min(x_points)
    min_y = min(y_points)
    return [[min_x, min_y],[min_x, max_y],[max_x, max_y],[max_x, min_y]]
    
if __name__ == '__main__':
    main()

#!/usr/bin/env python

import re
import xmltodict
import sys
import json
import os
import math
from svg.path import Path, Line, Arc, CubicBezier, QuadraticBezier, parse_path

total_width = 0
total_height = 0
scale_width = 1.0
scale_height = 1.0
max_x = 0
max_y = 0
min_x = 0
min_y = 0
otu_level = 0

points = []

def main():
    global points
    filename = sys.argv[1]
    # convert input_image -threshold 50% raw.pbm
    # potrace -s -k 0.8 -W 10 -H 10 -o output.svg raw.pbm
    file_name, extension = os.path.splitext(sys.argv[1])
    outfile = file_name
    if extension != '.svg':
        print "can't open this file"
        return
    f = open(filename, 'r')
    xmldict = ''
    xml = f.read()
    xmldict = xmltodict.parse(xml)['svg']
    global total_width, total_height, scale_width, scale_height
    s = re.sub('[a-zA-Z]+', '', xmldict['@width'])
    total_width = float(s)
    s = re.sub('[a-zA-Z]+', '', xmldict['@height'])
    total_height = float(s)
    xmlpaths = []
    style = {}
    transform = ''
    paths = [] 
    circles = []

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
    raw_polygons = []
    for path in xmlpaths:
        polygon = path_to_polygon(path['@d'])
        raw_polygons.append(polygon)

    global max_x, max_y, min_x, min_y
    global scale_width, scale_height
    max_x = abs(max_x * scale_width)
    max_y = abs(max_y * scale_height)
    polygons = []
    for polygon in raw_polygons:
        polygon = scale(polygon, scale_width, scale_height)
        polygon = translate(polygon, 0, max_y)
        if threshold_area(bounding_box(polygon), 0.6):
            path = {}
            path['@d'] = nodes_to_path(polygon)
            path['@style'] = "fill:#CCCCCC; stroke:#999999; stroke-width:1"
            paths.append(path) 
            polygons.append(polygon)
        else:
            path = {}
            path['@d'] = nodes_to_path(polygon)
            path['@style'] = "fill:#EEEEEE; stroke:#EEEEEE; stroke-width:1"
#             paths.append(path) 

    segments = []   
    polygon_total = []
    for polygon in polygons:
        polygon = straighten_polygon(polygon)
        polygon = even_out_polygon(polygon)
        # this path is for the cleaned-up lines
        path = {}
        path['@d'] = nodes_to_path(polygon)
        path['@name'] = "cleaned path"
        path['@style'] = "fill:none; stroke:#FF0000; stroke-width:1"
        paths.append(path) 

    segments.extend(normalize_polygon_to_lines(polygon))
    # generate raw svg first-pass, in case something fails during tree building:
    circles.extend(nodes_to_circles(points))
    lines = []
    lines.extend(segments_to_lines(segments, 'blue', 4))
    svgdict = {}
    svgdict['svg'] = {}
    svgdict['svg']['@width'] = xmldict['@width']
    svgdict['svg']['@height'] = xmldict['@height']
    svgdict['svg']['g'] = [{'line':lines, 'path':paths, 'circle':circles}]

    outf = open(outfile+'_raw.svg','w')
    outf.write(xmltodict.unparse(svgdict, pretty=True))
    outf.close()
    
    print "printed raw svg"
    
    # make the raw tree for making nexml:
    (nodes, edges, otus) = make_tree(segments)
    
    # generate nexml:
    nodedict = {}
    otudict = {}
    index = 1
    for otu in otus:
        otudict[str(otu)] = 'otu%d' % index
        index = index+1

    nodes.extend(otus)
    index = 1  
    for node in nodes:
        nodedict[str(node)] = 'node%d' % index
        index = index+1
    
    nexmldict = {}
    nexmldict['nex:nexml'] = {'@xmlns:nex':'http://www.nexml.org/2009'}
    nexmldict['nex:nexml']['@xmlns']="http://www.nexml.org/2009"
    nexmldict['nex:nexml']['@xmlns:rdf']="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
    nexmldict['nex:nexml']['@xmlns:xsd']="http://www.w3.org/2001/XMLSchema#"
    nexmldict['nex:nexml']['@xmlns:xsi']="http://www.w3.org/2001/XMLSchema-instance" 
    nexmldict['nex:nexml']['@version'] = '0.9'
    nexmldict['nex:nexml']['otus'] = {'@about':'#otus', '@id':'otus','@label':'taxa'}
    nexmldict['nex:nexml']['otus']['otu'] = []
    for otu in otus:
        nexml_otu = {'@id':otudict[str(otu)]}
        nexml_otu['@about'] = '#%s' % otudict[str(otu)]
        nexml_otu['@label'] = otudict[str(otu)]
        nexmldict['nex:nexml']['otus']['otu'].append(nexml_otu)
    
    nexmldict['nex:nexml']['trees']= {'@about':'#trees1','@id':'trees1','@label':'trees','@otus':'otus'}
    currtree = {'@id':'tree1', '@about':'#tree1', '@label':'tree', '@xsi:type':'nex:FloatTree'}
    nexmldict['nex:nexml']['trees']['tree'] = [currtree]
    currtree['node'] = []
    for node in nodes:
        nexml_node = {'@id':nodedict[str(node)]}
        if str(node) in otudict:
            nexml_node['@otu'] = otudict[str(node)]
        currtree['node'].append(nexml_node)   

    nexmldict['nex:nexml']['trees']['edge'] = []
    index = 1
    currtree['edge'] = []
    for edge in edges:
        nexml_edge = {}
        nexml_edge['@id'] = 'edge%d' % index
        nexml_edge['@length'] = str(edge[2]-edge[0])
        nexml_edge['@source'] = nodedict[str([edge[0],edge[1]])]
        nexml_edge['@target'] = nodedict[str([edge[2],edge[3]])]
        currtree['edge'].append(nexml_edge)  
        index = index+1 

    outf = open(outfile+'.xml','w')
    outf.write(xmltodict.unparse(nexmldict, pretty=True))
    outf.close()
    
    # generate nexus:
    nexus_str = "#NEXUS\n"
    nexus_str += "begin TREES;\n"
    nexus_str += "Tree tree=\n"
    nexus_str += tree_to_nexus(otus, nodes, edges)
    nexus_str += ";\nEnd;\n"
    
    outf = open(outfile+'.nex','w')
    outf.write(nexus_str)
    outf.close()
                

    textdict = {}
    textdict = []
    # for each otu, add a text label:
    for k in otudict.keys():
        coordmatch = re.match('\[(\d+), (\d+)\]',k)
        if coordmatch is not None:
            x = coordmatch.group(1)
            y = coordmatch.group(2)
            textnode = {'@x':str(int(x) + 10), '@y':y,'@fill':'black','#text':otudict[k]}
            textdict.append(textnode)
            
    # generate svg:
    circles.extend(nodes_to_circles(nodes))
    lines.extend(segments_to_lines(edges, "green", 3))
    svgdict = {}
    svgdict['svg'] = {}
    svgdict['svg']['@width'] = xmldict['@width']
    svgdict['svg']['@height'] = xmldict['@height']
    svgdict['svg']['g'] = [{'path':paths, 'line':lines, 'circle':circles, 'text':textdict}]
    

    outf = open(outfile+'_raw.svg','w')
    outf.write(xmltodict.unparse(svgdict, pretty=True))
    outf.close()
    print "reprinted raw svg"

def tree_to_nexus(otus, nodes, edges):
    otudict = {}
    nodedict = {}
    for i in range(len(otus)):
        nodedict[str(otus[i])] = str('otu%d' % i)
    # find the root: it's the node that is not the y2 of any edge
    root = None
    for node in nodes:
        targetnum = 0
        sourcenum = 0
        if str(node) not in nodedict:
            nodedict[str(node)] = []
        for edge in edges:
            if edge[1] == node[1] and edge[0] == node[0]:
                targetnum += 1
                edgenex = str([edge[2],edge[3]])
                nodedict[str(node)].append(edgenex)
            if edge[3] == node[1] and edge[2] == node[0]:
                sourcenum += 1
        if sourcenum == 0:
            root = node
    return replace_nodes(nodedict, str(root))

def replace_nodes(nodedict, newick):
    # find the nodes in newick:
    nodematcher = re.findall('\[\d+, \d+\]',newick)
    if len(nodematcher) is 0:
        return newick
    else:
        for node in nodematcher:
            # replace the node with the children of the node
            child = str(nodedict[node])
            if '[' in child:
                child = '(' + str(', '.join(nodedict[node])) + ')'
            newick = newick.replace(node, child)                
    return replace_nodes(nodedict, newick)

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
    vert_lines = []
    horiz_lines = []
    for seg in segments:
        if seg[0] == seg[2]:
            vert_lines.append(seg)
        elif seg[1] == seg[3]:
            horiz_lines.append(seg)
    
    # create a dictionary of nodes: each node has vertical endpoints and is indexed by its x value
    node_dict = {}
    sorted_verts = sort_lines(vert_lines, 0, 1)

    for bin in sorted_verts:
        if bin[0][0] not in node_dict:
            node_dict[bin[0][0]] = []
        for line in bin:
            x = line[0]
            y1 = line[1]
            y2 = line[3]
            node_dict[x].append([y1,y2])
        
    # okay, now we know what the nodes are. Match up the edges.
    edges = set()
    otus = []
    nodes = set()
    for line in horiz_lines:
        x1 = line[0]
        x2 = line[2]
        y1 = line[1]
        y2 = line[3]
        # we want to make the y1 equal to the y1 of the node
        # look for the node that this x1 is in:
        if x1 in node_dict:
            for node in node_dict[x1]:
                if y1 >= node[0] and y1 <= node[1]:
                    y1 = node[0]
            if x2 not in node_dict:
                otus.append([x2, y2])
            else:
                for node in node_dict[x2]:
                    if y2 >= node[0] and y2 <= node[1]:
                        y2 = node[0]
                nodes.add('%d %d' % (x1, y1))
                nodes.add('%d %d' % (x2, y2))
            edges.add('%d %d %d %d' % (x1, y1, x2, y2))

    final_nodes = []
    for node in nodes:
        coords = re.split(' ',node)
        final_nodes.append([int(coords[0]), int(coords[1])])
    final_edges = []
    for edge in edges:
        coords = re.split(' ',edge)
        final_edges.append([int(coords[0]), int(coords[1]), int(coords[2]), int(coords[3])])
    
    otus.sort(cmp=lambda x,y: cmp(x[1], y[1]))
    
    return (final_nodes, final_edges, otus)

def segments_to_lines(segments, color, width):
    lines = []
    for seg in segments:
        lines.append({'@x1':str(seg[0]), '@y1':str(seg[1]), '@x2':str(seg[2]), '@y2':str(seg[3]), '@stroke-width':str(width), '@stroke':color})
    return lines

def normalize_polygon_to_lines(polygon):    
    global max_x, min_x
    lines = lineify_path(polygon)

    root_level = max_x
    otu_level = min_x
    horiz_line_set = set()
    vert_line_set = set()
    for line in lines:
        forward_line = (line[0],line[1],line[2],line[3])
        if (line[0] == line[2]): # vertical lines
            if (line[1] < line[3]):
                forward_line = (line[0],line[3],line[2],line[1])
            else:
                vert_line_set.add('%d %d %d %d' % forward_line)
        elif (line[1] == line[3]): #horiz line
            if (line[0] > line[2]):
                forward_line = (line[2],line[1],line[0],line[3])
            else:
                horiz_line_set.add('%d %d %d %d' % forward_line)
            
            if (forward_line[0] < root_level):
                root_level = forward_line[0]
            if (forward_line[2] > otu_level):
                otu_level = forward_line[2]
    
    # sort all the horiz lines by y-value:
    horiz_lines_list = []
    for line in horiz_line_set:
        coord = line.split(' ')
        x1 = int(coord[0])
        y1 = int(coord[1])
        x2 = int(coord[2])
        y2 = int(coord[3])
        
        horiz_lines_list.append([x1, y1, x2, y2])
        
    horiz_lines_binned = sort_lines(horiz_lines_list, 1, 0)

    vert_lines_list = []
    # adjust the y2 of the vertical lines to match a horizontal line
    for line in vert_line_set:
        coord = line.split(' ')
        x = int(coord[0])
        y1 = int(coord[1])
        y2 = int(coord[3])

        my_nodeline = None
        for bin in horiz_lines_binned:
            if (y1 > bin[0][1]):
                continue
            else:
                for nodeline in bin:
                    if (nodeline[0] <= x) and (x <= nodeline[2]):
                        my_nodeline = nodeline
                        y1 = nodeline[1]
                        break
                if my_nodeline is not None:
                    break
        line = [x, y1, x, y2]

        vert_lines_list.append(line)

    # sort all the vertical lines by x-value:
    vert_lines_binned = sort_lines(vert_lines_list, 0, 1)

    # coalesce the vertical nodes
    for i in range(len(vert_lines_binned)):
        bin = vert_lines_binned[i]
        new_bin = []
        curr_node = bin[0]
        for j in range(1,len(bin)):
            line = bin[j]
            if (curr_node[1] == line[3]):
                curr_node[1] = line[1]
            else:
                new_bin.append(curr_node)
                curr_node = bin[j]
        new_bin.append(curr_node)
        vert_lines_binned[i] = new_bin
        
    lines = []
    for bin in vert_lines_binned:
        for line in bin:
            lines.append([line[0],line[3],line[2],line[1]])
    
    for line in horiz_line_set:
        coord = line.split(' ')
        x1 = int(coord[0])
        y = int(coord[1])
        x2 = int(coord[2])

        # if x2 is at otu_level, we don't need to worry about that end.
        if x2 < otu_level:        
            my_nodeline = None
            for bin in vert_lines_binned:
                if x2 > bin[0][0]:
                    continue
                else:
                    for nodeline in bin:
                        if (nodeline[3] <= y) and (y <= nodeline[1]):
                            my_nodeline = nodeline
                            x2 = nodeline[0]
                            break
                    if my_nodeline is not None:
                        break
        # if x1 is at root_level, we don't need to worry about that end.
        if x1 > root_level:    
            my_nodeline = None
            i = 0
            while i < len(vert_lines_binned):
                bin = vert_lines_binned[i]
                if (x1 > bin[0][0]):
                    i += 1
                    continue
                else:                    
                    for nodeline in vert_lines_binned[i]:
                        if (nodeline[3] <= y) and (y <= nodeline[1]):
                            my_nodeline = nodeline
                            x1 = nodeline[0]
                            break
                    if my_nodeline is not None:
                        break
                    i += 1
            line = [x1,y,x2,y]

        lines.append([x1, y, x2, y])
        
    return lines

def sort_lines(lines, key1, key2):
    bin_by_x1 = {}
    for line in lines:
        if line[key1] not in bin_by_x1:
            bin_by_x1[line[key1]] = []
        bin_by_x1[line[key1]].append(line)
    
    bin_keys = list(bin_by_x1.keys())
    bin_keys.sort()
    
    final_lines = []
    for bin in bin_keys:
        bin_by_x1[bin].sort(cmp=lambda x,y: cmp(x[key2], y[key2]))
        final_lines.append(bin_by_x1[bin])
    return final_lines


def lines_to_polygon(lines):
    polygon = []
    coord = lines[0].split(' ')
    polygon.append([int(coord[0]),int(coord[1])])
    for line in lines:
        coord = line.split(' ')
        polygon.append([int(coord[2]),int(coord[3])])
    return polygon
    

def lineify_path(polygon):
    lines = []
    lines.append([polygon[len(polygon)-1][0], polygon[len(polygon)-1][1], polygon[0][0], polygon[0][1]])
    last_node = polygon[0]
    for i in range(1, len(polygon)-1):
        node = polygon[i]
        lines.append([last_node[0], last_node[1], node[0], node[1]])
        last_node = node
    return lines
       
# remove all in-between singletons from a cleaned-up polygon
def straighten_polygon(polygon):
    # for convenience:
    x = 0
    y = 1
    global points
    points = []
    changes_made = False
    # we need to make sure we start with the last thing in polygon
    new_polygon = []
    new_polygon.insert(0,polygon.pop())
    new_polygon.append(polygon.pop(0))
    new_polygon.append(polygon.pop(0))
    
    tips = 0
    
    while len(polygon) >= 0:
        node3 = new_polygon.pop()
        node2 = new_polygon.pop()

        if len(new_polygon) > 0:
            node1 = new_polygon.pop()
        else:
            node1 = polygon.pop()
        
        keep_node = True

        #### Remove duplicate or near-dup points.
        # if node1 and node2 are peculiarly close together on the y-axis, we should drop node 2 and try again:
        if (node1[x] == node2[x]):
            if ((node2[y] - 2 <= node1[y]) and (node1[y] <= node2[y] + 2)) or ((node1[y] - 2 <= node2[y]) and (node2[y] <= node1[y] + 2)):
                polygon.insert(0,node3)
                new_polygon.append(node1)
                new_polygon.append(polygon.pop(0))
                continue
        
        if (node2[y] != node1[y]) and (node2[x] != node1[x]):
            if (math.fabs(node2[y] - node1[y]) <= 2): # so small that angles can't detect them
                node1[y] = node2[y] 
                changes_made = True
            if (math.fabs(node2[x] - node1[x]) <= 2): # so small that angles can't detect them
                node1[x] = node2[x] 
                changes_made = True
        
        #### FIRST: normalize the tips
        # if node2[x] is greater than either node1[x] or node3[x]
        if ((node3[x] < node2[x]) and (node1[x] < node2[x])) or ((node3[x] > node2[x]) and (node1[x] > node2[x])):
            node1[y] = node2[y]
            node3[y] = node2[y]
            tips += 1

        #### SECOND: normalize knees
        # if node2[x] is between node1[x] and node3[x]
        if ((node1[x] <= node2[x]) and (node2[x] <= node3[x])) or ((node3[x] <= node2[x]) and (node2[x] <= node1[x])):
            theta = math.degrees(math.atan2(math.fabs((node3[y]-node2[y])),math.fabs((node3[x]-node2[x]))))
            if (node2[x]==node1[x]): # vertical, so snap node 3's x into line
                if (theta >= 45): # theta == 90 means this is a straight knee
                    node3[x] = node2[x]
                    keep_node = False
                    if (theta != 90):
                        changes_made = True                
            elif (node2[y]==node1[y]): #horizontal
                if (theta <= 45): # theta == 0 means this is a straight knee
                    node3[y] = node2[y]
                    keep_node = False
                    if (theta != 0):
                        changes_made = True
            else:
                # calculate an angle between n1 and n2:
                theta = math.degrees(math.atan2(math.fabs((node1[y]-node2[y])),math.fabs((node1[x]-node2[x]))))
                if (node3[x] == node2[x]): # vertical, so snap node 3's x into line
                    if (theta >= 45): # theta == 90 means this is a straight knee
                        node1[x] = node2[x]
                        if (theta != 90):
                            changes_made = True
                elif (node2[y]==node3[y]): #horizontal
                    if (theta <= 45): # theta == 0 means this is a straight knee
                        node1[y] = node2[y]
                        if (theta != 0):
                            changes_made = True
        #### FINALLY: append nodes, without node2 if it's a straight knee
        new_polygon.append(node1)
        if keep_node:
            new_polygon.append(node2)
        new_polygon.append(node3)
        if len(polygon) == 0:
            break
        
        new_polygon.append(polygon.pop(0))
    
    if changes_made:
        new_polygon = straighten_polygon(new_polygon)

    return new_polygon
    
def even_out_polygon(polygon):
    # for convenience:
    x = 0
    y = 1
    polygon.insert(0,polygon[len(polygon)-1])
    polygon.insert(0,polygon[len(polygon)-2])
    polygon.insert(0,polygon[len(polygon)-3])

    global otu_level
    # find the maximum x-value
    for node in polygon:
        if (node[x] > otu_level):
            otu_level = node[x]
    
    global points
    points = []
    changes_made = False
    # we need to make sure we start with the last thing in polygon
    new_polygon = []
    new_polygon.insert(0,polygon.pop())
    new_polygon.append(polygon.pop(0))
    new_polygon.append(polygon.pop(0))
    new_polygon.append(polygon.pop(0))
    new_polygon.append(polygon.pop(0))
        
    while len(polygon) >= 0:
        node4 = new_polygon.pop()
        node3 = new_polygon.pop()
        node2 = new_polygon.pop()
        
        if (node3 == node2):
            node2 = new_polygon.pop()
                
        node1 = new_polygon.pop()
        if len(new_polygon) > 0:
            node0 = new_polygon.pop()
        else:
            node0 = polygon.pop()
        

        ### if node1[y] == node2[y] and node0[x] < node1[x] and node3[x] < node2[x]
        if node1[x] == node2[x] and node0[x] < node1[x] and node3[x] < node2[x]:
            node2[y] = node1[y]
            node3[y] = node1[y]            
        
        #### FIRST: normalize the tips
        # if node2[x] is greater than either node1[x] or node3[x]
        if (node1[y] == node2[y]) and (node2[y] == node3[y]):
            if ((node3[x] < node2[x]) and (node1[x] < node2[x])):
                new_x = node1[x]
                if (node3[x] > node1[x]):
                    new_x = node3[x]
                node0 = [new_x, node0[y]]
                node1 = [new_x, node2[y]]
                node2 = [otu_level, node2[y]]
                node3 = [new_x, node2[y]]
                node4 = [new_x, node4[y]]

        #### FINALLY: append nodes
        new_polygon.append(node0)
        new_polygon.append(node1)
        new_polygon.append(node2)
        new_polygon.append(node3)
        new_polygon.append(node4)
        if len(polygon) == 0:
            break
        
        new_polygon.append(polygon.pop(0))
    
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
    polygon.pop()
    return polygon
    
def nodes_to_path(nodes):
    path_points = []
    for point in nodes:
        path_points.append('%d %d' % (point[0],point[1]))
    return 'M' + 'L'.join(path_points) + 'Z'
          
def nodes_to_circles(nodes):
    circlelist = []
    for i in range(len(nodes)):
        circledict = {}
        coords = nodes[i]
        circledict['@r'] = '3'
        circledict['@stroke'] = 'black'
        circledict['@stroke-width'] = '1'
        circledict['@fill'] = 'yellow'
        circledict['@cx'] = str(coords[0])
        circledict['@cy'] = str(coords[1])
        circlelist.append(circledict)
    return circlelist

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

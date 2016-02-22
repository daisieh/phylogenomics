#!/usr/bin/env python

import re
import xmltodict
import sys
import json
import os

simplenodelist = []
otudict = {}
edgelookupdict = {}
otulookupdict = {}
root = ''
publication = ''
def main():
    filename = sys.argv[1]
    file_name, extension = os.path.splitext(sys.argv[1])
    outputfile = '%s.txt' % filename
    f = open(filename, 'r')
    xmldict = ''
    if (re.match('\.(ne)*xml$',extension)):
        xml = f.read()
        xmldict = xmltodict.parse(xml)['nex:nexml']
    else:
        print "can't parse this file"
        return;
    edgelist = xmldict['trees']['tree']['edge']
    nodelist = xmldict['trees']['tree']['node']
    otulist = []
    if 'otus' in xmldict:
        otulist = xmldict['otus']['otu']
    for otu in otulist:
        otulookupdict[otu['@id']] = otu['@label'].replace(' ','_')
    for x in xmldict['meta']:
        if '@property' in x:
            if x['@property'] == 'ot:studyPublicationReference':
                publication = x['#text']
    for node in nodelist:
        if '@otu' in node:
            otudict[node['@id']] = otulookupdict[node['@otu']]
        simplenodelist.append(node['@id'])
        if '@root' in node:
            root = node['@id']
    for edge in edgelist:
        source = edge['@source']
        target = edge['@target']
        if source not in edgelookupdict:
            edgelookupdict[source] = []
            
        edgelookupdict[source].append(target) 
    outf = open(outputfile, 'w')
    result = '%s\n%s' % (parse_pub_to_taxonomy(publication),make_clade_with_node(root))
    outf.write(result)
    outf.close()

def make_clade_with_node(node):
    result = ''
    if node not in edgelookupdict:
        return
    in_clade = edgelookupdict[node]
    named_clade = []
    for clade_member in in_clade:
        if clade_member in otudict:
            named_clade.append(otudict[clade_member])
        else:
            named_clade.append(clade_member)
    result += ('(%s %s)' % (node, ' '.join(named_clade)))
    for next_node in in_clade:
        res = make_clade_with_node(next_node)
        if res is not None:
            result += '\n' + res
    return result
        
def parse_pub_to_taxonomy(pubstring):
    result = ''
    authors = ''
    year = ''
    matcher = re.match("(.+)\s*,\s*(\d+).*", pubstring)
    # A. E. Marvaldi, R. G. Oberprieler, C. H. C. Lyal, T. Bradbury, R. S. Anderson, 2006, ' Phylogeny of the Oxycoryninae sensu lato (Coleoptera:Belidae) and evolution of host-plant associations ', Invertebrate Systematics, vol. 20, no. 4, p. 447
    if matcher is not None:
        authors = matcher.group(1)
        year = matcher.group(2)
    matcher = re.match("(.+)\s*\((\d+)\).*", pubstring)
    # Vanin, S.A. (1976) Taxonomic revision of the South American Belidae (Coleoptera). Arquivos de Zoologia, Vol. 28: 1-75
    if matcher is not None:
        authors = matcher.group(1)
        year = matcher.group(2)
    authorlist = re.split(',', authors)
    lastname = max(re.split('\s', authorlist[0]))
    
    result += ("taxonomy %s %s" % (year, lastname))
    return result
    
if __name__ == '__main__':
    main()

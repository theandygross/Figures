__author__ = 'agross'

import re
import itertools as itertools
import urllib

import pandas as pd
from matplotlib.colors import rgb2hex
from matplotlib.cm import RdBu

KEGG_PATH = 'http://www.kegg.jp/kegg-bin/'


from Figures.KEGG import *
def pull_pathway_info_from_kegg(kegg_id):
    o = urllib.urlopen('http://rest.kegg.jp/get/' + kegg_id).read()
    o = o.splitlines()

    '''need to parse out when new sections start'''
    sections = {}
    for i, n in enumerate(o):
        s = n.split()[0]
        if n[0] != ' ' and s not in sections:
            sections[s] = i
    sections = pd.Series(sections).order()
    o = [l[12:] for l in o]  # get rid of fixed-width section headings

    '''Pull out gene information, ec = enzyme, ko = complex'''
    start = sections['GENE']
    stop = [sections.iloc[i+1] for i, s in enumerate(sections.index)
            if s == 'GENE'][0]
    gene = o[start:stop]
    #return gene
    mapping = []
    for g in gene:
        try:
            g = g.split(';')
            gene_id = g[0].split()[0]
            gene_name = g[0].split()[1]
            desc = re.findall('\[(.*?)\]', g[1])  # stuff in [brackets]
            ko = [e for e in desc if 'KO:' in e]
            if len(ko) > 0:
                ko = ko[0][3:].split()
            else:
                ko = ['None']
            ec = [e for e in desc if 'EC:' in e]
            if len(ec) > 0:
                ec = ec[0][3:].split()
            else:
                ec = ['None']
            for i, j in itertools.product(ec, ko):
                mapping.append(pd.Series({'id': gene_id, 'gene': gene_name,
                                          'ec': i, 'ko': j}))
        except:
            print g
    mapping = pd.DataFrame(mapping)
    return mapping


def plot_data_on_pathway(kegg_id, mapping, dist):
    mapping = mapping[mapping.gene.isin(dist.index)]
    order = mapping.gene.map(mapping.groupby('gene').size()).order()
    mapping = mapping.ix[order.index]
    symbol_to_kegg = mapping.set_index('gene').id
    symbol_to_kegg = symbol_to_kegg.groupby(level=0).first()

    dist = pd.Series(dist, name='dist')
    ec = mapping.set_index('gene').join(dist).groupby('ko').median()
    ec = ec.dist.dropna().order()
    gm = pd.concat([mapping.groupby('ko').first().gene, ec], 1)
    gm = gm.set_index('gene').dist.groupby(level=0).first()

    cmap = gm.map(lambda v: rgb2hex(RdBu(1-v)).upper()[1:])
    s = '%0D%'.join(['hsa:{}+%23{}'.format(symbol_to_kegg.ix[i], cmap.ix[i])
                     for i in gm.index if i in symbol_to_kegg])
    s = '{}show_pathway?map={}&multi_query={}'.format(KEGG_PATH, kegg_id, s)
    print s


def parse_entry(e):
    d = dict(e.attrib.iteritems())
    components = [c.attrib['id'] for c in e.findall('component')]
    d['components'] = components
    d = pd.Series(d)
    return d

# remove multiple <?xml version="1.0" encoding="UTF-8"?> and add <T3DB></T3DB> to make a formal format
from io import StringIO
from lxml import etree
import csv
xml = 'drugbank.xml'
ns = {'drugbank':'http://www.drugbank.ca'}
context = etree.iterparse(xml, tag='{http://www.drugbank.ca}drug')

csvfile = open('drugbank.csv', 'w')
fieldnames = ['drugbankid','name', 'cas_registry_number', 'logpexp', 'chemical_formula', 'kingdom', 'direct_parent', 'super_class', 'class', 'sub_class']
writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
writer.writeheader()

for event, elem in context:

    drugbankid = elem.xpath('drugbank:drugbank-id/text()', namespaces=ns)
    try:
        name = elem.xpath('drugbank:name/text()', namespaces=ns)[0]
    except:
        name = 'NA'
    try:
        cas_registry_number = elem.xpath('drugbank:cas-number/text()', namespaces=ns)[0]
    except:
        cas_registry_number = 'NA'
    try:
        kingdom = elem.xpath('drugbank:classification/drugbank:kingdom/text()', namespaces=ns)[0]
    except:
        kingdom = 'NA'
    try:
        direct_parent = elem.xpath('drugbank:classification/drugbank:direct_parent/text()', namespaces=ns)[0]
    except:
        direct_parent = 'NA'
    try:
        super_class = elem.xpath('drugbank:classification/drugbank:super_class/text()', namespaces=ns)[0]
    except:
        super_class = 'NA'
    try:
        classorg = elem.xpath('drugbank:classification/drugbank:class/text()', namespaces=ns)[0]
    except:
        classorg = 'NA'
    try:
        sub_class = elem.xpath('drugbank:classification/drugbank:sub_class/text()', namespaces=ns)[0]
    except:
        sub_class = 'NA'
    try:
        logpexp = elem.xpath('drugbank:experimental-properties/drugbank:property[kind = "Hydrophobicity"]/value/text()', namespaces=ns)[0]
    except:
        logpexp = 'NA'
    try:
        chemical_formula = elem.xpath('drugbank:experimental-properties/drugbank:property[kind = "Molecular Formula"]/value/text()', namespaces=ns)[0]
    except:
        chemical_formula = 'NA'

    writer.writerow({'drugbankid':drugbankid, 'name':name, 'cas_registry_number': cas_registry_number,'logpexp': logpexp, 'chemical_formula': chemical_formula, 'kingdom': kingdom, 'direct_parent': direct_parent, 'super_class': super_class, 'class': classorg, 'sub_class': sub_class})
    # It's safe to call clear() here because no descendants will be
    # accessed
    elem.clear()
# Also eliminate now-empty references from the root node to elem
    for ancestor in elem.xpath('ancestor-or-self::*'):
        while ancestor.getprevious() is not None:
            del ancestor.getparent()[0]
del context

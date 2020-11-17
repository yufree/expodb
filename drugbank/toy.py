from io import StringIO
from lxml import etree
import csv
xml = 'drugbank.xml'
ns = {'drugbank': 'http://www.drugbank.ca'}
context = etree.iterparse(xml, tag='{http://www.drugbank.ca}drug')

csvfile = open('drugbank.csv', 'w')
fieldnames = ['drugbankid','cas_registry_number', 'kingdom']
writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
writer.writeheader()

for event, elem in context:

    drugbankid = elem.xpath('drugbank:drugbank-id/text()', namespaces=ns)

    try:
        cas_registry_number = elem.xpath('drugbank:cas-number', namespaces=ns)[0]
    except:
        cas_registry_number = 'NA'
    try:
        kingdom = elem.xpath('drugbank:classification/drugbank:kingdom/text()', namespaces=ns)[0]
    except:
        kingdom = 'NA'

    writer.writerow({'drugbankid':drugbankid, 'cas_registry_number':cas_registry_number, 'kingdom':kingdom})
    # It's safe to call clear() here because no descendants will be
    # accessed
    elem.clear()
# Also eliminate now-empty references from the root node to elem
    for ancestor in elem.xpath('ancestor-or-self::*'):
        while ancestor.getprevious() is not None:
            del ancestor.getparent()[0]
del context

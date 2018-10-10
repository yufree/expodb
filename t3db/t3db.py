# remove multiple <?xml version="1.0" encoding="UTF-8"?> and add <T3DB></T3DB> to make a formal format
from io import StringIO
from lxml import etree
import csv
xml = 'toxins.xml'

context = etree.iterparse(xml, tag='compound')

csvfile = open('t3db.csv', 'w')
fieldnames = ['omim', 'smpdb', 'keggmap', 'categories', 'route_of_exposure', 'accession', 'origin', 'logpexp', 'carcinogenicity', 'lethaldose', 'classes', 'monisotopic_molecular_weight', 'iupac_name', 'name', 'chemical_formula', 'InChIKey', 'cas_registry_number', 'hmdb', 'drugbank', 'kegg', 'biocyc', 'pubchem', 'chemspider', 'smiles', 'kingdom', 'direct_parent', 'super_class', 'class', 'sub_class', 'molecular_framework']
writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
writer.writeheader()

for event, elem in context:

    accession = elem.xpath('accession/text()')[0]
    try:
        monisotopic_molecular_weight = elem.xpath('monisotopic_moleculate_weight/text()')[0]
    except:
        monisotopic_molecular_weight = 'NA'
    try:
        categories = elem.xpath('categories/category/text()')
    except:
        categories = 'NA'
    try:
        smpdb = elem.xpath('pathways/pathway/smpdb_id/text()')
    except:
        smpdb = 'NA'
    try:
        keggmap = elem.xpath('pathways/pathway/kegg_map_id/text()')
    except:
        keggmap = 'NA'
    try:
        route_of_exposure = elem.xpath('route_of_exposure/text()')[0]
    except:
        route_of_exposure = 'NA'
    try:
        carcinogenicity = elem.xpath('carcinogenicity/text()')[0]
    except:
        carcinogenicity = 'NA'
    try:
        lethaldose = elem.xpath('lethaldose/text()')[0]
    except:
        lethaldose = 'NA'
    try:
        iupac_name = elem.xpath('iupac_name/text()')[0].encode('utf-8')
    except:
        iupac_name = 'NA'
    try:
        name = elem.xpath('common_name/text()')[0].encode('utf-8')
    except:
        name = 'NA'
    try:
        chemical_formula = elem.xpath('chemical_formula/text()')[0]
    except:
        chemical_formula = 'NA'

    try:
        inchikey = elem.xpath('inchikey/text()')[0]
    except:
        inchikey = 'NA'

    try:
        classes = elem.xpath('class/text()')[0]
    except:
        classes = 'NA'

    try:
        cas_registry_number = elem.xpath('cas_registry_number/text()')[0]
    except:
        cas_registry_number = 'NA'
    try:
        hmdb = elem.xpath('hmdb_id/text()')[0]
    except:
        hmdb = 'NA'
    try:
        drugbank = elem.xpath('drugbank_id/text()')[0]
    except:
        drugbank = 'NA'
    try:
        kegg = elem.xpath('kegg_id/text()')[0]
    except:
        kegg = 'NA'
    try:
        omim = elem.xpath('omim_id/text()')[0]
    except:
        omim = 'NA'
    try:
        biocyc = elem.xpath('biocyc_id/text()')[0]
    except:
        biocyc = 'NA'
    try:
        pubchem = elem.xpath('pubchem_compound_id/text()')[0]
    except:
        pubchem = 'NA'
    try:
        chemspider = elem.xpath('chemspider_id/text()')[0]
    except:
        chemspider = 'NA'
    try:
        smiles = elem.xpath('smiles/text()')[0]
    except:
        smiles = 'NA'
    try:
        logpexp = elem.xpath('experimental_properties/property[kind = "logp"]/value/text()')[0]
    except:
        logpexp = 'NA'
    try:
        kingdom = elem.xpath('taxonomy/kingdom/text()')[0]
    except:
        kingdom = 'NA'
    try:
        direct_parent = elem.xpath('taxonomy/direct_parent/text()')[0]
    except:
        direct_parent = 'NA'
    try:
        super_class = elem.xpath('taxonomy/super_class/text()')[0]
    except:
        super_class = 'NA'
    try:
        classorg = elem.xpath('taxonomy/class/text()')[0]
    except:
        classorg = 'NA'
    try:
        sub_class = elem.xpath('taxonomy/sub_class/text()')[0]
    except:
        sub_class = 'NA'
    try:
        molecular_framework = elem.xpath('taxonomy/molecular_framework/text()')[0]
    except:
        molecular_framework = 'NA'
    try:
        origin = elem.xpath('origin/text()')[0]
    except:
        origin = 'NA'

    writer.writerow({'omim': omim, 'smpdb': smpdb, 'keggmap': keggmap, 'route_of_exposure': route_of_exposure, 'categories': categories, 'accession': accession, 'logpexp': logpexp, 'origin': origin, 'carcinogenicity': carcinogenicity, 'lethaldose': lethaldose, 'classes': classes, 'monisotopic_molecular_weight': monisotopic_molecular_weight, 'iupac_name': iupac_name, 'name': name, 'chemical_formula': chemical_formula, 'InChIKey': inchikey, 'cas_registry_number': cas_registry_number, 'hmdb': hmdb, 'drugbank': drugbank, 'kegg': kegg, 'biocyc': biocyc, 'pubchem': pubchem, 'chemspider': chemspider, 'smiles': smiles, 'kingdom': kingdom, 'direct_parent': direct_parent, 'super_class': super_class, 'class': classorg, 'sub_class': sub_class, 'molecular_framework': molecular_framework})
    # It's safe to call clear() here because no descendants will be
    # accessed
    elem.clear()
# Also eliminate now-empty references from the root node to elem
    for ancestor in elem.xpath('ancestor-or-self::*'):
        while ancestor.getprevious() is not None:
            del ancestor.getparent()[0]
del context

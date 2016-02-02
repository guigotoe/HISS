#!/usr/bin/python
# Guillermo Torres MSc,
# ggtorrese@unal.edu.co, guigotoe@gmail.com
# September 2014
# This script parses Uniprot tabular protien list, using the Cross-reference (EMBL) IDs to
# extract the DNA sequences from NCBI data base
# the input file must have the following columns:
#   Entry   Entry_name  Cross-reference_EMBL    Gene_names
# ********************* How to use ***************************
# 
# 
# ************************************************************
###################
import os,re,urllib2,sys
from Bio import Entrez,SeqIO
from Bio.Seq import Seq

## global variables:

infile = 'gdh_Sept2014.txt'#sys.argv[1]#
prefix = 'gdh_'#sys.argv[2]#
try:
    proxyon = 'p'#sys.argv[3]
except IndexError:
    proxyon = None

## opening files to write ##
fa_file = open("%sdna.fa"%prefix,'w')
metadata_file = open("%smetadata.txt"%prefix,'w')
notfound_file = open("%sNotFound.txt"%prefix,'w')


def proxy(http):
    """
    Set the proxy url access. Set for Universidad Nacional de Colombia proxy.
    """
    os.environ["http_proxy"] = "http://cenbio_nal:cenbio_nalx*2010@proxyapp.unal.edu.co:8080"
    #proxy_url = "http://cenbio_nal:cenbio_nalx*2010@proxyapp.unal.edu.co:8080"
    #proxy_support = urllib2.ProxyHandler({'http': proxy_url})
    #opener = urllib2.build_opener(proxy_support)
    #urllib2.install_opener(opener)

def main():
    """
    Retrieve DNA sequences using the "Cross-reference (EMBL)" IDs 
    and "Gene names" of the uniprot tab-format file.
    """
    # If the proxy is needed, uncoment the next line 
    if proxyon:
        proxy(proxyon)
    else:
        pass
    
    # === * === #   
    
    head = "UniprotID\tGI\tEMBL\tGene_Name\tLocus_Tag\tProduct\tTaxonomyID\tMol_Type\tOrganism\tTaxonomy ->\n"
    metadata_file.write(head)
    
    
    
    lc = 0 
    for line in open(infile,'r'):
        line = line.strip("\n")
        if lc != 0 :
            qual = line.split("\t")
            uniprotid = qual[0]
            try: 
                emblid = qual[2].split(";")[0]
            except IndexError:
                notfound_file.write("%s\n"%(uniprotid))
                pass
            gene_name = qual[4].split(";")[0]
            protein_name = qual[5]
            seqinfo = extractinfo(uniprotid,emblid,gene_name,protein_name)
            try:
                write(seqinfo)
            except AttributeError:
                try:
                    emblid = qual[2].split(";")[1]
                    seqinfo = extractinfo(uniprotid,emblid,gene_name,protein_name)
                    write(seqinfo)
                except:
                    notfound_file.write("%s\n"%(uniprotid))         
            except KeyError:
                notfound_file.write("%s\n"%(uniprotid))
                
        lc+=1
def write(seqinfo):
    """
    Write fasta file and metadata file
    """
    # <== Fasta ==> #
    fa_file.write(">%s|%s|%s\n%s\n"%(seqinfo.getgene(),seqinfo.getuniprot(),seqinfo.getembl(),seqinfo.getdna()))
    # <== Metadata ==> #
    source_values = "%s\t%s\t%s"%(seqinfo.gettaxonomyid(),seqinfo.getmol_type(),seqinfo.getorganims())
    cds_values = "%s\t%s\t%s\t%s\t%s"%(seqinfo.getgid(),seqinfo.getembl(),seqinfo.getgene(),seqinfo.getlocus(),seqinfo.getproduct())
    taxonomy = "*".join(seqinfo.gettaxonomy())
    metadata_file.write(str(seqinfo.getuniprot())+"\t"+cds_values+"\t"+source_values+"\t"+taxonomy+"\n")
    
def extractinfo(uniprotid,emblid,gene_name,pname):
    """
    Parse genbank-format file of the DNA sequence.
    """
    Eseq = ExtractedSequence(uniprotid,emblid)
    Entrez.email = 'ggtorrese@unal.edu.co'
    handle = Entrez.efetch(db="nucleotide", id=emblid, rettype="gb", retmode="text")
    try:
        record = SeqIO.read(handle,"gb")
    except ValueError:
        return {}
    Eseq.addtaxonomy(record.annotations['taxonomy'])
    for feature in record.features:
        source = {}
        cds = {}
        if feature.type == 'source':
            info_id = ['mol_type','organism','taxonomy','isolation_source']
            for idx in range(len(info_id)):
                try:
                    if info_id[idx] == 'taxonomy':
                        source[info_id[idx]] = re.findall("\d+",feature.qualifiers['db_xref'][0])[0]
                    else:
                        source[info_id[idx]] = feature.qualifiers[info_id[idx]][0]
                except KeyError:
                    source[info_id[idx]] = ''
            Eseq.addsource(source)
        if feature.type == 'CDS':
            try:
                # the gene is there?? 
                is_there = False
                try:
                    if feature.qualifiers['gene'][0] in gene_name or feature.qualifiers['locus_tag'][0] in gene_name:
                        is_there = True
                except KeyError:
                    try: 
                        if feature.qualifiers['locus_tag'][0] in gene_name:
                            is_there = True
                    except KeyError:
                        if feature.qualifiers['product'][0] in pname:
                            is_there = True
                        elif feature.qualifiers['product'][0].upper() in pname.upper():
                            is_there = True
                        
                if is_there:
                #if feature.qualifiers['gene'][0] in gene_name or feature.qualifiers['locus_tag'][0] in gene_name:
                    info_id = ['gi_id','dna_seq','EC_number','protein_id','translation','gene','locus_tag','product']
                    for idx in range(len(info_id)):
                        try:
                            if info_id[idx] == 'gi_id':
                                cds[info_id[idx]] = re.findall("\d+",feature.qualifiers['db_xref'][0])[0]
                            elif info_id[idx] == 'dna_seq':
                                cds[info_id[idx]] = feature.location.extract(record.seq)
                            elif info_id[idx] == 'EC_number':
                                cds[info_id[idx]] = feature.qualifiers['EC_number'][0]
                            else:
                                cds[info_id[idx]] = feature.qualifiers[info_id[idx]][0]
                        except KeyError:
                            cds[info_id[idx]] = ''
                    Eseq.addcds(cds)
                    break
                Eseq.addcds(cds)
            except KeyError:
                pass
        else:
            pass
    return Eseq
    
class ExtractedSequence():
    """
    Store all relevant information of the gene from GeneBank file format
    """
    def __init__(self,uniprotid,emblid):
        self.embl = emblid
        self.uniprot = uniprotid
    def addsource(self,sourceinfo):
        assert type(sourceinfo) == dict
        self.source = sourceinfo
    def addcds(self,cdsinfo):
        assert type(cdsinfo) == dict
        self.cds = cdsinfo
    def addtaxonomy(self,taxonomy):
        assert type(taxonomy) == list
        self.taxonomy = taxonomy
    def gettaxonomy(self):
        return self.taxonomy
    def getsource(self):
        return self.source
    def getcds(self):
        return self.cds
    def getuniprot(self):
        return self.uniprot
    def getembl(self):
        return self.embl
    def getaa(self):
        self.dna_seq = self.source['translation']
        return self.aa
    def getdna(self):
        self.dna_seq = self.cds['dna_seq']   
        return self.dna_seq
    def getgene(self):
        self.gene = self.cds['gene']   
        return self.gene
    def getlocus(self):
        self.locus_tag = self.cds['locus_tag']   
        return self.locus_tag
    def getgid(self):
        self.gi_id = self.cds['gi_id']   
        return self.gi_id
    def getpid(self):
        self.protein_id = self.cds['protein_id']   
        return self.protein_id
    def getproduct(self):
        self.product = self.cds['product']   
        return self.product
    def getorganims(self):
        self.organism = self.source['organism']   
        return self.organism
    def gettaxonomyid(self):
        self.taxonomyid = self.source['taxonomy']   
        return self.taxonomyid
    def getmol_type(self):
        self.mol_type = self.source['mol_type']   
        return self.mol_type
    def __str__(self):
        return self.id
    
    
    
    
if __name__ == "__main__": main()

#    self.source = {}
#        self.cds = {} 
#        #* Source info ->
#        self.mol_type = ''
#        self.organism = ''
#        self.taxonomyid = ''
#        self.isolation_source = ''        
#        #* CDS info ->
#        self.gene_name = ''
#        self.gi_id = ''
#        self.locus_tag = ''
#        self.product = ''
#        self.dna_seq = ''
#        self.protein_id = ''
#        self.aa_seq = ''
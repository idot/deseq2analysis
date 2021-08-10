#!/usr/bin/env python

## converts an ensemble transcripts fa file to a table 


import gzip,re,sys,collections,os,subprocess
from Bio import SeqIO,SeqUtils

def gtfattr(att, line):
    m = re.search(att+" \"([0-9a-zA-Z-\.]*)\"", line)
    if m:
       return m.group(1)
    else:
       #raise Exception(att+" "+line)
       return None


def attr(att, line):
    m = re.search(att+":(\\w*)", line)
    if m:
       return m.group(1)
    else:
       raise Exception(att+" "+line)
       return None

#the combined downloaded fa files misses many genes
def read_noncoding(inpath):
    seqs = {}
    with open(inpath) as infile:
         for line in infile:
             if line.startswith(">") :
                gene_biotype = attr("gene_biotype", line)
                gene = attr("gene",line)
                name = attr("gene_name", line)
                transcript_biotype = attr("transcript_biotype", line)
                description = attr("description", line)
                tid = line.split()[0][1:]
                seqs[tid] = [gene, str(gene_biotype), str(transcript_biotype), str(description), str(name)]
    return seqs

def addItem(tid, seqs, parsed):
    stored = seqs.get(tid, [])
    stored.append(parsed)
    seqs[tid] = stored    

def read_gtf(inpath):
    seqs = {}
    with open(inpath) as infile:
         for line in infile:
             if not line.startswith("#"):
                gene_biotype = gtfattr("gene_biotype", line)
                gene = gtfattr("gene_id",line)
                name = gtfattr("gene_name", line)
                transcript_biotype = gtfattr("transcript_biotype", line)
                description = gtfattr("description", line)
                tid = gtfattr("transcript_id", line)
                parsed = [gene, str(gene_biotype), str(transcript_biotype), str(description), str(name)]                           
                if tid != None:
                   addItem(tid, seqs, parsed)
    return seqs


#gene:ENSG00000282431.1 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene gene_symbol:TRBD1 description:T cell receptor beta diversity 1 [Source:HGNC Symbol;Acc:HGNC:12158]
def read_fa(inpath):
    seqs = {}
    count = 0
    with open(inpath) as infile:
      for record in SeqIO.parse(infile, "fasta"):
        count = count + 1
        descfield = record.description
        tid = record.id
        #print(descfield)
        gene_biotype = attr("gene_biotype", descfield)
        gene = attr("gene", descfield)
        name = attr("gene_symbol", descfield)
        transcript_biotype = attr("transcript_biotype", descfield)
        #description = attr("description", descfield)
        #description = description.split("[")[0].strip()
        gc = SeqUtils.GC(str(record))
        length = len(record)
        parsed = [gene, str(gene_biotype), str(transcript_biotype), str(name), str(gc), int(length)]
        addItem(tid, seqs, parsed)
    return seqs



def write_transcript_table(gmap, outpath):
    with open(outpath, "w") as outfile:
         outfile.write("transcript\tgene\tgene_biotype\ttranscripts_biotype\tdescription\tgene_name\n")
         for key, values in gmap.items():
             outfile.write(key+"\t"+"\t".join(values)+"\n")

def write_gene_table(genes, outpath):
    with open(outpath, "w") as outfile:
         outfile.write("geneid\tgene_biotype\tgene_name\n")
         for key, values in genes.items():
             outfile.write("\t".join(values)+"\n")

def write_gene_gc_length_table(genes, outpath):
    with open(outpath, "w") as outfile:
         outfile.write("geneid\tgene_biotype\tsome_biotype\tgene_name\tgc\tlength\n")
         for g,values in genes.items(): 
             #print(values)
             outfile.write("\t".join([str(v) for v in values])+"\n")
         
def to_gene_table(trans):
    genes = {}
    for t,values in trans.items():
        g = values[0]
        genes[g] = [values[0],values[1],values[4]]
    return genes

#we take the mean length and mean GC of all transcripts
def to_gene_table_with_gc(trans):
    genes = []
    for t, values in trans.items():
        gene = values[0]
        #print(gene)
        genes[gene] = [values[0],values[1],values[2],values[3],values[4],values[5]]
    return genes

def convertGTF(inpath, outpath):
    seqs = read_gtf(inpath)
    genes = to_gene_table(seqs)
    write_gene_table(genes, outpath)


def convertFasta(inpath, outpath):
    transinfo = read_fa(inpath)
    geneinfo = to_gene_table_with_gc(transinfo)
    write_gene_gc_length_table(geneinfo, outpath) 
   

#convert("Mus_musculus.GRCm38.genes.fa", "ensembl.genes.tab")
#convert("Mus_musculus.GRCm38.79.chr.gtf", "ensembl.genes.tab")

convertFasta("Homo_sapiens.GRCh38.rna.fa", "ensembl.GRCh38.genes.tab")


#what i did before python:"
#cat Mus_musculus.GRCm38.79.chr.gtf | perl -p -i -e 's/.*gene_id "(\w*)".* gene_biotype "(\w*)".*/\1\t\2/' | grep ENSMUSG | sort | uniq > ensembl.gene_biotype.tab 

 


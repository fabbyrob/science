doc = """
Reads in either (1) an annotation file, an (2) file listing genes and expression level, or (3) a summary file. 

(1) adds to the database the position of every gene and exon listed in the annotation
(2) adds expression data to genes already in the database
(3) adds SNP info for exons already in the database
(3) also adds snp info for introns, associates introns with genes already found in the database

"""

import sys
import sqlite3
import getopt
import summary
import annotation

_t = "a"#nnotation, Summary, or Expression

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
    processArgs(3)
    
    con = checkDB()
    
    if _t == "a":
        processAnnotation(con)
    elif _t == "e":
        processExpression(con)
    else:
        processSummary(con)

def processAnnotation(con):
    reader = annotation.Reader(open(sys.argv[2], 'rb'))
    
    genes = []
    for gene in reader:
        gene.sortExons()
        vals = "\'%s\', \'%s\',\'%\s\', %s, %s, 0" % (gene.name, gene.scaf, gene.direction, gene.start, gene.end)
        genes.append(vals)
    con.executemany("INSERT INTO genes (name, scaffold, direction, start, stop, expression) VALUES (?, ?, ?, ?, ?, ?)", genes)

def processExpression(con):
    file = open(sys.argv[2], 'r')
    
    expressions = []
    for line in file:
        line = line.rstrip()
        sline = line.split()
        
        gene = sline[0]
        expression = sline[1]
        
        try:
            ID = con.execute("SELECT gene_id FROM genes WHERE name=\'%s\'" % (gene))
            expressions.append((expression, ID[0]))
        except sqlite3.Error, e:
            sys.stderr.write("Database error: %s\n" % e.args[0])
            
    con.executemany("UPDATE genes SET expression=? WHERE gene_id=?", expressions)
        

def processSummary(con):
    myReader = summary.Reader(open(sys.argv[2], 'rb'))
    
    gene = None
    ffold = 0
    snps = 0
    introns = []
    exons = []
    exon = []
    intron = []
    for site in myReader:
        if isCoding(myReader, site):
            if gene == None:#get the gene ID
                try:
                    res = con.execute("SELECT gene_id FROM genes WHERE scaf=%s AND (start=%s OR stop=%s)"  % (site.CHROM, site.POS, site.POS))
                    gene = res[0]
                except sqlite3.Error, e:
                    sys.stderr.write("Database Error: %s\nCould not get gene ID for gene at %s %s." % (e.args[0], site.CHROM, site.POS))
           
            if site.ALT_NUM > 0:
                snps += 1
                if myReader.TypeToCode['4fold'] in site.Types:
                    ffold += 1

def isCoding(reader, site):
    if reader.TypeToCode['4fold'] in site.Types:
        return True
    if reader.TypeToCode['0fold'] in siteTypes:
        return True
    if reader.TypeToCode['exon'] in siteTypes:
        return True
    
    return False

def checkDB():
    con = None
    try:
        con = sqlite3.connect(sys.argv[1])
        
        #build the tables, if necessary
        cur.execute('CREATE TABLE IF NOT EXISTS genes(gene_id integer primary key, name string, scaffold string, direction string, start integer, stop integer, expression real)')
        cur.execute('CREATE TABLE IF NOT EXISTS exons(exon_id integer primary key, snps integer, ffold_snps integer, zfold_snps integer, gc integer, afs string)')
        cur.execute('CREATE TABLE IF NOT EXISTS gene_exons(exon_id integer, gene_id integer, FOREIGN KEY (exon_id) REFERENCES exons(exon_id), FOREIGN KEY (gene_id) REFERENCES genes(gene_id))')
        cur.execute('CREATE TABLE IF NOT EXISTS introns(intron_id integer primary key, len integer, start integer, stop integer)')
        cur.execute('CREATE TABLE IF NOT EXISTS gene_introns(intron_id integer, gene_id integer, ord integer, FOREIGN KEY (intron_id) REFERENCES introns(intron_id), FOREIGN KEY (gene_id) REFERENCES genes(gene_id))')
        cur.execute('CREATE TABLE IF NOT EXISTS neighbors(intron_id integer, exon_id integer, FOREIGN KEY (exon_id) REFERENCES exons(exon_id), FOREIGN KEY (intron_id) REFERENCES introns(intron_id))')
        cur.execute('CREATE TABLE IF NOT EXISTS intron_bins(intron_id integer, bin integer, gc integer, snps integer, cncs integer, afs string, FOREIGN KEY (intron_id) REFERENCES intron(intron_id))')
    
        return con
    
    except sqlite3.Error, e:
        if con:
            con.rollback()
        
        sys.stderr.write("Database Error: %s\n" % e.args[0])
        sys.exit(1)

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"t:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-t":
            global _t 
            if arg not in "aes":
                sys.stderr.write("Inappropriate input mode specified (%s). Exiting." % arg)
                sys.exit()
            _t = arg
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()

    sys.stderr.write("database: %s infile: %s -t %s\n" % (sys.argv[1], sys.argv[2], _t))
   
use = "python "+__file__.split("/")[-1]+" database myInputFile.txt type"
def usage():
    print (use)
    sys.exit()
    
def details():
    print(doc)
    print (use)
    print("The database will be created if id does not exist. You should only give it a new database file if you are reading in an annotation.")
    print("The input file can either be an annotation, a expression file, or a summary file.")
    print("______________________________________")
    print("option - argument type - default - description")
    print("t - STR - %s - defines the type of file being read in must be one of:\n\ta - annotation file\n\te - expression file\n\ts-summary file")
    
if __name__ == "__main__":   
    __main__()
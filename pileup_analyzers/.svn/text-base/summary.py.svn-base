doc = '''
A parser for summary files. Makes a reader that returns one line at a time.
    >>> import summary
    >>> summary_reader = summary.Reader(open('mySummary.txt','rb'))
    >>> for record in annotation_reader:
    >>>     print record
    
    Summary.Types {STR:STR}
    Symmary.typeToCode {STR:STR}
    Summary.Fields [STR]
    Summary.Genotypes {STR:STR}
    Summary.Samples [STR]
    Summary.TypeNames [STR]
    Summary.TypeCodes [STR]
    
    Site.CHROM    STR
    Site.POS    INT
    Site.REF    STR
    Site.ALT    STR
    Site.REF_NUM    INT
    Site.ALT_NUM    INT
    Site.TOTAL    INT
    Site.Types    [STR]
    Site.DIVERGENCE    INT s:{-1,0,1}
    Site.GENE STR
    Site.DIR STR {0,-,+}
    Site.Genos [STR]
    Site.Genotypes {STR:STR}
    
'''

import sys

class Summary(object):
    def __init__(self):
        self.Fields = []
        self.Types = {}
        self.typeToCode = {}
        self.TypeNames = []
        self.TypeCodes = []
        self.Genotypes = {}
        self.Samples = {}
        self.hadGenes = False
        
    def parse(self, line):
        sline = line.split()
        line = line.rstrip()
        if line.startswith("#site types"):#this is for backwards compatibility with old summaries
            #type line
            newStr = line[line.index("{")+1:-1]
            ls = newStr.split(", ")
            myDict = {}
            for entry in ls:
                key, val = entry.split(": ") 
                key = key[1:-1]
                val = val.replace("'","")
                myDict[key] = val
                self.TypeNames.append(key)
                self.TypeCodes.append(val)
            self.typeToCode = myDict
            self.Types = {v:k for k, v in myDict.items()}
        elif line.startswith("#Genotypes:"):#this is for backwards compatibility with old summaries
            #genotypes line
            newStr = line[line.index(":")+1:]
            ls = newStr.split(",")
            myDict = {}
            for entry in ls:
                sentry = entry.split(" = ") 
                if len(sentry) != 2:
                    sys.stderr.write("Incorrectly formatted Genotype entry: \n\t%s\n"%entry)
                    continue
                key, val = sentry
                key = key.lstrip()
                myDict[key] = val
            self.Genotypes = myDict
        elif line.startswith("#TYPE"):##TYPE 4fold 4 
            sline = line.split()
            self.Types[sline[2]] = sline[1]
            self.typeToCode[sline[1]] = sline[2]
            self.TypeNames.append(sline[1])
            self.TypeCodes.append(sline[2])
        elif line.startswith("#GENOTYPE"):##GENOTYPE;homozygote alternate;A
            sline = line.split(";")
            self.Genotypes[sline[2]] = sline[1]
        elif line.startswith("#CHROM"):
            sline[0] = sline[0][1:]
            self.Fields = sline
            if not "GENE" in line:
                self.Samples = sline[9:]
            else:
                self.hadGenes = True
                self.Samples = sline[11:]#two more fields when gene name is present
        else:
            sys.stderr.write("Non-standard comment line in summary.\n\t%s\n" % line)
            
    def typesAsInt(self):
        newTypes = {}
        for t in self.Types.keys():
            newTypes[int(t)] = self.Types[t]
        return newTypes
    
    def prettyHeader(self):
        #myStr = "#site types: "+str(self.typeToCode)+"\n"
        myStr = ""
        for type, code in self.typeToCode.items():
            myStr += "#TYPE\t%s\t%s\n" % (type, code)
        for code, type in self.Genotypes.items():
            myStr += "#GENOTYPE;%s;%s\n" % (type, code)
        myStr = myStr[0:-1]
        myStr += "\n#"+"\t".join(self.Fields)
        return myStr
    
    def addGenes(self):
        if not "GENE" in self.Fields:
            if len(self.Fields) == 9:
                self.Fields.append("GENE")
                self.Fields.append("DIR")
            else:#individual samples are present, insert gene name and direction before samples
                self.Fields.insert(9, "GENE")
                self.Fields.insert(10, "DIR")
    
    def __str__(self):
        return "Summary\n\tTypes = %s\n\tSamples = %s" % (self.Types, self.Samples)
        
class Site(object):
    def __init__(self, summary, line):
        self.summary = summary
        
        #process line
        sline = line.split()
        
        self.CHROM = sline[0]
        self.POS = int(sline[1])
        self.REF = sline[2]
        self.ALT = sline[3]
        self.REF_NUM = int(sline[4])
        self.ALT_NUM = int(sline[5])
        self.TOTAL = int(sline[6])
        self.Types = sline[7].split(",")
        self.DIVERGENCE = int(sline[8])
        
        if self.summary.hadGenes:
            self.GENE = sline[9]
            self.DIR = sline[10]
            self.Genos = sline[11:]
        else:
            self.GENE = None
            self.DIR = None
            self.Genos = sline[9:]
            
        self.Genotypes = {}
        
        for samp, geno in zip(self.summary.Samples, self.Genos):
            self.Genotypes[samp] = geno
        
    def typeNames(self):
        myTypes = []
        for t in self.Types:
            myTypes.append(self.summary.Types[t])
        return myTypes
    
    def __str__(self):
        myStr = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % (self.CHROM, str(self.POS), self.REF, self.ALT, self.REF_NUM, self.ALT_NUM, self.TOTAL, ",".join([str(t) for t in self.Types]), self.DIVERGENCE)
        if "GENE" in self.summary.Fields:
            myStr += "%s\t%s\t" % (self.GENE, self.DIR)
        if self.Genos:
            myStr += "\t".join(self.Genos)
        return myStr
    
    def prettyStr(self):
        return "Site(CHROM=%s, POS=%s, REF=%s, Types=%s)" % (self.CHROM, self.POS, self.REF, self.typeNames())
    
class Reader(object):
    def __init__(self, filename):
        self.filename = filename
        self.summary = Summary()
        self.stopIter = False
        
        #get the summary info
        self.line = self.filename.next()
        while self.line.startswith("#"):
            self.summary.parse(self.line)
            self.line = self.filename.next()
        
    def __iter__(self):
        return self
    
    def next(self):
        if self.stopIter:
            raise StopIteration
        else:
            while self.line.startswith("#"):
                self.summary.parse(line)
                self.line = self.filename.next()
                
            record = Site(self.summary, self.line)
            self.line = self.filename.next()
            return record
        
    @property
    def Types(self):
        return self.summary.TypeNames
    
    @property
    def Codes(self):
        return self.summary.TypeCodes
    
    @property
    def TypeToCode(self):
        return self.summary.typeToCode
    
    @property
    def CodeToType(self):
        return self.summary.Types
    
    @property
    def Samples(self):
        return self.summary.Samples
    
    @property
    def Genotypes(self):
        return self.summary.Genotypes
    
    @property
    def Header(self):
        return self.summary.prettyHeader()
    
    def addGenes(self):
        self.summary.addGenes()
    
    def addCode(self, code, type):
        self.summary.TypeNames.append(type)
        self.summary.TypeCodes.append(str(code))
        self.summary.Types[str(code)] = type
        self.summary.typeToCode[type] = str(code)
   
if __name__ == "__main__":  
    #TEST
    myReader = Reader(open('io/summary.txt', 'rb'))
    
    print(myReader.Header)
    #
    #print(myReader.summary)
    #print(myReader.summary.Genotypes)
    for site in myReader:
        print(site)
        
    #for site in myReader:
    #    if myReader.summary.getTypeCode('4fold') in site.Types:
    #        maf = min(site.REF_NUM, site.ALT_NUM)
    #        print("%s %s %s" % (site.CHROM, site.POS, maf))
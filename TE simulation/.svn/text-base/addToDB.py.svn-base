doc = """
Takes in a database name where TE simulation data is to be stored and output files from a TE simulation, and inserts the data into the database.

If the database does not exist it is created.

Outputs the ID number of each file after it is inserted, along with the number of data entries that were added. 

If database is committed after each file is inserted. If any database errors occur any changes made since the last file commited will be rolled back.
(i.e. all output for the simulation that caused the error, and the entry in the simulation table will be deleted)

OUTPUT:

>python addToDB.py test.db sim*.txt 

Inserted 11: 50, 10, 0, "R", "[0]", 200, 0.1, 0.001000000, "ALL", 1, 200, 0, 100, 0, 0, 0.000000000, 0.000000000, 0.000000000
Weird line, skipping.
    
Committed file /Users/wiliarj/Desktop/temp/sim.txt with ID number 11 with 17 entries.

Inserted 12: 50, 10, 0, "R", "[0]", 200, 0.1, 0.001000000, "ALL", 1, 200, 0, 100, 0, 0, 0.000000000, 0.000000000, 0.000000000
Weird line, skipping.
    
Committed file /Users/wiliarj/Desktop/temp/sim2.txt with ID number 12 with 17 entries.

Inserted 13: 50, 10, 0, "R", "[0]", 200, 0.1, 0.001000000, "ALL", 1, 200, 0, 100, 0, 0, 0.000000000, 0.000000000, 0.000000000
Weird line, skipping.
    
Committed file /Users/wiliarj/Desktop/temp/sim3.txt with ID number 13 with 17 entries.

Inserted 14: 50, 10, 0, "R", "[0]", 200, 0.1, 0.001000000, "ALL", 1, 200, 0, 100, 0, 0, 0.000000000, 0.000000000, 0.000000000
Weird line, skipping.
    
Committed file /Users/wiliarj/Desktop/temp/sim4.txt with ID number 14 with 17 entries.
"""

import sys
import getopt
import sqlite3

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
    con = None
    
    try:
        con = sqlite3.connect(sys.argv[1])
        files = sys.argv[2:]
        
        cur = con.cursor()
        
        #make the tables if they don't exist
        cur.execute("CREATE TABLE IF NOT EXISTS simulation(ID integer primary key, teNum string, l text, N integer, G integer, u real, x real, calc text, t real, cross integer, s real, r real, K integer, sel text, c real, b integer, e real, z real, v real, I integer)")
        cur.execute("CREATE TABLE IF NOT EXISTS output(SimID integer, generation integer, class integer, fitness real, tenum real, silence real, silencediff real, ne integer, FOREIGN KEY (SimID) REFERENCES simulation(ID), PRIMARY KEY (SimID, generation))")
        cur.execute("CREATE TABLE IF NOT EXISTS finalFreqs(SimID integer, generation integer, freq, FOREIGN KEY (SimID) REFERENCES simulation(ID), PRIMARY KEY (SimID))")
        
        #read in each file, and put in the data
        for f in files:
            infile = open(f)
            
            data = []
            ID = None
            for line in infile:
                #read in all the data as tuples
                line = line.rstrip()
                sline = line.split(",")
                
                if line.startswith("FREQ"):
                    sline = line.split()
                    #get the info for the final frequency
                    cmd = "INSERT INTO finalFreqs (SimID, generation, freq) VALUES (%s, %s, %s)" % (ID, sline[1], sline[2])
                    cur.execute(cmd)
                elif len(line) > 0 and line[0] == "-":
                    #input the simulation
                    sline= line.split()
                    vals = [x for i, x in enumerate(sline) if i%2 == 1]

                    if vals[-1] == "None":
                        vals[-1] = 0

                    valString = "%s, \"%s\", %s, \"%s\", \"%s\", %s, %s, %s, \"%s\", %s, %s, %s, %s, %s, %s, %s, %s, %s, %s" % tuple(vals)
                    cmd = "INSERT INTO simulation (N, tenum, r, sel, l, G, u, x, calc, t, cross, s, K, b, c, e, z, v, I) VALUES ("+valString+")"
                    print(cmd)
                    cur.execute(cmd)
                    
                    #get the ID
                    ID = cur.lastrowid
                    print("\nInserted %s: %s" % (ID, valString))
                else:#data
                    if len(sline) == 7:
                        if ID == None:#oops...
                            sys.stderr.write("No ID in file %s before a data line. Skipping file.\n" % f)
                        data.append([ID]+sline)
                    else:
                        sys.stderr.write("Weird line, skipping.\n\t%s\n" % line)
            print("starting execute manys")
            con.executemany("INSERT INTO output (SimID, generation, class, fitness, tenum, silence, silencediff, ne) values ("+",".join(["?"]*8)+")", data)
            
            
            infile.close()
            con.commit()#commit this file after it is done
            print("Committed file %s with ID number %s with %s entries." % (f, ID, len(data)))
        
        con.close()
    except sqlite3.Error, e:
        
        if con:
            con.rollback()
        
        sys.stderr.write("Database Error: %s\n" % e.args[0])
        sys.exit(1)
        

def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        continue
#        if opt == "-o":
#            global _o 
#            _o = arg
#        else:
#            print ("Unrecognized option: "+opt+"\n")
#            usage()

    sys.stderr.write("infile: %s -arg %s\n" % (sys.argv[1], 1))
   
use = "python "+__file__.split("/")[-1]+""
def usage():
    print (use)
    sys.exit()
    
def details():
    print(doc)
    print (use)
    print("______________________________________")
    #print("OPTION - TYPE - DESCRIPTION (default)")

    
if __name__ == "__main__":   
    __main__()
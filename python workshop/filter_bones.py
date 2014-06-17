#this program should read in a VCF file, the output sites that pass some filters
#filter options should be:
#    minimum depth
#    maximum depth
#    quality
#
#there should also be another option that makes the program print ONLY variable sites


#TODO import any libraries you're going to need

#TODO put your default options here
#these variables typically start with an underscore (_)
#for example _D might be the maximum depth

#since they are declared up here they will be "global" variables
#ask me what that means when i forget to tell you

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
    
    if len(sys.argv) < 3:
        usage()
    
    #TODO read in all the options and process them
    
    #TODO open the input file
        
    #TODO define your regular expression pattern and compile it, if you are going to use a regex.
         
   
    #TODO  read the file
    
    #TODO match this line to your regular expression
    
    #TODO check if the line matches your regular expression (what should we do if it doesnt?)
    
    #TODO output sites that pass the filters


   
#TODO
use = "make this string say how the program should look"
def usage():
    print (use)
    sys.exit()
    
def details():
    print (use)
    #TODO explain each option on one line, state what the default is
    
if __name__ == "__main__":   
    __main__()
doc = """
A script that e-mails a message to a given e-mail.
"""
import sys
import getopt
import smtplib
import os
from email.mime.text import MIMEText

you = ''
message = "you've been notified of something."

def __main__():
    #check aruguments
    if len(sys.argv) == 1:
        details()
        sys.exit()
        
    processArgs(1)
    
    global message
    message += "\n\n---------------------------------------------\nSent from one of the sciencey servers using notifier.py by %s\nDon't try replying to the e-mail it won't get you anywhere. I promise." % (os.environ['USER'])
    
    me = os.environ['USER']+"@"+os.environ['HOSTNAME']
    
    msg = MIMEText(message)
    
    msg['Subject'] = 'Message from %s' % (os.environ['HOSTNAME'])
    msg['From'] = me
    msg['To'] = you
    
    
    s = smtplib.SMTP('localhost')
    s.sendmail(me, [you], msg.as_string()) 
    s.quit()


def processArgs(num):
    try: 
        opts, args = getopt.getopt(sys.argv[num:],"m:r:")
    except getopt.GetoptError:
        usage()
    
    for opt, arg in opts:
        if opt == "-m":
            global message 
            message = arg
        elif opt == "-r":
            global you
            you = arg
        else:
            print ("Unrecognized option: "+opt+"\n")
            usage()

    sys.stderr.write("-m %s -r %s\n" % (message, you))
   
use = "python "+__file__.split("/")[-1]+""
def usage():
    print (use)
    sys.exit()
    
def details():
    print(doc)
    print (use)
    print("______________________________________")
    print("option - argument type - default - description")
    print("m - STR - %s - the body of the message you want to e-mail." % message)
    print("r - STR - %s - the recipient of the e-mail" % you)

    
if __name__ == "__main__":   
    __main__()
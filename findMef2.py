#require for regular expression use
import re
import numpy
def findMEF2():
    #opens and reads sequence.txt file
    with open('sequence.txt') as sequenceFile:
        fasta=sequenceFile.read()
    #closes sequence reader
    sequenceFile.closed
    fasta=fasta.replace(" ", "")
    #removes end of line
    fasta=fasta.replace('\n', "")
    fasta=fasta.replace('\r', "")
    #converts to uppercase
    fasta=fasta.upper()
    #searches for YTAWWWWTAR consensus site
    print fasta
    tfBind=re.compile('[CT]TA[AT][AT][AT][AT]TA[AG]')
    f=tfBind.findall(fasta)
    iterableF=tfBind.finditer(fasta)
    print f
    for match in iterableF:
        mef2Location=match.span()
        print mef2Location
        a=-1*(numpy.subtract((200, -200), match.span()))
        start=mef2Location[0]+a[0]
        end=mef2Location[1]+a[1]
        print start
        print end
    if f:
        print "match found"
        print f
    else:
       print "no match"

if __name__=='__main__':
    findMEF2();

#!/usr/bin/env python
import sys
import os
import re
import csv
import warnings
import numpy
import csv
import scipy
warnings.simplefilter('error')

from optparse import OptionParser
from collections import OrderedDict
#sys.path[0:0]=['/home/sylv/Lagrange-Devel']
sys.path[0:0]=['/Users/sdukesy/Lagrange-Devel/LagrangePlus']
#sys.path[0:0]=['..']
import lagrange



def config():
    parser = OptionParser();
    parser.add_option("-w", "--workdir",    type="string",       dest="workDir",    help="an alternative directory where the model should be run from", default="."      )
    parser.add_option("-m", "--model",      type="string",       dest="modelFilename", help="name of the model file to use",                       )
    parser.add_option("-o", "--outfile",    type="string",       dest="outfile",    help="file to write output to, default = /dev/null",             default="/dev/null")
    parser.add_option("-v", "--verbose",    action="store_true", dest="verbose",    help="creates all the standard output from lagrange",            default=False) 

    options, args = parser.parse_args()

    os.chdir(options.workDir)

    modelFile           = open(options.modelFilename, 'r')
    modelTemplateString = modelFile.read()
    modelFile.close()

    outputFilename = options.outfile
    outputFile     = open(outputFilename, 'w')

    csvWriter = csv.writer(outputFile, delimiter=',', lineterminator='\n')

    return outputFile, csvWriter, modelTemplateString



def main():

    outputFile, csvWriter, modelTemplateString = config()
  
    ############################################
    # This translates the results of the configuration above
    # into something that can be handled in 
    # an automated way.
    #
    ###########################################

    modelTemplate = eval(modelTemplateString)

    matrixTemplate = modelTemplate['area_dispersal']
    initialGuess   = modelTemplate['initial_guess']
    numReps        = modelTemplate['numReps']

    err=0
    paramsDict = OrderedDict()
    paramRegex = re.compile( '<\w+>')
    for t, m in enumerate(matrixTemplate):
        for i, row in enumerate(m):
            for j, elt in enumerate(row):
                if ( type(elt) is str):
                    if ( paramRegex.match(elt) ):
                        if not(elt in paramsDict): paramsDict[elt] = list()
                        paramsDict[elt].append( (t,i,j) )

                    else:
                        print 'badly formed matrix entry',elt
                        err += 1

    if ( err ):
        print 'Something is wrong with the matrix. Fix that and restart'
        return err
                

    print paramsDict.keys()

    csvHeader = ['nLogL','<e>']
    csvHeader[2:] = paramsDict.keys()
    csvWriter.writerow(csvHeader)

    if ( 'runSimple' in modelTemplate ):
        print "Run Simple"
        # From the code below in this file.
        modelTemplate = eval(modelTemplateString)
        initialGuess = modelTemplate['initial_guess']

        modelString = modelTemplateString
        for key, eltList in paramsDict.iteritems():
            for elt in eltList:
                initValue = str( initialGuess[key] )
                modelString = re.compile( key ).sub( initValue, modelString)

        model, tree, data, nodelabels, base_rates = lagrange.input.eval_decmodel(modelString)
        
        #compute likelihoods
        
        # from optimize_FullML
        initialGuessList = [ initialGuess['<e>'] ] #Need to provide an initial guess for e, the extinction rate
        for key in paramsDict.keys():
            initialGuessList.append( initialGuess[ key ] )

        initialGuessList = scipy.log(initialGuessList)

        params = initialGuessList

        nll = lagrange.optimize.likelihood_FullML(params, model, tree, paramsDict)
        
        #print tree
        
        lagrange.output.ascii_tree(outputFile, tree, model, tee=True)
        lagrange.output.ancsplits_fullML(outputFile, tree, model, params, paramsDict)
        
        return 0

    elif ( 'runSingle' in modelTemplate ):
        modelTemplate = eval(modelTemplateString)
        initialGuess = modelTemplate['initial_guess']

        modelString = modelTemplateString
        for key, eltList in paramsDict.iteritems():
            for elt in eltList:
                initValue = str( initialGuess[key] )
                modelString = re.compile( key ).sub( initValue, modelString)

        model, tree, data, nodelabels, base_rates = lagrange.input.eval_decmodel(modelString)

        initialGuessDict = initialGuess

        initialGuessList = [ initialGuessDict['<e>'] ] #Need to provide an initial guess for e, the extinction rate
        for key in paramsDict.keys():
            initialGuessList.append( initialGuessDict[ key ] )

        initialGuessList = scipy.log(initialGuessList)

        nLogL = lagrange.optimize.likelihood_FullML(initialGuessList, model, tree, paramsDict)

        print nLogL

    else:        
        for rep in range(numReps):
            modelTemplate = eval(modelTemplateString)
            initialGuess = modelTemplate['initial_guess']

            modelString = modelTemplateString
            for key, eltList in paramsDict.iteritems():
                for elt in eltList:
                    initValue = str( initialGuess[key] )
                    modelString = re.compile( key ).sub( initValue, modelString)

            model, tree, data, nodelabels, base_rates = lagrange.input.eval_decmodel(modelString)

            r = lagrange.optimize.optimize_FullML(tree, model, paramsDict, initialGuess)

            csvWriter.writerow(r.values())

    outputFile.close()

    return 0

if __name__ == "__main__":
    sys.exit(main())

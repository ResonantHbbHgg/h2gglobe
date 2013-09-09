#!/usr/bin/env python
import os
import sys

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-M","--method",dest="methods",default=[],action="append")
parser.add_option("-d","--datacard",dest="datacard")
parser.add_option("-q","--qqbarcard",dest="qqbarcard")
parser.add_option("-m","--mass",dest="mass",type=float)
parser.add_option("-D","--directory",dest="directory",default="jobs")
parser.add_option("-f","--fitNuis",dest="fitNuis",type="int",default=0)
parser.add_option("-n","--njobs",dest="njobs",type=int)
parser.add_option("-t","--toysperjob",dest="toysperjob",type=int)
parser.add_option("-w","--workspaceName",dest="wsname",default="floatMu.root")
parser.add_option("-s","--skipWorkspace",dest="skipWorkspace",action="store_true",default=False)
parser.add_option("--dryRun",dest="dryRun",action="store_true",default=False)
(options,args)=parser.parse_args()

import ROOT as r
def getFittedSignal(filename):
  filename = filename.replace('.0','')
  tf = r.TFile(filename)
  tree = tf.Get('limit')
  tree.GetEntry(0)
  return tree.limit

def doChannelCompability():
  print 'Running channel compatibility..'
  cardnameSM = options.datacard.replace('.txt','_justSM.txt')
  cardnameGravGG = options.datacard.replace('.txt','_justGravGG.txt')
  cardnameGravQQ = options.datacard.replace('.txt','_justGravQQ.txt')
  f = open('%s/%s/chcomp.sh'%(os.getcwd(),options.directory),'w')
  f.write('#!/bin/bash\n')
  f.write('cd %s/%s\n'%(os.getcwd(),options.directory))
  f.write('rm -f %s.done\n'%f.name)
  f.write('rm -f %s.fail\n'%f.name)
  f.write('touch %s.run\n'%f.name)
  f.write('eval `scramv1 runtime -sh`\n')

  # figure out signal to throw given pre-fit or post fit
  if options.fitNuis==0:
    expectSignalSM=1.
    expectSignalGravGG=1.
    expectSignalGravQQ=1.
  elif options.fitNuis==1:
    print 'Extracting post fit values for channel compatibility toys'
    os.system('combine %s/%s -M MaxLikelihoodFit -m %3.1f --rMin=-5. -n SMRef'%(os.getcwd(),cardnameSM,options.mass))
    os.system('combine %s/%s -M MaxLikelihoodFit -m %3.1f --rMin=-5. -n GravGGRef'%(os.getcwd(),cardnameGravGG,options.mass))
    os.system('combine %s/%s -M MaxLikelihoodFit -m %3.1f --rMin=-5. -n GravQQRef'%(os.getcwd(),cardnameGravQQ,options.mass))
    expectSignalSM = getFittedSignal('higgsCombineSMRef.MaxLikelihoodFit.mH%3.1f.root'%(options.mass))
    expectSignalGravGG = getFittedSignal('higgsCombineGravGGRef.MaxLikelihoodFit.mH%3.1f.root'%(options.mass))
    expectSignalGravQQ = getFittedSignal('higgsCombineGravQQRef.MaxLikelihoodFit.mH%3.1f.root'%(options.mass))
  else:
    sys.exit('Invalid option for fitNuis: %d'%options.fitNuis)

  # throw SM asimov
  f.write('if ( combine %s/%s -M GenerateOnly -m %3.1f -t -1 --expectSignal=%1.3f -n SMAsimov --saveToys -s 0 '%(os.getcwd(),cardnameSM,options.mass,expectSignalSM)) 
  # throw grav_gg asimov
  f.write(' && combine %s/%s -M GenerateOnly -m %3.1f -t -1 --expectSignal=%1.3f -n GravGGAsimov --saveToys -s 0 '%(os.getcwd(),cardnameGravGG,options.mass,expectSignalGravGG))
  # throw grav_qq asimov
  f.write(' && combine %s/%s -M GenerateOnly -m %3.1f -t -1 --expectSignal=%1.3f -n GravQQAsimov --saveToys -s 0 '%(os.getcwd(),cardnameGravQQ,options.mass,expectSignalGravQQ))

  # fit SM asimov 
  f.write(' && combine %s/%s -M ChannelCompatibilityCheck -m %3.1f -t -1 --toysFile higgsCombineSMAsimov.GenerateOnly.mH%3.1f.root --rMin=-5. -g spinCat0 -g spinCat1 -g spinCat2 -g spinCat3 -g spinCat4 -n SM --saveFitResult'%(os.getcwd(),cardnameSM,options.mass,options.mass))
  # fit grav_gg asimov
  f.write(' && combine %s/%s -M ChannelCompatibilityCheck -m %3.1f -t -1 --toysFile higgsCombineGravGGAsimov.GenerateOnly.mH%3.1f.root --rMin=-5. -g spinCat0 -g spinCat1 -g spinCat2 -g spinCat3 -g spinCat4 -n GravGG --saveFitResult'%(os.getcwd(),cardnameSM,options.mass,options.mass))
  # fit grav_qq asimov
  f.write(' && combine %s/%s -M ChannelCompatibilityCheck -m %3.1f -t -1 --toysFile higgsCombineGravQQAsimov.GenerateOnly.mH%3.1f.root --rMin=-5. -g spinCat0 -g spinCat1 -g spinCat2 -g spinCat3 -g spinCat4 -n GravQQ --saveFitResult'%(os.getcwd(),cardnameSM,options.mass,options.mass))
  # fit data
  f.write(' && combine %s/%s -M ChannelCompatibilityCheck -m %3.1f --rMin=-5. -g spinCat0 -g spinCat1 -g spinCat2 -g spinCat3 -g spinCat4 -n Data --saveFitResult '%(os.getcwd(),cardnameSM,options.mass)) 

  f.write(' )\n')
  f.write('\tthen touch %s.done\n'%f.name)
  f.write('\telse touch %s.fail\n'%f.name)
  f.write('fi\n')
  f.write('rm -f %s.run\n'%f.name)
  os.system('chmod +x %s'%f.name)
  f.close()
  if not options.dryRun: 
      os.system('bsub -q 8nh -o %s.log %s'%(f.name,f.name))

def doSeparation():

  if not options.skipWorkspace:
    print 'Using combine tool to make workspace for separation from card', options.datacard
    os.system('text2workspace.py -m %3.1f %s -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs --PO=muFloating -o %s'%(options.mass,options.datacard,options.wsname))

  print 'Making jobs...'

  for j in range(options.njobs):
    f = open('%s/%s/sub%dSep.sh'%(os.getcwd(),options.directory,j),'w')
    f.write('#!/bin/bash\n')
    f.write('cd %s/%s\n'%(os.getcwd(),options.directory))
    f.write('rm -f %s.done\n'%f.name)
    f.write('rm -f %s.fail\n'%f.name)
    f.write('touch %s.run\n'%f.name)
    f.write('eval `scramv1 runtime -sh`\n')
    
    f.write('if ( combine %s -m %3.1f -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 --fitNuis=%d --singlePoint 1 --saveHybridResult -T %d -i 1 --clsAcc 0 --fullBToys -s 0 -n job%dsep );\n'%(options.wsname,options.mass,options.fitNuis,options.toysperjob,j))
    f.write('\tthen touch %s.done\n'%f.name)
    f.write('\telse touch %s.fail\n'%f.name)
    f.write('fi\n')
    f.write('rm -f %s.run\n'%f.name)
    os.system('chmod +x %s'%f.name)

    if not options.dryRun:
      os.system('bsub -q 8nh -o %s.log %s'%(f.name,f.name))
  
def doqqbar():

  options.wsname = options.wsname.replace('.root','_qqbar.root')
  if not options.skipWorkspace:
    print 'Using combine tool to make workspace for qqbar scan from card', options.qqbarcard
    os.system('text2workspace.py -m %3.1f %s -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs --PO=fqqFloating -o %s'%(options.mass,options.qqbarcard,options.wsname))

  print 'Making jobs...'

  # do best fit nll scan
  bf = open('%s/%s/subBFqqbar.sh'%(os.getcwd(),options.directory),'w')
  bf.write('#!/bin/bash\n')
  bf.write('cd %s/%s\n'%(os.getcwd(),options.directory))
  bf.write('rm -f %s.done\n'%bf.name)
  bf.write('rm -f %s.fail\n'%bf.name)
  bf.write('touch %s.run\n'%bf.name)
  bf.write('eval `scramv1 runtime -sh`\n')
  bf.write('if ( combine %s -M MultiDimFit -m %3.1f --redefineSignalPOIs fqq --setPhysicsModelParameters x=1.,MH=%3.1f --setPhysicsModelParameterRanges fqq=0.,1. --algo grid --points 100 )\n'%(options.wsname,options.mass,options.mass)) 
  bf.write('\tthen touch %s.done\n'%bf.name)
  bf.write('\telse touch %s.fail\n'%bf.name)
  bf.write('fi\n')
  bf.write('rm -f %s.run\n'%bf.name)
  os.system('chmod +x %s'%bf.name)
  if not options.dryRun:
    os.system('bsub -q 8nh -o %s.log %s'%(bf.name,bf.name))

  # do likelihood ratio at fqq points
  for j in range(options.njobs):
    f = open('%s/%s/sub%dqqbar.sh'%(os.getcwd(),options.directory,j),'w')
    f.write('#!/bin/bash\n')
    f.write('cd %s/%s\n'%(os.getcwd(),options.directory))
    f.write('rm -f %s.done\n'%f.name)
    f.write('rm -f %s.fail\n'%f.name)
    f.write('touch %s.run\n'%f.name)
    f.write('eval `scramv1 runtime -sh`\n')
    
    f.write('if ( combine %s -m %3.1f -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 --fitNuis=%d --singlePoint 1 --saveHybridResult -T %d -i 1 --clsAcc 0 --fullBToys --setPhysicsModelParameters fqq=0.0 --freezeNuisances fqq -s 0 -n job%dfqq0.00 '%(options.wsname,options.mass,options.fitNuis,options.toysperjob,j))
    f.write('&& combine %s -m %3.1f -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 --fitNuis=%d --singlePoint 1 --saveHybridResult -T %d -i 1 --clsAcc 0 --fullBToys --setPhysicsModelParameters fqq=0.25 --freezeNuisances fqq -s 0 -n job%dfqq0.25 '%(options.wsname,options.mass,options.fitNuis,options.toysperjob,j))
    f.write('&& combine %s -m %3.1f -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 --fitNuis=%d --singlePoint 1 --saveHybridResult -T %d -i 1 --clsAcc 0 --fullBToys --setPhysicsModelParameters fqq=0.50 --freezeNuisances fqq -s 0 -n job%dfqq0.50 '%(options.wsname,options.mass,options.fitNuis,options.toysperjob,j))
    f.write('&& combine %s -m %3.1f -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 --fitNuis=%d --singlePoint 1 --saveHybridResult -T %d -i 1 --clsAcc 0 --fullBToys --setPhysicsModelParameters fqq=0.75 --freezeNuisances fqq -s 0 -n job%dfqq0.75 '%(options.wsname,options.mass,options.fitNuis,options.toysperjob,j))
    f.write('&& combine %s -m %3.1f -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 --fitNuis=%d --singlePoint 1 --saveHybridResult -T %d -i 1 --clsAcc 0 --fullBToys --setPhysicsModelParameters fqq=1.00 --freezeNuisances fqq -s 0 -n job%dfqq1.00 '%(options.wsname,options.mass,options.fitNuis,options.toysperjob,j))
    f.write(' )\n')

    f.write('\tthen touch %s.done\n'%f.name)
    f.write('\telse touch %s.fail\n'%f.name)
    f.write('fi\n')
    f.write('rm -f %s.run\n'%f.name)
    os.system('chmod +x %s'%f.name)

    if not options.dryRun:
      os.system('bsub -q 8nh -o %s.log %s'%(f.name,f.name))

if len(options.methods)==0: options.methods=['ChannelCompatibility','Separation','qqbar']
if options.fitNuis==0 and 'qqbar' in options.methods:
  print 'qqbar scan is ill-defined for pre-fit so i wont run it.'
  options.methods.remove('qqbar')
options.wsname = os.getcwd()+'/'+options.directory+'/'+options.wsname
print options.wsname
os.system('mkdir -p %s'%options.directory)
if 'ChannelCompatibility' in options.methods: doChannelCompability()
if 'Separation' in options.methods: doSeparation()
if 'qqbar' in options.methods: doqqbar()

# fqq scan 
#os.system('text2workspace.py -m %3.1f %s -P HiggsAnalysis.CombinedLimit.HiggsJPC:twoHypothesisHiggs --PO=fqqFloating -o %s'%(options.mass,options.datacard,options.wsname))
#f.write('if ( combine %s -m %3.1f -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 --fitNuis=1 --singlePoint 1 --saveHybridResult -T %d -i 1 --clsAcc 0 --fullBToys --setPhysicsModelParameter fqq=0.025 --freezeNuisances fqq -s 0 -n job%d );\n'%(options.wsname,options.mass,options.toysperjob,j))

# compatibility
# combineCards.py -xc "ALT" datacard_spinanalysis
#combine combine datacard_spinanalysis_justSM.txt -M GenerateOnly -m 125 -t 4000 --expectSignal=1 -n SMToys --saveToys -s -1 
#combine datacard_spinanalysis_justSM.txt -M ChannelCompatibilityCheck -m 125 -t 4000 --toysFile -g spinCat0 -g spinCat1 -g spinCat2 -g spinCat3 -g spinCat4 
#combine datacard_spinanalysis_justGrav.txt -M ChannelCompatibilityCheck -m 125 -t 4000 --toysFile -g spinCat0 -g spinCat1 -g spinCat2 -g spinCat3 -g spinCat4 

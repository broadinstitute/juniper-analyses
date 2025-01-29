import sys
import os
import random
import math
import argparse
#import re
import numpy as np



parser = argparse.ArgumentParser()
parser.add_argument('--inputAl',"-a", help='file containing alignment. It must have a first row with the names of the samples separated by tabs. Then, each of the following rows represents a genomic position (both SNP and non-variable positions), each row has the nucleotide count data for the different samples separated by tabs. The count data must be in the form: 0-10-5-0 (this means that there are 0A, 10C, 5G and 0T at the considere position and sample).')
parser.add_argument('--inputEpi',"-e", help='file containing epi data. It must be a series of rows, one row per host, each row containing the three values separated by tabs: host name, time of the beginning of the exposure of host, time of the end of the exposure of the host.')
parser.add_argument('--inputSample',"-s", help='file containing samples data. It must be a series of rows, one row per sample, each row containing the three values separated by tabs: host name, sample name, time of collection of samples.')
parser.add_argument('--outputF',"-o", help='name of output xml file.')
parser.add_argument('--Ne',"-N", help='number of individuals in the virtual population (default=15, larger values come at higher computational expense).', type=int, default=15)
parser.add_argument('--downSample',"-d", help='Downsample to a certain maximum coverage (by default 100)', type=int, default=100)
#parser.add_argument('--GenomeSize',"-G", help='genome length (default=100000).', type=int, default=100000)
parser.add_argument('--sampleEvery',"-f", help='how often to collect samples from the posterior (default=100, meaning 1 every 100 steps).', type=int, default=100)
parser.add_argument('--mcmc',"-m", help='number of mcmc iterations (default=500000, larger values come at higher computational cost).', type=int, default=5000000)
parser.add_argument('--theta',"-t", help='do you want to estimate theta, the ancestral proportion of polymorphic sites? (default=False).', type=bool, default=False)
parser.add_argument('--seqError',"-E", help='do you want to estimate and account for the sequencing error? (default=False).', type=bool, default=False)
args = parser.parse_args()

virtualN=args.Ne
#coverage=args.SampleSize
Nmcmc=args.mcmc
theta=args.theta
#L=args.GenomeSize
sampleEvery=args.sampleEvery
seqError=args.seqError
downSample=args.downSample

if(args.inputAl=="" or args.outputF=="" or args.inputEpi=="" or args.inputSample==""):
	print("Error, input and output files must be specified with -a, -e, -s and -o options.")
	exit()
    
## Read input files
Hosts=[]
Starts=[]
Ends=[]
Infectors=[]
SampleNames=[]
SamplingDates=[]
SampleSeqs=[]
inpF=open(args.inputEpi)
line=inpF.readline()
while line!="" and line!="\n":
    linelist=line.split()
    Hosts.append(linelist[0])
    Starts.append(float(linelist[1]))
    Ends.append(float(linelist[2]))
    #Infectors.append(int(linelist[3]))
    SampleNames.append([])
    SamplingDates.append([])
    SampleSeqs.append([])
    line=inpF.readline()
inpF.close()
    
inpF=open(args.inputSample)
line=inpF.readline()
while line!="" and line!="\n":
    linelist=line.split()
    SampleNames[Hosts.index(linelist[0])].append(linelist[1])
    SamplingDates[Hosts.index(linelist[0])].append(linelist[2])
    line=inpF.readline()
inpF.close()

inpF=open(args.inputAl)
line=inpF.readline()
Samples=line.split()
line=inpF.readline()
SampleSeqs=[]
for i in range(len(Samples)):
	SampleSeqs.append([])
while line!="" and line!="\n":
	linelist=line.split()
	for i in range(len(Samples)):
		counts=linelist[i].split("-")
		for n in range(4):
			counts[n]=int(counts[n])
		tot=sum(counts)
		if tot>downSample:
			sampledDownArray=np.random.choice(tot, downSample, replace=False)
			newCounts=[0,0,0,0]
			cumulativeSamples=[counts[0],counts[0]+counts[1],counts[0]+counts[1]+counts[2],counts[0]+counts[1]+counts[2]+counts[3]]
			for d in range(downSample):
				if sampledDownArray[d]<cumulativeSamples[0]:
					newCounts[0]+=1
				elif sampledDownArray[d]<cumulativeSamples[1]:
					newCounts[1]+=1
				elif sampledDownArray[d]<cumulativeSamples[2]:
					newCounts[2]+=1
				else:
					newCounts[3]+=1
			linelist[i]=str(newCounts[0])+"-"+str(newCounts[1])+"-"+str(newCounts[2])+"-"+str(newCounts[3])

		SampleSeqs[i].append(linelist[i])

	#Host=int(line.replace(">","").split("_")[0].replace("Case",""))

	line=inpF.readline()
#    while line!="" and line!="\n" and line[0]!=">":
#        linelist=line.split()
#        counts=[int(linelist[0]),int(linelist[1]),int(linelist[2]),int(linelist[3])]
#        sum=0
#        for i in range(4):
#            sum+=counts[i]
#        freqs=[]
#        for i in range(4):
#            freqs.append(float(counts[i])/sum)
#        sample=np.random.multinomial(coverage, freqs, size=1)
#        SampleSeqs[Hosts.index(Host)].append(str(sample[0][0])+"-"+str(sample[0][1])+"-"+str(sample[0][2])+"-"+str(sample[0][3]))
#        line=inpF.readline()
inpF.close()

#NSNPs=len(SampleSeqs[0])
#
#for i in range(len(SampleSeqs)):
#    for k in range((L-NSNPs)/4):
#        SampleSeqs[i].append(str(coverage)+"-0-0-0")
#        SampleSeqs[i].append("0-"+str(coverage)+"-0-0")
#        SampleSeqs[i].append("0-0-"+str(coverage)+"-0")
#        SampleSeqs[i].append("0-0-0-"+str(coverage))

for i in range(len(SampleSeqs)):
    SampleSeqs[i]=", ".join(SampleSeqs[i])

averageLife=0
for i in range(len(Starts)):
    averageLife+=(Ends[i]-Starts[i])
averageLife=averageLife/len(Starts)


            

outF=open(args.outputF+".xml","w")
outF.write("<beast version='2.7' namespace='beast.base.core:beast.base.util:beast.base.math:\n")
outF.write("                     beast.base.inference:beast.base.inference.parameter:beast.base.inference.util:\n")
outF.write("                     beast.base.inference.operator:beast.base.evolution.likelihood:\n")
outF.write("                     beast.base.evolution.sitemodel:beast.base.evolution.alignment:\n")
outF.write("                     beast.base.evolution.substitutionmodel:beast.base.inference.distribution:\n")
outF.write("                     beast.base.evolution.speciation:beast.base.evolution.operator:\n")
outF.write("                     badtrip.evolution.alignment:badtrip.evolution.datatype:badtrip.evolution.sitemodel:\n")
outF.write("                     badtrip.evolution.operators:badtrip.evolution.substitutionmodel:\n")
outF.write("                     badtrip.evolution.likelihood:badtrip.evolution.tree'>\n\n\n")

outF.write("<alignmentPoMo spec=\"AlignmentPoMo\" id=\"alignment\" dataTypePoMo=\"PoModata\">\n\n")
for n in range(len(Samples)):
    outF.write("	<sequence spec=\"SequencePoMo\" id=\"s"+str(Samples[n])+"\" taxon=\"t"+str(Samples[n])+"\" value=\""+SampleSeqs[n]+"\"/>\n\n")
outF.write("</alignmentPoMo>\n\n\n")

outF.write("	<TraitSetPoMo spec=\'TraitSetPoMo\' id=\'traitSet\' direction=\"forward\" units=\"day\"\n        samplesHostsValue=\"")
nS=1
for n in range(len(Hosts)):
	for n2 in range(len(SampleNames[n])):
		outF.write("t"+str(SampleNames[n][n2])+"=H"+Hosts[n])
		if nS<len(Samples):
			outF.write(", ")
			nS+=1
		else:
			outF.write("\"\n        samplesDatesValue=\"")
nS=1
for n in range(len(Hosts)):
	for n2 in range(len(SamplingDates[n])):
		outF.write("t"+SampleNames[n][n2]+"="+str(SamplingDates[n][n2]))
		if nS<len(Samples):
			outF.write(", ")
			nS+=1
		else:
			outF.write("\"\n        HostDatesStartValue=\"")
#outF.write("t"+str(SampleNames[len(Hosts)-1])+"=H"+Hosts[len(Hosts)-1]+"\"\n        samplesDatesValue=\"")
#for n in range(len(Hosts)-1):
#    outF.write("t"+str(SampleNames[n])+"="+str(SamplingDates[n])+", ")
#outF.write("t"+str(SampleNames[len(Hosts)-1])+"="+str(SamplingDates[len(Hosts)-1])+"\"\n        HostDatesStartValue=\"")
for n in range(len(Hosts)-1):
    outF.write("H"+str(Hosts[n])+"="+str(Starts[n])+", ")
outF.write("H"+str(Hosts[len(Hosts)-1])+"="+str(Starts[len(Hosts)-1])+"\"\n        HostDatesEndValue=\"")
for n in range(len(Hosts)-1):
    outF.write("H"+str(Hosts[n])+"="+str(Ends[n])+", ")
outF.write("H"+str(Hosts[len(Hosts)-1])+"="+str(Ends[len(Hosts)-1])+"\">\n      <taxa spec=\'TaxonSet\' alignment=\'@alignment\'/>\n   </TraitSetPoMo>\n\n\n")

#outF.write("     <input spec=\'YuleModel\' id=\"yule\">\n        <parameter spec=\'birthDiffRate\' id=\"birthRate\" value=\"1\"/>\n            <tree spec=\'RandomTree\' id=\'tree\' estimate=\'true\'>\n                <taxa spec=\'Alignment\' idref=\'alignment\'/>\n                <populationModel id=\'ConstantPopulation0\' spec=\'beast.base.evolution.tree.coalescent.ConstantPopulation\'>\n            		<popSize id=\'randomPopSize\' spec=\'parameter.RealParameter\' value=\'1\'/>\n	            </populationModel>\n            </tree>\n     </input>\n\n\n")
    
outF.write("    <treePoMo spec=\'TreePoMo\' id=\'tree\' estimate=\'true\'>\n        <traits spec=\'TraitSetPoMo\' idref=\'traitSet\'/>\n    </treePoMo>\n\n\n")

#outF.write("   <tree estimate=\"true\" id=\"tree\" name=\"stateNode\">\n       <taxonset id=\"TaxonSet.alignment\" spec=\"TaxonSet\">\n           <data idref=\"alignment\" name=\"alignment\"/>\n       </taxonset>\n   </tree>\n\n\n")
   
#outF.write("<input spec=\'YuleModel\' id=\"yule\">\n        <birthDiffRate idref=\"birthRate\"/>\n        <tree idref=\'tree\'/>\n    </input>\n    <parameter id=\"birthRate\" value=\"2.0\" lower=\"0.0\" upper=\"100.0\"/>\n\n\n")

outF.write("<!-- Substitution model (PoMo + HKY), theta (if used) is the frequency of polymorphisms at the root. If root nuc freqs are not estimated, the ones from mutation rates are used. -->\n    <siteModel spec=\"SiteModelPoMo\" id=\"siteModel\">\n     <mutationRate spec=\'RealParameter\' id=\"mutationRate\" value=\"0.001\" lower=\"0.0\"/>\n     <substModel spec=\"PoMoGeneral\" virtualPop=\""+str(virtualN)+"\" estimateRootFreqs=\"false\" useTheta=\""+str(theta)+"\" theta=\"0.01\" id=\"PoMo.substModel\">\n       <fitness spec=\'RealParameter\' id=\"PoMo.fitness\" value=\"1.0 1.0 1.0 1.0\" lower=\"0.8\" upper=\"1.2\"/>\n       <rootNucFreqs spec=\'RealParameter\' id=\"PoMo.rootFreqs\" value=\"0.25 0.25 0.25 0.25\"/>\n       <mutModel spec=\'MutationModel\' id=\"PoMo.mutModel\">\n       		<rateVector spec=\'RealParameter\' id=\"PoMo.mutRates\" value=\"0.00002\" upper=\"0.3\"/>\n       		<nucFreqs spec=\'RealParameter\' id=\"PoMo.nucFreqs\" value=\"0.25 0.25 0.25 0.25\"/>\n       </mutModel>\n     </substModel>\n   </siteModel>\n\n\n")

errStrin=""
if seqError:
	errStrin="     <distribution spec=\'Prior\' x=\"@PoMo.seqError\">\n       <distr spec=\'LogNormalDistributionModel\' M=\"-4.0\" S=\"7.0\"/>\n     </distribution>\n"
outF.write("  <input spec=\'CompoundDistribution\' id=\'parameterPriors\'>\n    <distribution spec=\'Prior\' x=\"@PoMo.mutRates\">\n       <distr spec=\'LogNormalDistributionModel\' M=\"-4.0\" S=\"20.0\"/>\n     </distribution>\n     <distribution spec=\'Prior\' x=\"@mutationRate\">\n       <distr spec=\'LogNormalDistributionModel\' M=\"0.0\" S=\"2.0\"/>\n     </distribution>\n             <distribution spec=\'Prior\' x=\"@PoMo.bottleneck\">\n       <distr spec=\'LogNormalDistributionModel\' M=\"0.0\" S=\"100.0\"/>\n     </distribution>\n"+errStrin+"    </input>\n\n\n")#


#<distribution spec=\'Prior\' x=\"@mutationRate\">\n       <distr spec=\'LogNormalDistributionModel\' M=\"0.0\" S=\"4.0\"/>\n     </distribution>\n\n     <distribution id=\"yule.prior\" idref=\'yule\'/>\n
if seqError:
	valE=0.00000001
else:
	valE=0.0
outF.write("   <input spec=\'TreeLikelihoodPoMoTransmissionSlow\' id=\"treeLikelihood\" treePoMo=\"@tree\">\n     <seqError spec=\'RealParameter\' id=\"PoMo.seqError\" value=\""+str(valE)+"\"/>\n     <bottleneck spec=\'RealParameter\' id=\"PoMo.bottleneck\" value=\"1.0\" upper=\"10000\"/>\n     <data idref=\"alignment\"/>\n     <siteModel idref=\'siteModel\'/>\n   </input>\n\n\n")

#outF.write("   <input spec=\'TreeLikelihoodPoMo\' id=\"treeLikelihood1\">\n     <seqError spec=\'RealParameter\' id=\"PoMo.seqError\" value=\"0.0\"/>\n     <data idref=\"alignment\"/>\n     <tree idref=\"tree\"/>\n     <siteModel idref=\'siteModel\'/>\n   </input>\n\n\n")

outF.write("   <run spec=\"MCMC\" id=\"mcmc\" chainLength=\""+str(Nmcmc)+"\" storeEvery=\"10000\">\n     <state>\n       <stateNode idref=\"tree\"/>\n       <stateNode idref=\"PoMo.bottleneck\"/>\n       <stateNode idref=\"mutationRate\"/>\n       <stateNode idref=\"PoMo.fitness\"/>\n       <stateNode idref=\"PoMo.mutRates\"/>\n       <stateNode idref=\"PoMo.seqError\"/>\n       <stateNode idref=\"PoMo.nucFreqs\"/>\n       <stateNode idref=\"PoMo.rootFreqs\"/>\n     </state>\n\n     <distribution spec=\'CompoundDistribution\' id=\'posterior\'>\n       <distribution idref=\"treeLikelihood\"/>\n       <distribution idref=\"parameterPriors\"/>\n     </distribution>\n\n\n")#       <stateNode idref=\'birthRate\'/>\n

            
outF.write("    <operator spec=\'ScaleOperator\' id=\'mutScaler\'\n 	      parameter=\"@PoMo.mutRates\"\n 	      scaleFactor=\"0.95\" weight=\"5\"/>\n\n    <operator spec=\'ScaleOperator\' id=\'muRateScaler\'\n 	      parameter=\"@mutationRate\"\n 	      scaleFactor=\"0.95\" weight=\"5\"/>\n\n    <operator spec=\'ScaleOperator\' id=\'bottleneckScaler\'\n 	      parameter=\"@PoMo.bottleneck\"\n 	      scaleFactor=\"0.9\" weight=\"5\"/>\n\n    <operator spec=\'UpDownOperator\' id=\'RelativeScaler\' scaleFactor=\"0.95\" weight=\"5\" down=\"@mutationRate\" up=\"@PoMo.mutRates\"/>\n\n        <operator spec=\'UniformTransmission\' weight=\"50\" tree=\"@tree\"/>\n        <operator spec=\'WilsonBaldingTransmission\' weight=\"50\" tree=\"@tree\"/>\n       <operator spec=\'WilsonBaldingNewHeight\' weight=\"50\" tree=\"@tree\"/>\n       <operator spec=\'ParentChildExchange\' weight=\"50\" tree=\"@tree\"/>\n\n        <operator spec=\'SubtreeSlideTransmission\' weight=\"50\" tree=\"@tree\"/>\n\n    <operator spec=\"DeltaExchangeOperator\" id=\"freqExchanger\"\n	       parameter=\"@PoMo.nucFreqs\"\n	       delta=\"0.01\" weight=\"0.5\"/>")#
if seqError:
	outF.write("<operator spec=\'ScaleOperator\' id=\'errRateScaler\'\n 	      parameter=\"@PoMo.seqError\"\n 	      scaleFactor=\"0.95\" weight=\"5\"/>\n\n ")

#	 <operator spec=\"DeltaExchangeOperator\" id=\"freqExchanger\"\n	       parameter=\"@PoMo.nucFreqs\"\n	       delta=\"0.0002\" weight=\"5\"/>\n\n

#outF.write("    <operator spec=\'ScaleOperator\' id=\'mutScaler\'\n 	      parameter=\"@PoMo.mutRates\"\n 	      scaleFactor=\"0.98\" weight=\"5\"/>\n\n	 <operator spec=\"DeltaExchangeOperator\" id=\"freqExchanger\"\n	       parameter=\"@PoMo.nucFreqs\"\n	       delta=\"0.0002\" weight=\"5\"/>\n\n		<operator id=\'treeScaler\' spec=\'ScaleOperator\' scaleFactor=\"0.75\" weight=\"3\" tree=\"@tree\"/>\n        <operator spec=\'Uniform\' weight=\"30\" tree=\"@tree\"/>\n        <operator spec=\'SubtreeSlide\' weight=\"15\" gaussian=\"true\" size=\"0.0077\" tree=\"@tree\"/>\n        <operator id=\'narrow\' spec=\'Exchange\' isNarrow=\'true\' weight=\"15\" tree=\"@tree\"/>\n        <operator id=\'wide\' spec=\'Exchange\' isNarrow=\'false\' weight=\"3\" tree=\"@tree\"/>\n        <operator spec=\'WilsonBalding\' weight=\"1\" tree=\"@tree\"/>\n        <operator spec=\'UpDownOperator\' scaleFactor=\"0.75\" weight=\"3\" down=\"@tree\"/>\n\n\n")#     <operator spec=\'ScaleOperator\' id=\'birthScaler\'\n	       parameter=\"@birthRate\"\n	       scaleFactor=\"0.8\" weight=\"0.3\"/>\n     <operator spec=\"ScaleOperator\" id=\"muRateScaler\"\n	       parameter=\"@mutationRate\"\n	       scaleFactor=\"0.8\" weight=\"0.3\"/>\n\n
    
outF.write("    <logger logEvery=\""+str(sampleEvery)+"\" fileName=\""+args.outputF+".log\">\n       <model idref=\'posterior\'/>\n       <log idref=\"posterior\"/>\n       <log idref=\"parameterPriors\"/>\n       <log idref=\"treeLikelihood\"/>\n       <log idref=\"mutationRate\"/>\n       <log idref=\"PoMo.bottleneck\"/>\n       <log idref=\"PoMo.seqError\"/>\n       <log idref=\"PoMo.mutRates\"/>\n       <log idref=\"PoMo.nucFreqs\"/>\n     </logger>\n\n       <logger logEvery=\""+str(sampleEvery)+"\" fileName=\""+args.outputF+".trees\" mode=\"tree\">\n            <log idref=\"tree\"/>\n        </logger>\n\n     <logger logEvery=\""+str(sampleEvery*10)+"\">\n       <model idref=\'posterior\'/>\n       <log idref=\"posterior\"/>\n       <log idref=\"parameterPriors\"/>\n       <log idref=\"treeLikelihood\"/>\n       <log idref=\"mutationRate\"/>\n       <log idref=\"PoMo.bottleneck\"/>\n       <log idref=\"PoMo.seqError\"/>\n       <log idref=\"PoMo.mutRates\"/>\n       <log idref=\"PoMo.nucFreqs\"/>\n       <ESS spec=\'ESS\' name=\'log\' arg=\"@posterior\"/>\n       <ESS spec=\'ESS\' name=\'log\' arg=\"@treeLikelihood\"/>\n     </logger>\n\n   </run>\n </beast>\n")#       <log idref=\"yule.prior\"/>\n        <log idref=\"birthRate\"/>\n       <log spec=\'TreeHeightLogger\' tree=\'@tree\'/>\n

outF.close()

exit()
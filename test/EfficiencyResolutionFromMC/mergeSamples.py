import sys, os, subprocess
sys.path.append('../')
def main():

	#from histos import samples
	samples = ['ana_datamc_Run2016MuonsOnly_SingleMuonRun2016B-ReReco-v3','ana_datamc_Run2016MuonsOnly_SingleMuonRun2016C-ReReco-v1','ana_datamc_Run2016MuonsOnly_SingleMuonRun2016D-ReReco-v1','ana_datamc_Run2016MuonsOnly_SingleMuonRun2016E-ReReco-v1','ana_datamc_Run2016MuonsOnly_SingleMuonRun2016F-ReReco-v1','ana_datamc_Run2016MuonsOnly_SingleMuonRun2016G-ReReco-v1','ana_datamc_Run2016MuonsOnly_SingleMuonRun2016H-ReReco-v2','ana_datamc_Run2016MuonsOnly_SingleMuonRun2016H-ReReco-v3']
	for sample in samples:
		dirName = sample
		#dirName = sample[1].split("/")[1]
		if os.path.isdir(dirName):
			
			fileList  = [dirName+"/"+f for f in os.listdir(dirName) if os.path.isfile(dirName + "/" + f)]
			print "merging %d files for %s"%(len(fileList),sample)
			command = ["hadd","-f","ana_datamc_%s.root"%sample]			
			command += fileList
			subprocess.call(command,stdout=open(os.devnull, 'wb'))

		else:
			print "no output for sample ", dirName
main()

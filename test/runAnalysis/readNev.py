import json


evFile=open("5TeVZprime_v1.json")
evDict=json.load(evFile)

for i in range(len(evDict)):
	if evDict[i]['file'][0]['nevents']!=100:
		print evDict[i]['file'][0]['nevents']
		print evDict[i]['file'][0]['name']

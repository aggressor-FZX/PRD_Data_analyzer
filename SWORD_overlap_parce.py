import fileinput 
G4Exception = []
Estart = False
CurrentException = None
Outfile = "overLap.txt"

for line in fileinput.input():
	print(line.strip())
	if "G4Exception-START" in line:
		Estart = True
		CurrentException = ""
	elif "G4Exception-END" in line:
		Estart = False
		G4Exception.append(CurrentException)
		with open(Outfile, "a") as fout:
			fout.write(CurrentException + "\n")
		CurrentException = None 

	if Estart:
		CurrentException += line


print (len(G4Exception))

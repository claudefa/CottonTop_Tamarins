import sys,random

for line in sys.stdin:
	line = line.strip()
	fields = line.split()
#['chr1', '188900835', '0/0', '0/0', '0/0', '0/0', '0/0', '.', '0/0', '0/0', '0/0', '0/0', '0/0', '0/0', '0/0', '.', '0/0', '0/0', '.', '0/0', '0/1', '0/0', '.', '0/0', '0/0', '0/0', '0/0', '0/1', '0/0', '0/0',     

	chrom = fields[0]
	position = int(fields[1])
   
	newline = "" 
	for g in fields[2:]:
		if g == "." or g == "./.":
			newline += "9" 
		else:	
			newline += str(sum([int(x) for x in (g.split("/"))]))
		#print (newline)
	print(newline)

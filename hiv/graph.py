from Bio import SeqIO

contigs = {}


for record in SeqIO.parse("graph-contigs.fasta", "fasta"):
    if record.id not in contigs:
    	contigs[record.id] = {}
    	contigs[record.id] = len(str(record.seq))

print(len(contigs), "supercontigs + reference loaded")


data = {}

header = True
for line in open("smithwaterman.csv", "r"):
	if header:
		header = False
	else:
		source, target, dist = line.strip().split(",")
		dist = float(dist)

		minval = -1

		if contigs[source] > contigs[target]:
			minval = contigs[target]
		else:
			minval = contigs[source]
		

		if source not in data:
			data[source] = {}
			data[source][target] = dist / minval

			data[target] = {}
			data[target][source] = dist / minval
		else:
			data[source][target] = dist / minval

			if target not in data:
				data[target] = {}
				data[target][source] = dist / minval

print(len(data), "data loaded")


outputfile = open("edges.csv", "w")
outputfile.write("source,target,lenght\n")

for source in data:
	maxadj = ""
	maxdist = -1

	for adj in data[source]:
		if data[source][adj] >= maxdist:
			if source in contigs and adj in contigs:
				maxadj = adj
				maxdist = data[source][adj]

	if maxdist != -1:
		outputfile.write(source + "," + maxadj + "," + str(maxdist) + "\n")

outputfile.close()



outputfile = open("nodes.csv", "w")
outputfile.write("id,coverage\n")

inputfile = open("coverage.csv", "r")
header = True
maxcov = -1

for line in inputfile:
	if header:
		header = False
	else:
		id_, coverage_ = line.strip().split(",")

		outputfile.write(id_ + "," + coverage_ + "\n")

		if float(coverage_) > maxcov:
			maxcov = float(coverage_)

outputfile.write("reference," + str(maxcov) + "\n")
outputfile.close()
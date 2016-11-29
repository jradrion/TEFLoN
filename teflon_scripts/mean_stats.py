def average(lst):
    return int(sum(lst)/float(len(lst)))

def mean_stats_portal(samples):
    means=[[],[],[],[],[]]
    for sample in samples:
        statsFile = sample[0].replace(".bam", ".stats.txt")
        covFILE = sample[0].replace(".bam", ".cov.txt")
        with open(statsFile, 'r') as fIN:
            for line in fIN:
                if 'average length' in line:
                    means[0].append(int(float(line.split()[-1])))
                if 'insert size average' in line:
                    means[1].append(int(float(line.split()[-1])))
                if 'insert size standard deviation' in line:
                    means[2].append(int(float(line.split()[-1])))
        with open(covFILE, "r") as fIN:
            for line in fIN:
                if line.startswith("Av"):
                    means[3].append(int(float(line.split()[-1])))
                if line.startswith("St"):
                    means[4].append(int(float(line.split()[-1])))
    return [average(means[0]),average(means[1]),average(means[2]),average(means[3]),average(means[4])]


ns = []
with open("SC_1400.txt","r") as f:
    for line in f:
        temp = line.strip().split()
        # if int(temp[1]) > 1:
        ns.append(int(temp[0]))

with open("SC.txt","w+") as f:
    for n in ns:
        f.write(str(n)+"    ")

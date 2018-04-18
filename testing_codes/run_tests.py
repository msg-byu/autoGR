import os

for i in list(range(50,1001,50)):
    with open("struct_enum.in","w+") as f:
        f.write("System \n")
        f.write("bulk \n")
        f.write("1 0 0 \n")
        f.write("0 1 0 \n")
        f.write("0 0 2 \n")
        f.write("2 \n")
        f.write("1 \n")
        f.write("0.0 0.0 0.0  0/1 \n")
 
        f.write("1 "+str(i)+"\n")
        f.write("0.00001 \n")
        f.write("part list of labelings (including incomplete labelings) is used \n")
        f.write("# Concentration ranges \n")
        f.write("5 5 10 \n")
        f.write("5 5 10")

    print("N cells",i)
    os.system("./srHNF.x")

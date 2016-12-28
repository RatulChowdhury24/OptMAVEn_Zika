import os

keyword = 'experiment["Energies"]'

list = os.listdir(".")
for name in list:
    if name[-3:] != ".py":
        continue
    with open(name) as file:
        for line in file:
            if keyword in line:
                print line
                print name
                print "~~~~~~~~~~~"

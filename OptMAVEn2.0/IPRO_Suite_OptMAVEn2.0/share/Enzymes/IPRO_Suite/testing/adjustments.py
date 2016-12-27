import os

list = os.listdir("../modules/")

List = []
for i in list:
    if i[:2] == "IO" and i[-3:] == ".py":
        List.append(i)
List.append("EXPERIMENT.py")     

text = ""
for file in List:
    command = "diff ../modules/" + file + " ~/work/IPRO_Suite/modules/" + file
    os.system(command + " > output.out")
    text += file + "\n\n"
    with open("output.out") as file:
        for line in file:
            text += line

with open("modify.txt", "w") as file:
    file.write(text)

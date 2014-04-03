import os

dataset = []
names = set()
for filename in os.listdir("."):
    if filename.endswith(".csv"):
        with open(filename) as f:
            header= next(f)
            for line in f:
                if "Warning" in line: continue
                if "<br />" in line: continue
                name = line.split(",")[0].replace('"','')
                if name in names: continue # duplicate results available
                names.add(name)
                dataset.append(line)
                
with open("all.csv","w") as w:
    w.write(header)
    for _ in dataset:
        w.write(_)                
                
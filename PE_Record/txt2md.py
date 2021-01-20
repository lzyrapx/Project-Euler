
# read context from lzyrapx_history
history_filename = "lzyrapx_history.txt"

with open(history_filename) as f:
    content = f.readlines()
content = [x.strip() for x in content]

markdown_filename = "README.md"

# clear all content
with open(markdown_filename, 'r+') as f:
    res = f.readlines()
    f.seek(0)
    f.truncate()

# write new context
writer = open(markdown_filename,'w')
writer.write("|Problem ID|Date|time|\n")
writer.write("|------:|---:|------:|\n")
for line in content:
    problem_id = line.split(" ")[0]
    id = problem_id[0:len(problem_id)-1]
    date = line.split(" ")[1] + " " + line.split(" ")[2] + " " + line.split(" ")[3]
    time = line.split(" ")[4]
    writer.write("|" + id + "|" + date + "|" + time + "|\n")
writer.close()


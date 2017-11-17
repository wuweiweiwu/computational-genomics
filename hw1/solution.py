import re

if __name__ == "__main__":
    lines = None
    with open('output.txt') as f:
        lines = f.readlines()
    total = len(lines)
    count = 0
    for line in lines:
        if float(re.split(r"\t", line)[2]) >= 97:
            count += 1
    print "Percent >= 97:", count * 1.0 / total

    dictionary = {}
    for line in lines:
        name = re.split(r"\t", line)[0].split(".")[0]
        if name not in dictionary:
            dictionary[name] = 0
        dictionary[name] += 1

    max_key = max(dictionary, key=dictionary.get)
    print "Most occurance:", max_key, dictionary[max_key]

    total_percent = 0.0
    for line in lines:
        total_percent += float(re.split(r"\t", line)[2])
    print "Average percent:", total_percent / len(lines)

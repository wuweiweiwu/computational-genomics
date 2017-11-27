import collections
import numpy as np
import matplotlib.pyplot as plt

# function to read the fna file. Returns the ids and the dictionary
# of relationships of id -> sequence for fast lookups
def read_fasta_file(filename):
    dictionary = {}
    with open(filename) as f:
        content = f.read().splitlines()
    currId = None
    for index, line in enumerate(content):
        if index % 2 == 0:
            currId = line[1:]
        else:
            dictionary[currId] = line
    return dictionary


def main():
    dictionary = read_fasta_file('seqs.fna')
    y = []
    for i in xrange(7682):
        count = collections.defaultdict(int)
        for key, value in dictionary.items():
            if value[i] == '-':
                continue
            count[value[i]] += 1
        if count:
            total = 0
            for value in count.values():
                total += value
            maxKey = max(count, key=count.get)
            y.append(count[maxKey]/float(total))
        else:
            y.append(0)
    x = np.arange(1, 7683)
    y = np.array(y)
    plt.plot(x, y, "o")
    plt.show()

if __name__ == '__main__':
    main()

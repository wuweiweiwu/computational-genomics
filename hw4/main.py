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


def calculate_variability(dictionary):
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
            y.append(count[maxKey] / float(total))
    return y


def write_variability(variabilities):
    with open('variability.txt', 'w') as f:
        for variability in variabilities:
            f.write(str(variability) + '\n')


def read_variability():
    y = []
    with open('variability.txt') as f:
        content = f.read().splitlines()
    for index, line in enumerate(content):
        y.append(float(line))
    return y


def plot_variability(variability):
    x = np.arange(1, len(variability) + 1)
    y = np.array(variability)
    plt.plot(x, y, linewidth=1.0)
    plt.show()


def main():
    # dictionary = read_fasta_file('seqs.fna')
    # variabilities = calculate_variability(dictionary)
    # write_variability(variabilities)
    y = read_variability()
    plot_variability(y)


if __name__ == '__main__':
    main()

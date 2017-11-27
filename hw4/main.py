# function to read the fna file. Returns the ids and the dictionary
# of relationships of id -> sequence for fast lookups
def read_fasta_file(filename):
    dictionary = {}
    ids = []
    with open(filename) as f:
        content = f.read().splitlines()
    currId = None
    for index, line in enumerate(content):
        if index % 2 == 0:
            currId = line[1:]
            ids.append(currId)
        else:
            dictionary[currId] = line
    return ids, dictionary


def main():
    read_fasta_file('seqs.fna')


if __name__ == '__main__':
    main()

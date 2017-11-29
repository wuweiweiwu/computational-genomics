import collections
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline, LSQUnivariateSpline


'''
function to read the fna file. Returns the ids and the dictionary
of relationships of id -> sequence for fast lookups
'''
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


'''
calculate the variability
count up the # of each base there is (excluding gaps)
find the base that has the most count and divide by the total count for that position
'''
def calculate_variability(dictionary):
    y = []
    for i in xrange(7682):
        count = collections.defaultdict(int)
        # skip the gaps
        for key, value in dictionary.items():
            if value[i] == '-':
                continue
            count[value[i]] += 1
        # if the dictionary is not empty
        if count:
            # get max and divide by total
            maxKey = max(count, key=count.get)
            total = 0
            for value in count.values():
                total += value
            y.append(count[maxKey] / float(total))
    return y


'''
get the variable regions x1, x2 as a list
They are determined by eye.

Since there are a lot of local minima, I picked the locations
where they seemed to be the lowest with respect to the neighboring regions
and use matplotlib plots to estimate where the regions started and ended
'''
def get_variable_regions():
    variable_regions = [
    (0, 20),
    (140, 200),
    (400, 450),
    (538, 600),
    (668, 704),
    (780, 830),
    (953, 1000),
    (1074, 1123),
    (1200, 1241),
    (1379, 1422)
    ]
    return variable_regions


'''
write the variable regions file -> regions.txt
one on each line \t delimited
'''
def write_regions():
    variable_regions = get_variable_regions()
    with open('regions.txt', 'w') as f:
        for x1, x2 in variable_regions:
            f.write(str(x1) + '\t' + str(x2) + '\n')


'''
write the variability file -> variability.txt
one on each line
'''
def write_variability(variabilities):
    with open('variability.txt', 'w') as f:
        for variability in variabilities:
            f.write(str(variability) + '\n')


'''
read the variability file so we dont have to
do all the computations over and over again
'''
def read_variability():
    y = []
    with open('variability.txt') as f:
        content = f.read().splitlines()
    for index, line in enumerate(content):
        y.append(float(line))
    return y


'''
plot the variability over position number
I ignore all the positions that only consists of gaps

I use scipy.interpolate.UnivariateSpline to fit the curve
the commented out code is where i was experimenting with differnt s values

I also plot the variable regions as red horizontal lines to where I think they are
They are determined by eye. Since there are a lot of local minima, I picked the locations
where they seemed to be the lowest with respect to the neighboring regions
'''
def plot_variability(variability):
    # range from 1 - # samples
    x = xs = np.arange(1, len(variability) + 1)
    y = np.array(variability)

    # code to test the different smoothing factors s=
    # spl.set_smoothing_factor(1)
    # xs = np.linspace(0, len(y), 2000)
    # plt.plot(x, y, linewidth=1.0)
    # plt.figure(1)
    # plt.subplot(411)
    # spl = UnivariateSpline(x, y, s=10)
    # plt.plot(xs, y, 'o', ms=.5)
    # plt.plot(xs, spl(xs), 'b', lw=1)
    # plt.subplot(412)
    # spl = UnivariateSpline(x, y, s=20)
    # plt.plot(xs, y, 'o', ms=.5)
    # plt.plot(xs, spl(xs), 'b', lw=1)
    # plt.subplot(413)
    # spl = UnivariateSpline(x, y, s=30)
    # plt.plot(xs, y, 'o', ms=.5)
    # plt.plot(xs, spl(xs), 'b', lw=1)
    # plt.subplot(414)

    # using a spline from scipy to to fit the data points
    # smoothing factor = 40 seems to be the best for my data set
    spl = UnivariateSpline(x, y, s=40)


    # plot the intervals for the variability regions
    # format: [x1, x2], [y1, y2]
    # I used .6 for y because it was around the center of the graph
    variable_regions = get_variable_regions()
    for x1, x2 in variable_regions:
        plt.plot([x1, x2], [.6, .6], 'r-', lw=2)
    plt.plot(xs, y, 'o', ms=.5)
    plt.plot(xs, spl(xs), 'b', lw=1)
    plt.show()


'''
write fasta file (bonus) given a dictionary of id->sequence
'''
def write_fasta_file(dictionary, filename):
    with open(filename, 'w') as f:
        for key, value in dictionary.items():
            f.write('>' + key + '\n')
            f.write(value + '\n')


'''
code to do the bonus stuff
get region 1 and 4
and write the fasta Files
'''
def bonus():
    #select 100 sequences
    dictionary = read_fasta_file('seqs.fna')
    first100 = {key: dictionary[key] for key in dictionary.keys()[:100]}
    first100_no_gaps = collections.defaultdict(str)
    for i in xrange(7682):
        # skip the gaps
        gap_count = 0
        for key, value in first100.items():
            if value[i] == '-':
                gap_count += 1
        # if its all gaps then ignore it
        if gap_count == 100:
            continue
        for key, value in first100.items():
            first100_no_gaps[key] += value[i]
    #get regions 1 and 4
    regions = get_variable_regions()
    reg1 = regions[0]
    reg4 = regions[3]
    reg1_dict = {}
    reg4_dict = {}
    for key, value in first100_no_gaps.items():
        reg1_dict[key] = value[reg1[0]:reg1[1]+1]
    write_fasta_file(reg1_dict, 'reg1.fna')
    for key, value in first100_no_gaps.items():
        reg4_dict[key] = value[reg4[0]:reg4[1]+1]
    write_fasta_file(reg4_dict, 'reg4.fna')
    write_fasta_file(first100_no_gaps, '100-no-gaps.fna')



'''
main function

the first 4 lines calculate the variability
and write the variability file as well as the regions file

the last 2 lines read the values and plot the data
'''
def main():
    # dictionary = read_fasta_file('seqs.fna')
    # variabilities = calculate_variability(dictionary)
    # write_variability(variabilities)
    # write_regions()
    y = read_variability()
    plot_variability(y)



'''
comment out main() and uncomment bonus to run bonus code
'''
if __name__ == '__main__':
    main()
    # bonus()

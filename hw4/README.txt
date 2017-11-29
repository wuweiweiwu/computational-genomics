HW4
wuxx1045

NOTE: I ignore all positions that are consisted of entirely gaps (They are not considered bases)

I use scipy for curve fitting and matplotlib for plotting so there are some dependencies

virtualenv venv
source venv/bin/activate
pip install -r requirements.txt

Do the above to get all the requirements


How to use:

source venv/bin/activate
python main.py

The main function contains all the code to run everything

The first 4 lines calculate the variability and write the variability file as well as the regions file
The last 2 lines read the values and plot the data


So if you just want to plot the data just leave as is and do python main.py
If you also want to recalculate variability and write the files, uncomment the first 4 lines in main()


Files:

main.py - python code to run everything, every function is commented as specified in the assignment. Also discussed how the variable regions are found and marked
requirements.txt - requirements for this code
variability.txt - the file containing the identity values, one on each line
regions.txt - the file containing the start and end position of each variable region (tab delimited)
smoothing.png - plot showing the different smoothing factors
smoothing_w_overlay.png - plot showing the different smoothing factors but also with the data points overlaid
with_regions.png - plot showing the data points with the smoothing curve overlaid and also the variable regions drawn (IMPORTANT)


BONUS files (IMPORTANT):

reg1.fna - FASTA file for region 1
reg4.fna - FASTA file for region 4
100-no-gaps.fna - FASTA file for all 100 sequences with no gaps
BONUS.txt - discussion about the bonus solution
